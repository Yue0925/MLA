TOL = 0.00001



function sub_model(y_star::Array{Float64,1})
    SM = Model(CPLEX.Optimizer) 

    @variable(SM, v[1:M+N] >= 0)

    @objective(SM, Max, -sum(v[i] * y_star[i] * bnd for i in 1:M) + sum(D[i-M] * v[i] for i in M+1:N+M))

    @constraint(SM, sum(v) == 1) # vertex
    @constraint(SM, v[s+M] == 0) 

    for a in 1:M
        i = Arcs[a][1]
        j = Arcs[a][2]

        @constraint(SM, -v[a] - v[i+M] + v[j+M] <= 0) 
        @constraint(SM, -v[a] + v[i+M] - v[j+M] <= 0)            

    end

    set_silent(SM) # turn off cplex output
    return SM
end


function cp_ex1()
    start_time = time()
    ite = 0

    # --------------
    # step 1
    # --------------

    MP = Model(CPLEX.Optimizer)

    @variable(MP, y[1:M] >=0, Int)

    @objective(MP, Min, sum(y))

    set_silent(MP) # turn off cplex output
    optimize!(MP)
    statusMP = termination_status(MP)
    @info "masterPB status $statusMP"

    while (statusMP == MOI.DUAL_INFEASIBLE) || (statusMP == MOI.OPTIMAL)
        if ite > MAXITE
            optimize!(MP)
            solved_Time = round(time() - start_time, digits = 2)
            return (solved_Time, ite, round(Int, objective_value(MP)))
        end
        ite += 1

        optimize!(MP)
        # status of model
        statusMP = termination_status(MP)
        @info "masterPB status $statusMP"
        println("\n ---------------")
        println(" ite : ", ite)
        #@show value.(MP[:y])

        y_star = value.(MP[:y])
        obj_v = objective_value(MP)
        @show "objective_value(MP) = $obj_v"

        # ----------------
        # step 2
        # ----------------
        SM = sub_model(y_star)
        optimize!(SM)

        # status of model
        statusSM = termination_status(SM)
        @info "subproblem status $statusSM"
        @show objective_value(SM)
        # @show value.(SM[:v])

        # ------------------
        # step 3
        # ------------------
        if termination_status(SM) == MOI.INFEASIBLE_OR_UNBOUNDED
            error("shouldn't apprear here :(")
        elseif TOL < objective_value(SM)
            @info "feasibility cut added ! "
            v_star = value.(SM[:v])
            @constraint(MP, 0 >= -sum(v_star[i] * y[i] * bnd for i in 1:M) + sum(D[i-M] * v_star[i] for i in M+1:N+M))
        else
            println("Ending with subproblem objective = ", objective_value(SM), " total ite : ", ite)
            break
        end

    end

    println("Finally objective = ", objective_value(MP))
    # @show value.(MP[:y])

    solved_Time = round(time() - start_time, digits = 2)
    return (solved_Time, ite, round(Int, objective_value(MP)))

end




function cp_ex2()
    start_time = time()
    ite = 0

    # --------------
    # step 1
    # --------------

    MP = Model(CPLEX.Optimizer)

    @variable(MP, y[1:M] >= 0)

    @objective(MP, Min, sum(y))

    set_silent(MP) # turn off cplex output
    optimize!(MP)
    statusMP = termination_status(MP)
    @info "masterPB status $statusMP"

    function adding_cuts()
        while (statusMP == MOI.DUAL_INFEASIBLE) || (statusMP == MOI.OPTIMAL)
            if ite > MAXITE
                optimize!(MP)
                solved_Time = round(time() - start_time, digits = 2)
                return (solved_Time, ite, round(Int, objective_value(MP)))
            end
    
            ite += 1
    
            optimize!(MP)
            # status of model
            statusMP = termination_status(MP)
            @info "masterPB status $statusMP"
            println("\n ---------------")
            println(" ite : ", ite)
            #@show value.(MP[:y])
    
            y_star = value.(MP[:y])
            obj_v = objective_value(MP)
            @show "objective_value(MP) = $obj_v"
    
            # ----------------
            # step 2
            # ----------------
            SM = sub_model(y_star)
            optimize!(SM)
    
            # status of model
            statusSM = termination_status(SM)
            @info "subproblem status $statusSM"
            @show objective_value(SM)
            # @show value.(SM[:v])
    
            # ------------------
            # step 3
            # ------------------
            if termination_status(SM) == MOI.INFEASIBLE_OR_UNBOUNDED
                error("shouldn't apprear here :(")
            elseif TOL < objective_value(SM)
                @info "feasibility cut added ! "
                v_star = value.(SM[:v])
                @constraint(MP, 0 >= -sum(v_star[i] * y[i] * bnd for i in 1:M) + sum(D[i-M] * v_star[i] for i in M+1:N+M))
            else
                println("Ending with subproblem objective = ", objective_value(SM), " total ite : ", ite)
                break
            end
    
        end
    end

    set_integer.(MP[:y])
    # ite = 0
    # adding_cuts()

    optimize!(MP)

    statusMP = termination_status(MP)
    @info "masterPB status $statusMP"
    obj_v = objective_value(MP)
    @show "objective_value(MP) = $obj_v"

    solved_Time = round(time() - start_time, digits = 2)
    return (solved_Time, ite, round(Int, obj_v))
end



"""
For bande passante equals to 1 only !!
"""
function Dijkstra()
    # start a chronometer
    start = time()

    dist = [Inf for _ in 1:N]
    prec = [ 0 for _ in 1:N]
    Q = [v for v in 1:N]
    dist[s] = 0

    ite = 0

    while size(Q, 1) > 0
        ite += 1

        # the min dist u
        u = reduce((x, y) -> dist[x] ≤ dist[y] ? x : y, Q)
        # remove u from Q
        filter!(e -> e ≠ u, Q)

        # println(" u -> ", u)

        neighbours = filter(v -> matrixAdj[u, v] || matrixAdj[v, u], Q)

        # println("neighbours : ", neighbours)
        
        for v in neighbours
            alt = dist[u] + 1
            if alt < dist[v]
                dist[v] = alt
                prec[v] = u
            end
        end

    end
    solveTime = time() - start

    # @show dist
    # @show prec
    obj_v = sum(dist .* D)
    @show obj_v
    return (round(solveTime, digits = 2), round(Int, obj_v))
end