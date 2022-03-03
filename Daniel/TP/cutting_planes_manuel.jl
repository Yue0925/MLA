
TOL = 0.00001


function modelMaster()
    global MP = Model(CPLEX.Optimizer) 
    @variable(MP, y[1:N], Bin)
    @variable(MP, w >=0)

    @objective(MP, Min, sum(F .* y) + w)

    @constraint(MP, y[1] >= y[3])
    @constraint(MP, y[1] >= y[2])

end

function modelSubManuel(y::Array{Float64,1})

    SM = Model(CPLEX.Optimizer) 

    @variable(SM, v[1:N] >= 0)
    @variable(SM, b)

    @objective(SM, Max, -sum(v .* y) + D * b)

    @constraint(SM, [i in 1:N], b - v[i] <= C[i])

    set_silent(SM) # turn off cplex output

    return SM
end


function CP_manuel()
    start_time = time()
    
    ite = 1
    modelMaster()

    # --------------
    # step 1
    # --------------

    set_silent(MP) # turn off cplex output
    optimize!(MP)
    statusMP = termination_status(MP)
    @info "masterPB status $statusMP"


    while (statusMP == MOI.DUAL_INFEASIBLE) || (statusMP == MOI.OPTIMAL)
        optimize!(MP)
        # status of model
        statusMP = termination_status(MP)
        @info "masterPB status $statusMP"

        y_star = value.(MP[:y])
        w_star = JuMP.value(MP[:w])
        y = MP[:y]
        w = MP[:w]

        println("\n\n ite : ", ite)
        println("objective Master = ", objective_value(MP))
        println("y_star = ", y_star)
        println("w_star = ", w_star)


        # ----------------
        # step 2
        # ----------------
        obj_sub = 0

        if sum(y_star) == D
            @info " ∑y == D"
            obj_sub = sum(y_star .* C)
            b_star = maximum(C)
            v_star = [b_star - C[i] for i in 1:N]

        elseif sum(y_star) > D
            @info " ∑y > D"
            SM = modelSubManuel(y_star)
            optimize!(SM)
            # status of model
            statusSM = termination_status(SM)
            @info "subproblem status $statusSM"


            obj_sub = objective_value(SM)
            b_star = JuMP.value(SM[:b])
            v_star = JuMP.value.(SM[:v])
        else
            @info " ∑y < D"
            @info "Feasibility cut found"
            @constraint(MP, D - sum(y) <= 0)
            ite +=1
            continue
        end

        # ------------------
        # step 3
        # ------------------

        if  obj_sub - w_star > TOL
            @info "Optimality cut found"

            println("b_star : ", b_star)
            println("v_star : ", v_star)
            @constraint(MP, w >= D * b_star - sum(y .* v_star))

        else
            println("Subproblem objective = ", obj_sub)
            break
        end

        ite += 1
    end

    solved_Time = round(time() - start_time, digits = 2)
    return solved_Time
end
