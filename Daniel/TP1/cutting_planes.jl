
TOL = 0.00001


function modelMaster()
    global MP = Model(CPLEX.Optimizer) 
    @variable(MP, y[1:N], Bin)
    @variable(MP, w >=0)

    @objective(MP, Min, sum(F[i] * y[i] for i in 1:N) + w)

    @constraint(MP, y[1] >= y[3])
    @constraint(MP, y[1] >= y[2])

    # TODO : feasible cut
    @constraint(MP, D - sum(y) <= 0)
end

function modelSub(y::Array{Float64,1})

    SM = Model(CPLEX.Optimizer) 

    @variable(SM, v[1:N] >= 0)
    @variable(SM, b)

    @objective(SM, Max, -sum(v[i] * y[i] for i in 1:N) + D * b)

    @constraint(SM, [i in 1:N], b - v[i] <= C[i])

    set_silent(SM) # turn off cplex output

    return SM
end


function cutting_planes()
    include("instance1.txt")

    ite = 1
    modelMaster()


    while true
        # --------------
        # step 1
        # --------------

        set_silent(MP) # turn off cplex output
        optimize!(MP)
        # println(solution_summary(MP))
        # println(MP)

        # status of model
        statusMP = termination_status(MP)
        isOptimalMP = statusMP==MOI.OPTIMAL
        println("masterPB isOptimal? ", isOptimalMP)

        y_star = value.(MP[:y])
        w_star = JuMP.value(MP[:w])
        y = MP[:y]
        w = MP[:w]

        # ----------------
        # step 2
        # ----------------
        SM = modelSub(y_star)
        optimize!(SM)

        # status of model
        statusSM = termination_status(SM)
        isOptimalSM = statusSM == MOI.OPTIMAL
        println("subproblem isOptimal? ", isOptimalSM)


        b_star = JuMP.value(SM[:b])
        v_star = JuMP.value.(SM[:v])


        println("\n\n ite : ", ite)
        println("objective Master = ", objective_value(MP))
        println("w_star = ", w_star)
        println("Subproblem objective = ", objective_value(SM))
        println("y_star = ", y_star)

        # ---------
        # step 3
        # ---------
        if  objective_value(SM) - w_star > TOL
            println("b_star : ", b_star)
            println("v_star : ", v_star)
            @constraint(MP, w >= D * b_star - sum(y[i] * v_star[i] for i in 1:N))
            println("violated")
        else
            break
        end

        ite += 1
    end


end
