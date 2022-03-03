TOL = 0.00001

function cplexSolveMIP()
    include("instance1.txt")

    M = Model(CPLEX.Optimizer) 

    @variable(M, y[1:N], Bin)
    @variable(M, x[1:N] >= 0)

    @constraint(M, sum(x) == D)
    @constraint(M, [i in 1:N], x[i] <= y[i])
    @constraint(M, y[1] >= y[2])
    @constraint(M, y[1] >= y[3])

    @objective(M, Min, sum(F .* y) + sum(C .* x))

    optimize!(M)
    # status of model
    status = termination_status(M)
    isOptimal = status==MOI.OPTIMAL
    println("isOptimal ? ", isOptimal)

    obj_val = objective_value(M)
    println("obj_val : ", obj_val)
    
end





function CP_auto()
    start_time = time()

    M = Model(CPLEX.Optimizer) 
    set_optimizer_attribute(M, "CPXPARAM_Benders_Strategy", 3)

    @variable(M, y[1:N], Bin)
    @variable(M, x[1:N] >= 0)

    @constraint(M, sum(x) == D)
    @constraint(M, [i in 1:N], x[i] <= y[i])
    @constraint(M, y[1] >= y[2])
    @constraint(M, y[1] >= y[3])

    @objective(M, Min, sum(F .* y) + sum(C .* x))

    optimize!(M)
    # status of model
    status = termination_status(M)
    isOptimal = status==MOI.OPTIMAL
    println("isOptimal ? ", isOptimal)

    obj_val = objective_value(M)
    println("obj_val : ", obj_val)

    solved_Time = round(time() - start_time, digits = 2)
    return solved_Time
end