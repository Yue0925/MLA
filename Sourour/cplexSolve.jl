# This file contains methods to solve an instance with CPLEX

using Combinatorics

TOL = 0.00001


# distance euclidienne arrondie
distCouv(s, c) = round(Int64, sqrt((S[s, 1] - C[c, 1])^2 + (S[s, 2] - C[c,2])^2 + (S[s,3] - C[c,3])^2))
distCom(i, j) = round(Int64, sqrt((S[i, 1] - S[j, 1])^2 + (S[i, 2] - S[j,2])^2))


# function subsetGeneration(size::Int64, V::Array{Int64,1})
    
# end


function cplexSolveLocalisation(dir::String, fileName::String, BranchCut=false)
    countCuts = 0
    # modelization
    M = Model(CPLEX.Optimizer)
    set_optimizer_attribute(M, "CPXPARAM_TimeLimit", 500) # seconds


    if BranchCut
        MOI.set(M, MOI.NumberOfThreads(), 1) 
    end

    # variables
    @variable(M, y[1:N], Bin) # =1, if the site is opened
    @variable(M, z[1:K], Bin) # =1, if the client is covered by a site
    @variable(M, x[1:N, 1:N], Bin) # =1, if an arc between two opened sites belongs to a spanning tree
    @variable(M, q[1:N], Int) # the ordering of each opened site


    # at most P drones
    @constraint(M, sum(y[i] for i in 1:N) <= P)


    # coverage constraint
    for c in 1:K
        @constraint(M, z[c] <= sum(y[s] for s in 1:N if graphCouv[s, c]))
    end


    # -----------------------------------------
    # prefix values (in communication network)
    # -----------------------------------------
    @constraint(M, x[N, N] == 0) # no bucles

    for i in 1:N-1
        @constraint(M, x[i, i] == 0)

        for j in i+1:N
            # if arc x[ij]=1, then arc x[ji]=0
            @constraint(M, x[i, j] + x[j, i] <= 1) 

            # if i, j are not adjacent, then x[ij] = 0
            if !graphCom[i, j]
                @constraint(M, x[i, j] == 0)
                @constraint(M, x[j, i] == 0)
            end
        end
    end


    # --------------------------------
    # a tree covering all opened sites
    # --------------------------------
    
    # the number of arcs equals to the total sites -1
    @constraint(M, sum(x[i, j] for i in 1:N for j in 1:N if graphCom[i, j]) == (sum(y[i] for i in 1:N) -1))

    # if the sites is not opened, then no arc passes it
    # on the contrary, if a site is opened, then at least one arc passes 
    for v in 1:N
        @constraint(M, sum(x[i, v] + x[v, i] for i in 1:N) <= P * y[v])
        @constraint(M, sum(x[i, v] + x[v, i] for i in 1:N) >= y[v])
        # each vertex has at most one incoming arc
        @constraint(M, sum(x[i, v] for i in 1:N if graphCom[i, v]) <= y[v])
    end

    # in a spanning tree, exactly one vertex that has no incoming arc
    @constraint(M, sum(y[v] for v in 1:N) - 1 == sum(x[i, v] for v in 1:N for i in 1:N if graphCom[i, v]))

    # the q[v]=0, if v is not opened, otherwise, 1 <= q[v] <= P
    # more precisely, any q[v] cannot exceed to the total number of opened sites
    for v in 1:N
        @constraint(M, y[v] <= q[v])
        @constraint(M, q[v] <= P * y[v])
        @constraint(M, q[v] <= sum(y[i] for i in 1:N))
    end
    
    # the elimination sub-cycle constraint MTZ
    for u in 1:N
        for v in 1:N
            if graphCom[u, v]
                @constraint(M, q[v] >= q[u] + 1 - P * (1-x[u, v]))
            end
        end
    end

    # objective is to maximize the number of covered clients
    @objective(M, Max, sum(z[i] for i in 1:K))


    # ------------------------------------------------------------------
    # callback function for the cutting planes algorithm
    # vaild inegality constraint added when the actual sol is violated
    # ------------------------------------------------------------------
    function callback_cuttingPlanes(cb_data::CPLEX.CallbackContext)
        println("callback")

        # get current variables
        x_star = zeros(N, N)
        y_star = zeros((N))
        for i in 1:N
            y_star[i] = callback_value(cb_data, y[i])
            for j in 1:N
                x_star[i, j] = callback_value(cb_data, x[i, j])
            end
        end

        # oppened sites
        V = filter(id-> y_star[id]>TOL, [i for i in 1:N])
        subsets = reverse!(collect(powerset(V)))
        for subset in subsets
            if size(subset, 1) < 2
                break
            end
            # when the sub-tour constraint is violated
            if sum(x_star[i, j] for i in subset for j in subset) > sum(y_star[i] for i in subset)-1
                println("inegality added !")
                constr = @build_constraint(sum(x[i, j] for i in subset for j in subset) <= sum(y[i] for i in subset)-1)
                MOI.submit(M, MOI.UserCut(cb_data), constr)
                countCuts += 1
                break
            end
        end

    end

    if BranchCut
        MOI.set(M, MOI.UserCutCallback(), callback_cuttingPlanes)
    end

    # start a chronometer
    start = time()

    # solve the problem
    optimize!(M)

    computationTime = time() - start
    exploredNodes = MOI.get(backend(M), MOI.NodeCount())

    # status of model
    status = termination_status(M)
    isOptimal = status==MOI.OPTIMAL
    isFeasible = false

    # write solution
    outputFile = dir * fileName
    writer = open(outputFile, "a")

    if BranchCut
        println(writer, " --------------------------------------")
        println(writer, "Cutting Planes used with ", countCuts, " cuts appiled ! ")
    end
    println(writer, "isOptimal ? ", isOptimal)
    println(writer, "time(s) : ", computationTime)
    println(writer, "nodes : ", exploredNodes)


    if isOptimal
        println(writer, "objective value : ", objective_value(M))
        vertices = Array{Int64, 1}()
        orders = zeros(N)
        for i in 1:N
            if JuMP.value(y[i]) > TOL
                append!(vertices, i)
            end
            if JuMP.value(q[i]) > TOL
                orders[i] = round(Int64, JuMP.value(q[i]))
            end
    
        end
    
        isFeasible = isConnectedComponent(vertices, graphCom)
        println(writer, "Verification is feasible ? ", isFeasible)
        println(writer, "sites = ", vertices)

    end

    close(writer)
    return isFeasible, isOptimal
end


function isConnectedComponent(vertices::Array{Int64, 1}, graphCom::BitArray{2})
    isVisited = Dict(i=>false for i in vertices)
    todo = []

    # pick up a source
    s = vertices[1]
    append!(todo, s)

    while size(todo, 1) >0
        v = pop!(todo)
        if ! isVisited[v] 
            isVisited[v] = true
            for u in vertices
                if graphCom[v, u] && ! isVisited[u]
                    append!(todo, u)
                end
            end
        end
    end

    for v in vertices
        if ! isVisited[v]
            return false
        end
    end
    return true
end

function test()
    dir = "./data/"
    # fileName = "grid9_9_L2_Rcouv3_Rcom3_P40.txt"
    # isFeasible, isOptimal = cplexSolveLocalisation(dir, fileName)
    # cplexSolveLocalisation(dir, fileName, true)


    for file in readdir(dir)
        # charging data (global vars)
        include(dir * file)

        # matrix adjacecy 
        global graphCouv = falses(N, K)
        for s in 1:N 
            for c in 1:K
                if distCouv(s, c) <= Rcouv
                    graphCouv[s, c] = true
                end
            end
        end

        # matrix adjacecy 
        global graphCom = falses(N, N)
        for i in 1:N-1
            for j in i+1:N
                if distCom(i, j) <= Rcom
                    graphCom[i, j] = true
                    graphCom[j, i] = true
                end
            end
        end

        res = "./res/grid/"

        isFeasible, isOptimal = cplexSolveLocalisation(res, file)
        if isFeasible != isOptimal
            println("cplex solve error ! ")
            println("isFeasible : ", isFeasible, " and isOptimal : ", isOptimal)
            println(file)
            break
        end

        # test with cutting planes 
        isFeasible, isOptimal = cplexSolveLocalisation(res, file, true)
        if isFeasible != isOptimal
            println("cutting planes error ! ")
            println("isFeasible : ", isFeasible, " and isOptimal : ", isOptimal)
            println(file)
            break
        end
    end
end