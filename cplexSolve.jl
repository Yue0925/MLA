# This file contains methods to solve an instance with CPLEX

TOL = 0.00001


function cplexSolveLocalisation(dir::String, fileName::String)
    # charging data (global vars)
    include(dir * fileName)

    # distance euclidienne arrondie
    distCouv(s, c) = round(Int64, sqrt((S[s, 1] - C[c, 1])^2 + (S[s, 2] - C[c,2])^2 + (S[s,3] - C[c,3])^2))
    distCom(i, j) = round(Int64, sqrt((S[i, 1] - S[j, 1])^2 + (S[i, 2] - S[j,2])^2 + (S[i,3] - S[j,3])^2))

    # matrix adjacecy 
    graphCouv = falses(N, K)
    for s in 1:N 
        for c in 1:K
            if distCouv(s, c) <= Rcouv
                graphCouv[s, c] = true
            end
        end
    end

    # matrix adjacecy 
    graphCom = falses(N, N)
    for i in 1:N-1
        for j in i+1:N
            if distCom(i, j) <= Rcom
                graphCom[i, j] = true
                graphCom[j, i] = true
            end
        end
    end


    # modelization
    M = Model(CPLEX.Optimizer)

    # variables
    @variable(M, y[1:N], Bin) # =1, if the site is opened
    @variable(M, z[1:K], Bin) # =1, if the client is covered by a site
    @variable(M, x[1:N, 1:N], Bin) # =1, if an arc between two opened sites belongs to a spanning tree
    @variable(M, q[1:N], Int) # the ordering of each opened site


    # at most P drones
    @constraint(M, sum(y[i] for i in 1:N) <= P)


    # coverage constraint
    for s in 1:N
        for c in 1:K
            if graphCouv[s, c]
                @constraint(M, z[c] <= y[s])
            end
        end
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
    @constraint(M, sum(x[i, j] for i in 1:N for j in 1:N) == (sum(y[i] for i in 1:N) -1))

    # if the sites is not opened, then no arc passes it 
    for v in 1:N
        @constraint(M, sum(x[i, v] + x[v, i] for i in 1:N) <= P * y[v])
    end

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

    # start a chronometer
    start = time()

    # solve the problem
    optimize!(M)

    computationTime = time() - start
    exploredNodes = MOI.get(backend(M), MOI.NodeCount())

    # status of model
    status = termination_status(M)
    isOptimal = status==MOI.OPTIMAL

    # write solution
    outputFile = "./res/" * fileName
    writer = open(outputFile, "w")

    println(writer, "isOptimal ? ", isOptimal)
    println(writer, "objective value : ", objective_value(M))
    println(writer, "time(s) : ", computationTime)
    println(writer, "nodes : ", exploredNodes)

    vertices = Array{Int64, 1}()
    for i in 1:N
        if JuMP.value(y[i]) > TOL
            append!(vertices, i)
        end
    end

    isFeasible = isConnectedComponent(vertices, graphCom)
    println(writer, "Verification is feasible ? ", isFeasible)
    println(writer, "sites = ", vertices)

    close(writer)
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
    # fileName = "grid5_5_L1_Rcouv1_Rcom1_P12.txt"
    # cplexSolveLocalisation(dir, fileName)


    for file in readdir(dir)
        cplexSolveLocalisation(dir, file)
    end
end