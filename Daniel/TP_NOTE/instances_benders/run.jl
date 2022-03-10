include("ex.jl")

function reader(fileName::String)
    datafile = open(fileName)

    global N = parse(Int64, split(readline(datafile), " ")[1])
    global M = parse(Int64, split(readline(datafile), " ")[1])

    global matrixAdj = falses(N, N)
    global D = zeros(Int, ((N)))
    global Arcs = [[] for _ in 1:M]
    global bnd = 1
    global s = 1
    data = readlines(datafile)
    close(datafile)

    counter = 0
    for eachLine in data
        counter += 1
        line = split(eachLine, " ")

        if counter <= M
            u = parse(Int, line[1])+1
            v = parse(Int, line[2])+1
            matrixAdj[u, v] = true
            Arcs[counter] = [u, v]
        elseif counter<= M+N
            v = parse(Int, line[1])+1
            d = parse(Int, line[2])
            D[v] = d
        end

    end
end


function test()
    
    reader("benders2.txt")
    println("\n\n cp_ex1")
    cp_ex1()
    println("\n\n cp_ex2")
    cp_ex2()
    println("\n\n Dijkstra")
    Dijkstra()
    # benders1.txt
    # benders2.txt
    # benders3.txt
    # benders-graphe-hexagone.txt
end