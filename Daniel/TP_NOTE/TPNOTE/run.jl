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
    fout = open("res_bnd1", "w")
    for file in ["benders-graphe-hexagone.txt", "benders1.txt", "benders2.txt", "benders3.txt", "benders4.txt"]
        reader(file)
        println(fout)
        println(fout, "| ", file, " | nombre total d'itérations  | temps total (s) | valeur objective |")
        println(fout, "|----------------------------|----------------------------|-----------------|-----------------|")

        println("\n\n cp_ex1")
        (times, ite, obj_v) = cp_ex1()
        println(fout, "|EX1                         |", ite, "                         |", times, "             |", obj_v,
             "              |")

        println("\n\n cp_ex2")
        (times, ite, obj_v) = cp_ex2()
        println(fout, "|EX2                         |", ite, "                         |", times, "             |", obj_v,
        "              |")

        println("\n\n Dijkstra")
        (times, obj_v) = Dijkstra()
        println(fout, "|Dijkstra                    |-                           |", times, "             |", obj_v,
        "              |")
        
    end
    close(fout)

end