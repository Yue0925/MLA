
function reader(fileName::String)
    datafile = open(fileName)

    global N = parse(Int64, split(readline(datafile), " ")[1])+1
    global M = parse(Int64, split(readline(datafile), " ")[1])

    global matrixAdj = falses(N, N)
    global D = zeros(Int, ((N)))
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
        else
            v = parse(Int, line[1])+1
            d = parse(Int, line[2])
            D[v] = d
        end

    end
end


function test()
    reader("benders1.txt")
end