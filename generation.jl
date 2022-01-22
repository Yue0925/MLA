# This file contains methods to generate a data set of instances (for the drones localisation problem)

using Printf

"""
An Instance object:
- K: the number of clients
- N: the number of sites
- L: altitude
- S: a list of sites' positions [a o L;]
- C: a list of clients' positions [a o 0;]
- grid: we assume that drones and clients are deployed in a 3D (bord*bord*L) espace
"""
# struct Instance
#     K::Int64
#     N::Int64
#     L::Int64
#     S::Array{Int64,2}
#     C::Array{Int64,2}
#     bord::Int64
#     Rcouv::Int64
#     Rcom::Int64
# end



function generateInstance(L::Int64, bord::Int64, Rcouv::Int64, Rcom::Int64, P::Int64, outputFile::String)

    writer = open(outputFile, "w")
    println(writer, "Rcouv = ", Rcouv)
    println(writer, "Rcom = ", Rcom)
    println(writer, "L = ", L)
    println(writer, "bord = ", bord)
    println(writer, "P = ", P)

    println(writer, "N = ", bord * bord)
    println(writer, "K = ", bord * bord)

    println(writer, "S = [")
    for i in 1:bord
        for j in 1:bord
            println(writer, i, " ", j, " ", L, ";")

            if i==bord && j==bord
                println(writer, i, " ", j, " ", L, "]")
            end
        end
    end
    
    println(writer, "C = [")
    for i in 1:bord
        for j in 1:bord
            println(writer, i, " ", j, " ", 0, ";")

            if i==bord && j==bord
                println(writer, i, " ", j, " ", 0, "]")
            end
        end
    end

    close(writer)
end



function generateAllInstances()
    dir = "./data/"

    for bord in 5:10
        for L in 1:2
            for Rcouv in 1:3
                for Rcom in 1:3
                    P = round(Int64, bord^2/2)
                    fileName = dir * "grid$bord" * "_$bord" * "_L$L" * "_Rcouv$Rcouv" * "_Rcom$Rcom" * "_P$P" * ".txt"
                    generateInstance(L, bord, Rcouv, Rcom, P, fileName)
                end
            end
        end
    end
end