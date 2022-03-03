include("cutting_planes_beginning.jl")
include("cutting_planes_normal.jl")
include("cutting_planes_auto.jl")
include("cutting_planes_manuel.jl")

function loadData(n::Int )
    global D = div( n , 2 )
    global F = zeros( Int , n)
    global C = zeros( Int , n )
    F[ 1 ] = 7
    C[ 1 ] = 8
    for i in 1:n
        F[i] = rem(F[i]*F[ 1 ], 159 )
        C[i] = rem(C[i]*C[ 1 ], 61 )
    end
end


function test_1()

    include("instance1.txt")

    println("------------")
    println(" DB normal")
    println("------------")
    @time CP_normal()

    println("------------")
    println(" DB manuel")
    println("------------")
    @time CP_manuel()

    println("------------")
    println(" DB auto")
    println("------------")
    @time CP_auto()

    println("------------")
    println(" MIP ")
    println("------------")
    @time cplexSolveMIP()
end


function test_2()
    for n in 5:10
        global N = 10
        loadData(N)

    end
end