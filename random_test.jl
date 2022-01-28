include("cplexSolve.jl")
using Random


Random.seed!(1234) # initialisation du germe

#donnees
K = 1000 # nombre de points au sol
N = 50 # nombre de positions de drones
Rcouv=30 #rayon de couverture des drones
Rcom =15 #rayon de communication des drones
L=20 #altitude
P = 10 # nombre de drones à placer

#Les coordonnees d'un point i sont (abs_sol, ord_sol, 0)
abs_sol=[100*rand() for i in 1:K] 
ord_sol=[100*rand() for i in 1:K]

#Les coordonnees d'une position j sont (abs_alt, ord_alt, L)
abs_alt=[100*rand() for j in 1:N]
ord_alt=[100*rand() for j in 1:N]

d = [0 for i in 1:K, j in 1:N]

for i in 1:K, j in 1:N
    d[i,j]=round(sqrt( (abs_sol[i]-abs_alt[j])^2 + (ord_sol[i]-ord_alt[j])^2 + L^2 ))
end

d_com = [0 for j in 1:N, j1 in 1:N]

for j in 1:N, j1 in 1:N
    d_com[j,j1]=round(sqrt( (abs_alt[j]-abs_alt[j1])^2 + (ord_alt[j]-ord_alt[j1])^2 ))
end

graphCouv = falses(N, K)
for s in 1:N 
    for c in 1:K
        if d[c, s] <= Rcouv
            graphCouv[s, c] = true
        end
    end
end

graphCom = falses(N, N)
for i in 1:N-1
    for j in i+1:N
        if d_com[i, j] <= Rcom
            graphCom[i, j] = true
            graphCom[j, i] = true
        end
    end
end

# S = Array{Float64}(undef, 0, 3)
# C = Array{Float64}(undef, 0, 3)

# for i in 1:N
#     println([abs_alt[i] ord_alt[i] L])
#     S = vcat(S, [abs_alt[i] ord_alt[i] L])
# end

# for i in 1:K
#     C = vcat(C, [abs_sol[i] ord_sol[i] 0])
# end

res = "./res/"

isFeasible, isOptimal = cplexSolveLocalisation(res, "randomInstance.txt")
if isFeasible != isOptimal
    println("cplex solve error ! ")
    println("isFeasible : ", isFeasible, " and isOptimal : ", isOptimal)
end

# test with cutting planes 
isFeasible, isOptimal = cplexSolveLocalisation(res, "randomInstance.txt", true)
if isFeasible != isOptimal
    println("cutting planes error ! ")
    println("isFeasible : ", isFeasible, " and isOptimal : ", isOptimal)
end

