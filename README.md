# MLA TP noté
 Problème localisation avec décomposition benders.


# Exécution

Se déplacer dans le répertoire "TP_NOTE/TPNOTE", et taper ci-dessous dans le terminal ebvironnement Julia

```julia
using CPLEX 
using JuMP

include("run.jl")
test() 
```


# Résultat

## Exercice 1 & 2



| benders-graphe-hexagone.txt | nombre total d'itérations  | temps total (s) | valeur objective |
|----------------------------|----------------------------|-----------------|-----------------|
|EX1                         |9                           |54.51             |21              |
|EX2                         |9                           |0.35             |21              |
|Dijkstra                    |-                           |0.31705498695373535             |21.0              |

| benders1.txt | nombre total d'itérations  | temps total (s) | valeur objective |
|----------------------------|----------------------------|-----------------|-----------------|
|EX1                         |54                          |2.05             |688              |
|EX2                         |78                          |1.17             |688              |
|Dijkstra                     |-                           |0.0003311634063720703             |688.0              |

| benders2.txt | nombre total d'itérations  | temps total (s) | valeur objective |
|----------------------------|----------------------------|-----------------|-----------------|
|EX1                         |112                          |6.78             |390              |
|EX2                         |84                          |1.64             |390              |
|Dijkstra                     |-                           |0.00041604042053222656             |390.0              |

| benders3.txt | nombre total d'itérations  | temps total (s) | valeur objective |
|----------------------------|----------------------------|-----------------|-----------------|
|EX1                         |128                          |11.46             |579              |
|EX2                         |117                          |2.91             |579              |
|Dijkstra                     |-                           |0.0006771087646484375             |579.0              |

| benders4.txt | nombre total d'itérations  | temps total (s) | valeur objective |
|----------------------------|----------------------------|-----------------|-----------------|
|EX1                         |119                          |5.12             |93              |
|EX2                         |70                          |1.4             |93              |
|Dijkstra                     |-                           |0.00026297569274902344             |93.0              |
