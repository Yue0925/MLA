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

## Avec **bande passant de 1**, on obtient les résultats suivants : 


| benders-graphe-hexagone.txt | nombre total d'itérations  | temps total (s) | valeur objective |
|----------------------------|----------------------------|-----------------|-----------------|
|EX1                         |9                         |0.24             |21              |
|EX2                         |9                         |0.29             |21              |
|Dijkstra                    |-                           |0.52             |21              |

| benders1.txt | nombre total d'itérations  | temps total (s) | valeur objective |
|----------------------------|----------------------------|-----------------|-----------------|
|EX1                         |54                         |3.47             |688              |
|EX2                         |78                         |1.42             |688              |
|Dijkstra                    |-                           |0.0             |688              |

| benders2.txt | nombre total d'itérations  | temps total (s) | valeur objective |
|----------------------------|----------------------------|-----------------|-----------------|
|EX1                         |112                         |10.0             |390              |
|EX2                         |84                         |2.91             |390              |
|Dijkstra                    |-                           |0.02             |390              |

| benders3.txt | nombre total d'itérations  | temps total (s) | valeur objective |
|----------------------------|----------------------------|-----------------|-----------------|
|EX1                         |128                         |16.21             |579              |
|EX2                         |117                         |2.46             |579              |
|Dijkstra                    |-                           |0.0             |579              |

| benders4.txt | nombre total d'itérations  | temps total (s) | valeur objective |
|----------------------------|----------------------------|-----------------|-----------------|
|EX1                         |119                         |5.36             |93              |
|EX2                         |70                         |1.39             |93              |
|Dijkstra                    |-                           |0.0             |93              |

My important paragraph.
{: .alert .alert-info}
* Pour le petit instance "benders-graphe-hexagone.txt", algo Dijkstra est moins rapide que la décomposition benders. Par contre, pour les grands instances, Dijkstra trouve la solution optimale tout de suite.

* Pour les grands instances, EX2 avec le problème maître relaché est plus efficace que le EX1 DB classique.

## Avec **bande passant de 3**, on obtient les résultats suivants : 


| benders-graphe-hexagone.txt | nombre total d'itérations  | temps total (s) | valeur objective |
|----------------------------|----------------------------|-----------------|-----------------|
|EX1                         |9                         |0.39             |21              |
|EX2                         |9                         |0.21             |21              |

| benders1.txt | nombre total d'itérations  | temps total (s) | valeur objective |
|----------------------------|----------------------------|-----------------|-----------------|
|EX1                         |54                         |3.66             |688              |
|EX2                         |78                         |0.99             |688              |

| benders2.txt | nombre total d'itérations  | temps total (s) | valeur objective |
|----------------------------|----------------------------|-----------------|-----------------|
|EX1                         |112                         |6.8             |390              |
|EX2                         |84                         |3.21             |390              |

| benders3.txt | nombre total d'itérations  | temps total (s) | valeur objective |
|----------------------------|----------------------------|-----------------|-----------------|
|EX1                         |128                         |11.76             |579              |
|EX2                         |117                         |2.89             |579              |

| benders4.txt | nombre total d'itérations  | temps total (s) | valeur objective |
|----------------------------|----------------------------|-----------------|-----------------|
|EX1                         |119                         |5.47             |93              |
|EX2                         |70                         |1.38             |93              |



* En changeant à la bande passante 3, le nombre d'itérations sur les instances ne changent pas.

* Par contre, concernant du temps d'exécution, le EX2 avec problème maître relaché est beaucoup plus efficace que la DB classique.


# Optimalité du le plus court chemin

:::danger
** Claim ** : le problème original est équivalent au problème le plus court chemin quand la bande passante est fixée à 1.
:::

:::info
** Proof ** :
:::