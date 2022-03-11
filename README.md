# MLA TP noté
 Problème localisation avec décomposition benders, réalisé par Yue Zhang.


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

<!---
:::success
-->
```diff
+ * Pour le petit instance "benders-graphe-hexagone.txt", algo Dijkstra est moins rapide que la décomposition benders. Par contre, pour les grands instances, Dijkstra trouve la solution optimale tout de suite.

+ * Pour les grands instances, EX2 avec le problème maître relaché est plus efficace que le EX1 DB classique.
```
<!---
:::
-->

## Avec **bande passant de 3**, on obtient les résultats suivants : 

| benders-graphe-hexagone.txt | nombre total d'itérations  | temps total (s) | valeur objective |
|----------------------------|----------------------------|-----------------|-----------------|
|EX1                         |10                         |0.32             |10              |
|EX2                         |9                         |0.16             |10              |

| benders1.txt | nombre total d'itérations  | temps total (s) | valeur objective |
|----------------------------|----------------------------|-----------------|-----------------|
|EX1                         |219                         |39.78             |237              |
|EX2                         |62                         |0.92             |235              |

| benders2.txt | nombre total d'itérations  | temps total (s) | valeur objective |
|----------------------------|----------------------------|-----------------|-----------------|
|EX1                         |MAXITE                         |199.05             |135              |
|EX2                         |102                         |2.88             |133              |

| benders3.txt | nombre total d'itérations  | temps total (s) | valeur objective |
|----------------------------|----------------------------|-----------------|-----------------|
|EX1                         |MAXITE                         |751.07             |201              |
|EX2                         |90                         |2.59             |198              |

| benders4.txt | nombre total d'itérations  | temps total (s) | valeur objective |
|----------------------------|----------------------------|-----------------|-----------------|
|EX1                         |MAXITE                         |4472.67             |38              |
|EX2                         |87                         |2.12             |35              |



Particulièrement, sur certains instances, la méthode EX2 avec le problème relâché n'a pas l'aire à s'arrêter, donc on fixe une limite itération ```MAXITE = 300```.

<!---
:::success
-->
```diff
+ * En changeant à la bande passante 3, le nombre d'itérations sur les instances ne changent pas.

+ * Par contre, concernant du temps d'exécution, le EX2 avec problème maître relaché est beaucoup plus efficace que la DB classique.
```

<!---
:::
-->


# Optimalité du le plus court chemin

<!---
:::danger
-->
```diff
- Le problème original est équivalent au problème le plus court chemin quand la bande passante est fixée à 1.
```
<!---
:::
-->

<!---
:::info
-->
Quand $b_{nd} = 1$, alors la contrainte (1d) devient $y_{ij} \geq x_{ij} + x_{ji},\forall \{ij\} \in E, i<j$. Donc la function objective devient à $\min \sum_{ij\in E} x_{ij}$ c'est la quantité de flux circule dans le réseaux. Comme il n'y a pas de contrainte capacité sur les arcs, pour chaque chemin $P_t$ de $s$ à $t\in T$, chaque arc apporte exactement la demande $d_t$ au terminal $t$. Donc le objective peut être reécrit comme $\min \sum_{t \in T} |P_t| \cdot d_t$. Donc on cherche le plus court chemin de $s$ à tous les terminaux.
<!---
:::
-->