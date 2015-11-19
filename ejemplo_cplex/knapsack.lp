\ El backslash sirve para comentarios de una linea.
\ Instancia del Knapsack Problem.
\ Ejemplo 2.1 del libro Knapsack Problems, Toth y Martello, pag. 21.
\ n = 8.
\ (p_j) = {15, 100, 90, 60, 40, 15, 10, 1).
\ (w_j) = {2, 20, 20, 30, 40, 30, 60, 10}.
\ c = 102.
\ Modelo (considerando variables continuas):
\ max  \sum_{j = 1}^n p_j x_j
\ sujeto a:
\      \sum_{j = 1}^n w_j x_j <= c
\      x_j >= 0, j = 1,\dots,n.
\ Definimos la funcion objetivo. "obj" es el nombre que le ponemos a la funcion. Se lo puede omitir.
Maximize
obj: 15 x1 + 100 x2 + 90 x3 + 60 x4 + 40 x5 + 15 x6 + 10 x7 + 1 x8
\ Definimos la restriccion de capacidad. cap es el nombre que le ponemos a la restriccion. Tambien se lo puede omitir.
\ Para agregar mas de una restriccion, repetir la misma idea en una nueva linea.
Subject To
cap: 2 x1 + 20 x2 + 20 x3 + 30 x4 + 40 x5 + 30 x6 + 60 x7 + 10 x8 <= 102
\ Definimos las cotas sobre las variables. El formato es lb <= x <= ub. Cada definicion debe
\ hacerse en una nueva linea. Si no se pone nada, por default toma 0 <= x <= +inf.
Bounds
0 <= x1 <= 1
0 <= x2 <= 1
0 <= x3 <= 1
0 <= x4 <= 1
0 <= x5 <= 1
0 <= x6 <= 1
0 <= x7 <= 1
0 <= x8 <= 1
\ Linea para indicar el fin de archivo.
End
