\pi: cantidad de kilogramos de harina producidos en el mes i = 1, 2, 3
\ai: cantidad de kilogramos de harina almacendos al finalizar el mes i = 1, 2, 3

Min 8p1 + 8.1p2 + 8.3p3+          \costo de produccion
    0.4a1 + 0.405 a2 + 0.415 a3   \costo de almacenamiento


Subject to

\Balance mes 1:
p1 - a1 = 5000
\Balance mes 2:
p2 + a1 - a2 = 35000
\Balance mes 3:
p3 + a2 - a3 = 30000

Bounds
\No es necesario 0 <=, lo toma por defecto 
\Maxima produccion:
p1 <= 30000
p2 <= 30000
p3 <= 30000
\Maximo almacenaje:
a1 <= 10000
a2 <= 10000
a3 <= 10000
\Minimo inventario mes 3:
a3 >= 4000
End
