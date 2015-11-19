\xm: cantidad de porciones diarias de maiz de nuestra dieta
\xl: cantidad de porciones diarias de leche de nuestra dieta
\xp: cantidad de porciones diarias de pan de nuestra dieta

Min 1.8xm + 2.3xl + 1.5xp  \ costo de compra

Subject to

170xm + 50xl + 300xp >= 2000 \min. calorias
170xm + 50xl + 300xp <= 2300 \max. calorias
3xm + 400xl + 40xp >= 1200 \min. calcio

Bounds
xp <= 3 \max. porciones pan
xl >= 2 \min. porciones leche
\No es necesario 0 <=, lo toma por defecto 
End

