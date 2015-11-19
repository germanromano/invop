\XcA: metros de listones chicos fabricados en la maquina A
\XcB: metros de listones chicos fabricados en la maquina B
\XcC: metros de listones chicos fabricados en la maquina C
\XmA: metros de listones medianos fabricados en la maquina A
\XmB: metros de listones medianos fabricados en la maquina B
\XmC: metros de listones medianos fabricados en la maquina C
\XgA: metros de listones grandess fabricados en la maquina A
\XgB: metros de listones grandess fabricados en la maquina B
\XgC: metros de listones grandess fabricados en la maquina C
\XeA: metros de listones extragrandes fabricados en la maquina A
\XeB: metros de listones extragrandes fabricados en la maquina B
\XeC: metros de listones extragrandes fabricados en la maquina C

Min 0.1 XcA + 0.12 XmA + 0.15 XgA + 0.3 XeA +         \costo hs usadas de maquina A
    0.1 XcB + 0.125 XmB + 0.15625 XgB + 0.25 XeB +    \costo hs usadas de maquina B
    0.1 XcC + 0.125XmC + 0.16 XgC + 0.32 XeC          \costo hs usadas de maquina C

Subject to

\Restricciones de demanda:
XcA + XcB + XcC = 10000
XmA + XmB + XmC = 8000
XgA + XgB + XgC = 6000
XeA + XeB + XeC = 4000

\Restricciones de disponibilidad:
1 XcA + 1.2 XmA + 1.5 XgA + 3 XeA <=  15000
1 XcB + 1.25 XmB + 1.5625 XgB + 2.5 XeB <=  25000 
1 XcC + 1.25 XmC + 1.6 XgC + 3.2 XeC <=  24000

Bounds
\No es necesario 0 <=, lo toma por defecto 
XmA <= 7000
End
