Trabajo Práctico Final - Investigación Operativa

*Coloreo Particionado con Branch and Bound
	-ejecutar Makefile
	-parámetro: input_file.dot
	
*Coloreo Particionado con Planos de Corte o Cut and Branch
	-ejecutar Makefile
	-parámetros: input_file.dot tope_de_iteraciones
	
*Generador de grafos
	-parámetros: cant_nodos cant_particiones densidad output_file
	-con: 0 <= densidad <= 1

*Para plotear grafos .dot
	-instalar Graphviz
	-ejecutar comando: dot -Tpng input_file.dot -o output_file.png
