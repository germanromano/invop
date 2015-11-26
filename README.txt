*Generador de grafos
	-parámetros: cant_nodos cant_particiones densidad output_file (0 <= densidad <= 1)

*Para plotear grafos .dot
	-instalar Graphviz
	-ejecutar comando: dot -Tpng input_file.dot -o output_file.png
	
*Coloreo Particionado con Branch and Bound
	-ejecutar Makefile
	-parámetro: path al grafo en formato .dot
