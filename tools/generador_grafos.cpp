//GENERA GRAFOS NO DIRIGIDOS EN FORMATO DOT
//Uso: cant_nodos cant_particiones densidad output_file (0 <= densidad <= 1)

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <time.h>
#include <vector>

using namespace std;

vector< pair <int, int> > generar_grafo(int n, double p);
vector< pair <int, int> > generar_particion(int n, int k);

string colores[] = {"red", "green", "blue", "yellow", "black", "chocolate", "gold", "pink", 
	"purple", "steelblue", "tomato", "seagreen", "firebrick", "cyan", "brown", "forestgreen", "indigo"};

int main(int argc, char** argv){	
	if(argc < 5){
		cerr << "Uso: cant_nodos cant_particiones densidad output_file (0 <= densidad <= 1)" << endl;  
		exit(1);
	}
	
	int n = atoi(argv[1]);
	int k = atoi(argv[2]);         
	float p = atof(argv[3]);
	string archivoOut(argv[4]);
	srand (time(NULL));

	fstream ofs;
	ofs.open(archivoOut.c_str(), fstream::out);
	
    string graph_name;
    graph_name.append(archivoOut.substr(0, archivoOut.rfind(".")));
	
	vector< pair <int, int> > particion = generar_particion(n, k);
	
	vector< pair <int, int> > aristas = generar_grafo(n, p);
	int m = aristas.size();

	//FORMATO SNAP
	//ofs << "# Directed graph (each unordered pair of nodes is saved once):" << endl;
	//ofs << "# Example shown in the Bryan and Leise." << endl;
	//ofs << "# Nodes: " << n <<" Edges: " << m << endl;
	//ofs << "# FromNodeId	ToNodeId" << endl;
	
	//FORMATO DOT
	//Primeras tres lineas: cant nodos (n), cant aristas (m), cant particiones (k)
	ofs << "//" << n << endl << "//" << m << endl << "//" << k << endl;
	ofs << "graph " << graph_name << " {" << endl;
	
	bool vacia = true;
	for(unsigned int i = 0; i < k; i++){
		for(unsigned int j = 0; j < n; j++)
			if(particion[j].second == i){
				if(!vacia) ofs << ", " << j;
				else{
					vacia = false;
					ofs << endl << "\tsubgraph cluster_" << i << "{" << endl << "\t";
					ofs << j;
					}
			}
		
		if(!vacia){
			ofs << ";" << endl;
			ofs << "\tlabel = \"V" << i << "\";" << endl;
			ofs << "\tcolor = red;" << endl << "\tpenwidth = 3;" << endl << "\t}" << endl;
			vacia = true;
		}
	}
	
	ofs << endl;
		
	for(unsigned int i = 0; i < aristas.size(); i++)
		ofs << "\t" << aristas[i].first << " -- " << aristas[i].second << ";" << endl;
		
	ofs << "}" << endl;

	return 0;
}

vector< pair <int, int> > generar_particion(int n, int k){
	int dado;
	vector< pair <int, int> > v;
	
	// Los primeros k nodos van uno a cada subconjunto, para que ninguno quede vacío
	for(int i = 0; i < k; i++){
		pair <int, int> ubicacion (i, i);
		v.push_back(ubicacion);
		}
	
	for(int i = k; i < n; i++){
		dado = rand() % k;
		pair <int, int> ubicacion (i, dado);
		v.push_back(ubicacion);
		}
	
	return v;
	}

vector< pair <int, int> > generar_grafo(int n, double p){
	int dado;
	vector< pair <int, int> > v;
	
	for(int j = 0; j < n; j++){
		for (int i = (j+1); i < n; i++){
			dado = rand() % 1000 + 1; //random entre 1 y 1000
			if(dado <= p*1000){
				pair <int, int> link (j, i);
				v.push_back(link);
			}
		}
	}
	return v;
}

