//GENERA GRAFOS NO DIRIGIDOS EN FORMATO DOT
//Uso: cant_nodos densidad output_file (0 <= densidad <= 1)

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <time.h>
#include <vector>

using namespace std;

vector< pair <int, int> > generar_grafo(int n, double p);

int main(int argc, char** argv){	
	if(argc < 3) cout << "Uso: cant_nodos densidad output_file (0 <= densidad <= 1)" << endl;  

	int n = atoi(argv[1]);         
	float p = atof(argv[2]);
	string archivoOut(argv[3]);
	srand (time(NULL));

	fstream ofs;
	ofs.open(archivoOut.c_str(), fstream::out);
	
    string graph_name;
    graph_name.append(archivoOut.substr(0, archivoOut.rfind(".")));
	
	vector< pair <int, int> > v = generar_grafo(n, p);
	int m = v.size();

	//FORMATO SNAP
	//ofs << "# Directed graph (each unordered pair of nodes is saved once):" << endl;
	//ofs << "# Example shown in the Bryan and Leise." << endl;
	//ofs << "# Nodes: " << n <<" Edges: " << m << endl;
	//ofs << "# FromNodeId	ToNodeId" << endl;
	
	//FORMATO DOT
	ofs << "//" << n << endl << "//" << m << endl;
	ofs << "graph " << graph_name << " {" << endl;
	
	for(unsigned int i = 0; i < v.size(); i++)
		ofs << "\t" << v[i].first << " -- " << v[i].second << ";" << endl;
		
	ofs << "}" << endl;

	return 0;
}

vector< pair <int, int> > generar_grafo(int n, double p){
	int dado;
	vector< pair <int, int> > v;
	
	for(int j = 1; j <= n; j++){
		for (int i = (j+1); i<=n; i++){
			dado = rand() % 1000 + 1; //random entre 1 y 1000
			if(dado <= p*1000){
				pair <int, int> link (j, i);
				v.push_back(link);
			}
			
			dado = rand() % 1000 + 1; //random entre 1 y 1000
			if(dado <= p*1000){
				pair <int, int> link (i, j);
				v.push_back(link);
			}
		}
	}
	return v;
}

