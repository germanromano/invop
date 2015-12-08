#include <ilcplex/ilocplex.h>
#include <ilcplex/cplex.h>
ILOSTLBEGIN
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <utility>

#define TOL 1E-05

pair <int, pair<vector<vector<bool> >*, vector<vector<bool> >* > > parsear_entrada(string input_file);

int main(int argc, char *argv[]) {
	
	if(argc < 2){
		cerr << "Uso: input_file" << endl;
		exit(1);
	}
	
	string archivo_entrada(argv[1]);
	
	pair <int, pair<vector<vector<bool> >*, vector<vector<bool> >* > > grafo = parsear_entrada(archivo_entrada);
	int cant_ejes = grafo.first;
	vector<vector<bool> > *adyacencias = grafo.second.first; // matriz de adyacencia
	vector<vector<bool> > *particion = grafo.second.second;	// filas: subconjuntos de la particion. columnas: nodos.
	
	// Variables binarias:
	//		* X_n_j = nodo n pintado con el color j? (son cant_nodos * cant_colores_disp variables)
	//		* W_j	= hay algun nodo pintado con el color j? (son cant_colores_disp variables)
	//			=> TOTAL: (cant_nodos * cant_colores_disp + cant_colores_disp) variables
	//
	// Orden de las variables:
	//		X_0_0, X_0_1, ... , X_0_(cant_col_disp), X_1_0, ... , X_(cant_nodos)_(cant_col_disp), W_0, ... , W(cant_col_disp)

	
	int cant_nodos = adyacencias->size();
	int cant_subconj_particion = particion->size(); //cant de subconjuntos de la particion
	int cant_colores_disp = particion->size(); // cant colores usados <= cant de subconjuntos de la particion
	
	int n = cant_nodos * cant_colores_disp + cant_colores_disp; // n = cant de variables
	
	// Genero el problema de cplex.
	int status;
	CPXENVptr env; // Puntero al entorno.
	CPXLPptr lp; // Puntero al LP
	 
	// Creo el entorno.
	env = CPXopenCPLEX(&status);

		
	if (env == NULL) {
		cerr << "Error creando el entorno" << endl;
		exit(1);
	}
		
	// Creo el LP.
	lp = CPXcreateprob(env, &status, "Coloreo Particionado");

		
	if (lp == NULL) {
		cerr << "Error creando el LP" << endl;
		exit(1);
	}

	//TUNNING
	//Para que haga Branch & Cut:
	CPXsetintparam(env, CPX_PARAM_MIPSEARCH, CPX_MIPSEARCH_TRADITIONAL);
	//Para que no se adicionen planos de corte: ( => Branch & Bound)
	CPXsetintparam(env,CPX_PARAM_EACHCUTLIM, 0);
	CPXsetintparam(env, CPX_PARAM_FRACCUTS, -1);
	//Para facilitar la comparación evitamos paralelismo:
	CPXsetintparam(env, CPX_PARAM_THREADS, 1);
	//Para desactivar preprocesamiento
	CPXsetintparam(env, CPX_PARAM_PRESLVND, -1);
	CPXsetintparam(env, CPX_PARAM_REPEATPRESOLVE, 0);
	CPXsetintparam(env, CPX_PARAM_RELAXPREIND, 0);
	CPXsetintparam(env, CPX_PARAM_REDUCE, 0);
	CPXsetintparam(env, CPX_PARAM_LANDPCUTS, -1);

	double *ub, *lb, *objfun; // Cota superior, cota inferior, coeficiente de la funcion objetivo.
	char *xctype, **colnames; // tipo de la variable (por ahora son siempre continuas), string con el nombre de la variable.
	ub = new double[n]; 
	lb = new double[n];
	objfun = new double[n];
	xctype = new char[n];
	colnames = new char*[n];
	
	// Defino las variables X_n_j
	for (int i = 0; i < n - cant_colores_disp; i++) {
		ub[i] = 1;
		lb[i] = 0;
		objfun[i] = 0; // Estas var no figuran en la funcion objetivo
		xctype[i] = 'B';
		colnames[i] = new char[10];
		sprintf(colnames[i], "X_%d_%d", i / cant_colores_disp, i % cant_colores_disp);
	}

	// Defino las variables W_j
	for (int i = n - cant_colores_disp; i < n; i++) {
		ub[i] = 1;
		lb[i] = 0;
		objfun[i] = 1;
		xctype[i] = 'B';
		colnames[i] = new char[10];
		sprintf(colnames[i], "W_%d", i - (n - cant_colores_disp));
	}
	
	// Agrego las columnas.
	status = CPXnewcols(env, lp, n, objfun, lb, ub, xctype, colnames);
	
	if (status) {
		cerr << "Problema agregando las variables CPXnewcols" << endl;
		exit(1);
	}
	
	// Libero las estructuras.
	for (int i = 0; i < n; i++) {
		delete[] colnames[i];
	}
	
	delete[] ub;
	delete[] lb;
	delete[] objfun;
	delete[] xctype;
	delete[] colnames;

	// ccnt = numero nuevo de columnas en las restricciones.
	// rcnt = cuantas restricciones se estan agregando.
	// nzcnt = # de coeficientes != 0 a ser agregados a la matriz. Solo se pasan los valores que no son cero.

	// Restricciones:
	//	(1) Nodos adyacentes tienen distinto color (cant_ejes * cant_colores_disp restricciones por <=)
	//	(2) Cada nodo tiene a lo sumo un color (cant_nodos restricciones por <=)
	//	(3) Solo un nodo de cada subconj. de la particion tiene color (cant. de subconj. de la particion restricciones por =)
	//	(4) W_j = 1 sii "X_i_j = 1 para algún i" (cant_colores_disp restricciones por >=)
	//		=> TOTAL: (cant_ejes * cant_colores_disp + cant_nodos + cant_subconj_particion + cant_colores_disp) restricciones

	int ccnt = 0;
	int rcnt = cant_ejes * cant_colores_disp + cant_nodos + cant_subconj_particion + cant_colores_disp;
	int nzcnt = 0; 

	char sense[rcnt]; // Sentido de la desigualdad. 'G' es mayor o igual y 'E' para igualdad.
	for(unsigned int i = 0; i < cant_ejes * cant_colores_disp; i++)
		sense[i] = 'L';
	for(unsigned int i = cant_ejes * cant_colores_disp; i < cant_ejes * cant_colores_disp + cant_nodos; i++)
		sense[i] = 'L';
	for(unsigned int i = cant_ejes * cant_colores_disp + cant_nodos; i < cant_ejes * cant_colores_disp + cant_nodos + cant_subconj_particion; i++)
		sense[i] = 'E';
	for(unsigned int i = cant_ejes * cant_colores_disp + cant_nodos + cant_subconj_particion; i < rcnt; i++)
		sense[i] = 'G';

	double *rhs = new double[rcnt]; // Termino independiente de las restricciones.
	int *matbeg = new int[rcnt]; //Posicion en la que comienza cada restriccion en matind y matval.
	int *matind = new int[rcnt*n]; // Array con los indices de las variables con coeficientes != 0 en la desigualdad.
	double *matval = new double[rcnt*n]; // Array que en la posicion i tiene coeficiente ( != 0) de la variable cutind[i] en la restriccion.

	//El termino indep. de restr (1), (2) y (3) es 1
	for(unsigned int i = 0; i < rcnt - cant_colores_disp; i++)
		rhs[i] = 1;
		
	//El termino indep. de restr (4) es 0
	for(unsigned int i = rcnt - cant_colores_disp; i < rcnt; i++)
		rhs[i] = 0;
	
	unsigned int indice = 0; //numero de restriccion actual
	
	//Restricciones (1)
	for(unsigned int i = 0; i < cant_nodos; i++) //itero nodo 1
		for(unsigned int j = i+1; j < cant_nodos; j++) //itero nodo 2
			if((*adyacencias)[i][j])
				for(unsigned int p = 0; p < cant_colores_disp; p++){ //itero color
					matbeg[indice] = nzcnt;
					indice++;
					//cargo una de las variables participantes de la restr.
					matind[nzcnt] = cant_colores_disp*i + p; //var1: X_nodo1_color
					matval[nzcnt] = 1;
					nzcnt++;
					//idem con la otra variable
					matind[nzcnt] = cant_colores_disp*j + p; //var2: X_nodo2_color
					matval[nzcnt] = 1;
					nzcnt++;
				}
				
	//Restricciones (2)
	for(unsigned int i = 0; i < cant_nodos; i++){ //itero nodo
		matbeg[indice] = nzcnt;
		indice++;
		for(unsigned int p = 0; p < cant_colores_disp; p++){ //itero color
			matind[nzcnt] = cant_colores_disp*i + p; //var: X_nodo_color
			matval[nzcnt] = 1;
			nzcnt++;
		}
	}
	
	//Restricciones (3)
	for(unsigned int v = 0; v < cant_subconj_particion; v++){ //itero subconjunto de la particion
		matbeg[indice] = nzcnt;
		indice++;
		for(unsigned int i = 0; i < cant_nodos; i++) //itero nodo
			if((*particion)[v][i])
				for(unsigned int p = 0; p < cant_colores_disp; p++){ //itero color
					matind[nzcnt] = cant_colores_disp*i + p; //var: X_nodo_color
					matval[nzcnt] = 1;
					nzcnt++;
				}
	}
	
	//Restricciones (4)
	for(unsigned int p = 0; p < cant_colores_disp; p++){ //itero color
		matbeg[indice] = nzcnt;
		indice++;
		matind[nzcnt] = cant_nodos * cant_colores_disp + p; //var: W_color
		matval[nzcnt] = cant_nodos;
		nzcnt++;
		for(unsigned int i = 0; i < cant_nodos; i++){ //itero nodo
			matind[nzcnt] = cant_colores_disp*i + p; //var: X_nodo_color
			matval[nzcnt] = -1;
			nzcnt++;
		}
	}
	
	// Esta rutina agrega la restriccion al lp.
	status = CPXaddrows(env, lp, ccnt, rcnt, nzcnt, rhs, sense, matbeg, matind, matval, NULL, NULL);
			
	if (status) {
		cerr << "Problema agregando restricciones." << endl;
		exit(1);
	}
			
	delete[] rhs;
	delete[] matbeg;
	delete[] matind;
	delete[] matval;

	// Seteo de algunos parametros.
	// Para desactivar la salida poner CPX_OFF.
	status = CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
		
	if (status) {
		cerr << "Problema seteando SCRIND" << endl;
		exit(1);
	}
		
	// Por ahora no va a ser necesario, pero mas adelante si. Setea el tiempo
	// limite de ejecucion.
	status = CPXsetdblparam(env, CPX_PARAM_TILIM, 3600);
	
	if (status) {
		cerr << "Problema seteando el tiempo limite" << endl;
		exit(1);
	}
 
	// Escribimos el problema a un archivo .lp.
	status = CPXwriteprob(env, lp, "output.lp", NULL);
		
	if (status) {
		cerr << "Problema escribiendo modelo" << endl;
		exit(1);
	}
	
	// Tomamos el tiempo de resolucion utilizando CPXgettime.
	double inittime, endtime;
	status = CPXgettime(env, &inittime);

	// Optimizamos el problema.
	status = CPXmipopt(env, lp);

	status = CPXgettime(env, &endtime);

	if (status) {
		cerr << "Problema optimizando CPLEX" << endl;
		exit(1);
	}

	// Chequeamos el estado de la solucion.
	int solstat;
	char statstring[510];
	CPXCHARptr p;
	solstat = CPXgetstat(env, lp);
	p = CPXgetstatstring(env, solstat, statstring);
	string statstr(statstring);
	cout << endl << "Resultado de la optimizacion: " << statstring << endl;
	
	if(solstat!=CPXMIP_OPTIMAL && solstat!=CPXMIP_OPTIMAL_TOL &&
		solstat!=CPXMIP_NODE_LIM_FEAS && solstat!=CPXMIP_TIME_LIM_FEAS){
			exit(1);
	}
	
	double objval;
	status = CPXgetmipobjval(env, lp, &objval);
		
	if (status) {
		cerr << "Problema obteniendo valor de mejor solucion." << endl;
		exit(1);
	}
		
	cout << "Optimo: " << objval << "\t(Time: " << (endtime - inittime) << " sec)" << endl; 

	// Tomamos los valores de la solucion y los escribimos a un archivo.
	std::string outputfile = "output.sol";
	ofstream solfile(outputfile.c_str());

	// Tomamos los valores de todas las variables. Estan numeradas de 0 a n-1.
	double *sol = new double[n];
	status = CPXgetmipx(env, lp, sol, 0, n - 1);

	if (status) {
		cerr << "Problema obteniendo la solucion del LP." << endl;
		exit(1);
	}

	// Solo escribimos las variables distintas de cero (tolerancia, 1E-05).
	solfile << "Status de la solucion: " << statstr << endl;
	// Imprimo var X_n_j
	for (int i = 0; i < n - cant_colores_disp; i++) {
		if (sol[i] > TOL) {
			solfile << "X_" << i / cant_colores_disp << "_" << i % cant_colores_disp << " = " << sol[i] << endl;
		}
	}
	// Imprimo var W_j
	for (int i = n - cant_colores_disp; i < n; i++) {
		if (sol[i] > TOL) {
			solfile << "W_" << i - (n - cant_colores_disp) << " = " << sol[i] << endl;
		}
	}

	delete [] sol;
	solfile.close();
	
	return 0;
}

pair <int, pair<vector<vector<bool> >*, vector<vector<bool> >* > > parsear_entrada(string input_file){
	unsigned int n, m, k;
	string aux, line;
	char* pch;
	ifstream file;
	
	file.open(input_file.c_str());

	//Cargo parametros
	file >> aux;
	n = atoi(aux.substr(2, aux.length()-2).c_str());
	file >> aux;
	m = atoi(aux.substr(2, aux.length()-2).c_str());
	file >> aux;
	k = atoi(aux.substr(2, aux.length()-2).c_str());
	
	vector<vector<bool> > *particion = new vector<vector<bool> >(k, vector<bool>(n, false));
	 
	//Cargo particion
	for(unsigned int i = 0; i < k; i++){
		getline(file, line); getline(file, line); getline(file, line); getline(file, line);
		aux = file.get();

		getline(file, line, ';');
		char line_no_const[line.length()];
		for(int j = 0; j < line.length(); j++) line_no_const[j] = line[j];
		
		pch = strtok(line_no_const, " ,");
		while (pch != NULL){
			(*particion)[i][atoi(pch)] = true;
			pch = strtok(NULL, " ,");
		}
		
		getline(file, line); getline(file, line); getline(file, line);
	}
	
	vector<vector<bool> > *adyacencias = new vector<vector<bool> >(n, vector<bool>(n, false));
	
	//Cargo adyacencias
	getline(file, line); getline(file, line);
	unsigned int x, y;
	
	for(unsigned int i = 0; i < m; i++){
		aux = file.get(); aux = file.get();
		getline(file, line, ';');
		char line_no_const[line.length()];
		for(int j = 0; j < line.length(); j++) line_no_const[j] = line[j];
		
		pch = strtok(line_no_const, " -");
		x = atoi(pch);
		pch = strtok(NULL, " -");
		y = atoi(pch);
		(*adyacencias)[x][y] = true;
		(*adyacencias)[y][x] = true;
	}
	
	pair <int, pair<vector<vector<bool> >*, vector<vector<bool> >* > > result(m, pair<vector<vector<bool> >*,
		vector<vector<bool> >* >(adyacencias, particion)); 
	
	file.close();
	return result;	
}
