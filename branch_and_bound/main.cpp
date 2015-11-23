#include <ilcplex/ilocplex.h>
#include <ilcplex/cplex.h>
ILOSTLBEGIN
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <utility>

#define TOL 1E-05

pair <vector<vector<bool> >, vector<vector<bool> > > parsear_entrada(string input_file);

int main(int argc, char *argv[]) {
  string archivo_entrada(argv[1]);
  
  pair <vector<vector<bool> >, vector<vector<bool> > > grafo = parsear_entrada(archivo_entrada);
  vector<vector<bool> > adyacencias = grafo.first; // matriz de adyacencia
  vector<vector<bool> > particion = grafo.second;  // filas: subconjuntos de la particion. columnas: nodos.

  // Datos de la instancia
  // Variables binarias:
  //		* X_n_j = nodo n pintado con el color j? (son cant_nodos * cant_colores_disp variables)
  //		* W_j	= hay algun nodo pintado con el color j? (son cant_colores_disp variables)
  //
  // Orden de las (cant_nodos * cant_col + cant_col_disp) variables:
  //		X_0_0, X_0_1, ... , X_0_(cant_col_disp), X_1_0, ... , X_(cant_nodos)_(cant_col_disp), W_0, ... , W(cant_col_disp)
  
  int cant_nodos = adyacencias.size();
  int cant_colores_disp = particion.size(); // cant colores usados <= cant de subconjuntos de la particion
  
  int n = cant_nodos * cant_colores_disp + cant_colores_disp; // n = cant de variables
  double calorias[] = {170,50,300};
  double calcio[] = {3,400,40};
  double minCalorias = 2000;
  double maxCalorias = 2300;
  double minCalcio = 1200;
  double maxPan = 3;
  double minLeche = 2;
    
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
  lp = CPXcreateprob(env, &status, "instancia dieta");

    
  if (lp == NULL) {
    cerr << "Error creando el LP" << endl;
    exit(1);
  }

    
  // Definimos las variables. No es obligatorio pasar los nombres de las variables, pero facilita el debug. La info es la siguiente:
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
    objfun[i] = 0; 		// Estas var no figuran en la funcion objetivo
    xctype[i] = 'B'; 	// 'C' es continua, 'B' binaria, 'I' Entera. Para LP (no enteros), este parametro tiene que pasarse como NULL. No lo vamos a usar por ahora.
    colnames[i] = new char[10];
    sprintf(colnames[i], "X_"); // TODO: concatenar num de nodo y de color
  }
  
  // Defino las variables W_j
  for (int i = n - cant_colores_disp; i < n; i++) {
    ub[i] = 1;
    lb[i] = 0;
    objfun[i] = 1;
    xctype[i] = 'B'; // 'C' es continua, 'B' binaria, 'I' Entera. Para LP (no enteros), este parametro tiene que pasarse como NULL. No lo vamos a usar por ahora.
    colnames[i] = new char[10];
    sprintf(colnames[i], "W_"); // TODO: concatenar num de color
  }
  
  // Agrego las columnas.
  status = CPXnewcols(env, lp, n, objfun, lb, ub, NULL, colnames);
  
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


  // CPLEX por defecto minimiza. Le cambiamos el sentido a la funcion objetivo si se quiere maximizar.
  // CPXchgobjsen(env, lp, CPX_MAX);

  // Generamos de a una las restricciones.
  // Estos valores indican:
  // ccnt = numero nuevo de columnas en las restricciones.
  // rcnt = cuantas restricciones se estan agregando.
  // nzcnt = # de coeficientes != 0 a ser agregados a la matriz. Solo se pasan los valores que no son cero.

  // Restricciones a agregar:
  //	* 

  int ccnt = 0, rcnt = 3, nzcnt = 0; 

  char sense[] = {'G','L','G'}; // Sentido de la desigualdad. 'G' es mayor o igual y 'E' para igualdad.

  double *rhs = new double[rcnt]; // Termino independiente de las restricciones.
  int *matbeg = new int[rcnt]; //Posicion en la que comienza cada restriccion en matind y matval.
  int *matind = new int[3*n]; // Array con los indices de las variables con coeficientes != 0 en la desigualdad.
  double *matval = new double[3*n]; // Array que en la posicion i tiene coeficiente ( != 0) de la variable cutind[i] en la restriccion.

  // Podria ser que algun coeficiente sea cero. Pero a los sumo vamos a tener 3*n coeficientes. CPLEX va a leer hasta la cantidad
  // nzcnt que le pasemos.


  //Restriccion de minimas calorias (asume que no hay coeficientes nulos; en caso de haberlos, agregar "if calorias[i] =! 0 ...")
  matbeg[0] = nzcnt;
  rhs[0] = minCalorias;
  for (int i = 0; i < n; i++) {
     matind[nzcnt] = i;
     matval[nzcnt] = calorias[i];
     nzcnt++;
  }

  //Restriccion de maximas calorias
  matbeg[1] = nzcnt;
  rhs[1] = maxCalorias;
  for (int i = 0; i < n; i++) {
     matind[nzcnt] = i;
     matval[nzcnt] = calorias[i];
     nzcnt++;
  }

  //Restriccion de minimo calcio
  matbeg[2] = nzcnt;
  rhs[2] = minCalcio;
  for (int i = 0; i < n; i++) {
     matind[nzcnt] = i;
     matval[nzcnt] = calcio[i];
     nzcnt++;
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
  status = CPXwriteprob(env, lp, "dieta2.lp", NULL);
    
  if (status) {
    cerr << "Problema escribiendo modelo" << endl;
    exit(1);
  }
    
  // Tomamos el tiempo de resolucion utilizando CPXgettime.
  double inittime, endtime;
  status = CPXgettime(env, &inittime);

  // Optimizamos el problema.
  status = CPXlpopt(env, lp);

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
  if(solstat!=CPX_STAT_OPTIMAL){
     exit(1);
  }  
    
  double objval;
  status = CPXgetobjval(env, lp, &objval);
    
  if (status) {
    cerr << "Problema obteniendo valor de mejor solucion." << endl;
    exit(1);
  }
    
  cout << "Datos de la resolucion: " << "\t" << objval << "\t" << (endtime - inittime) << endl; 

  // Tomamos los valores de la solucion y los escribimos a un archivo.
  std::string outputfile = "dieta.sol";
  ofstream solfile(outputfile.c_str());


  // Tomamos los valores de todas las variables. Estan numeradas de 0 a n-1.
  double *sol = new double[n];
  status = CPXgetx(env, lp, sol, 0, n - 1);

  if (status) {
    cerr << "Problema obteniendo la solucion del LP." << endl;
    exit(1);
  }

    
  // Solo escribimos las variables distintas de cero (tolerancia, 1E-05).
  solfile << "Status de la solucion: " << statstr << endl;
  for (int i = 0; i < n; i++) {
    if (sol[i] > TOL) {
      solfile << "x_" << i << " = " << sol[i] << endl;
    }
  }

  
  delete [] sol;
  solfile.close();

  return 0;
}

pair <vector<vector<bool> >, vector<vector<bool> > > parsear_entrada(string input_file){
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
	
	vector<vector<bool> > particion(k, vector<bool>(n, false));
	 
	//Cargo particion
	for(unsigned int i = 0; i < k; i++){
		getline(file, line); getline(file, line); getline(file, line); getline(file, line);
		aux = file.get();

		getline(file, line, ';');
		char line_no_const[line.length()];
		for(int j = 0; j < line.length(); j++) line_no_const[j] = line[j];
		
		pch = strtok(line_no_const, " ,");
		while (pch != NULL){
			particion[i][atoi(pch)] = true;
			pch = strtok(NULL, " ,");
		}
		
		getline(file, line); getline(file, line); getline(file, line);
	}

	/*cout << "Particion:" << endl;
	for(int i = 0; i < particion.size(); i++){
		for(int j = 0; j < particion[i].size(); j++)
			cout << particion[i][j] << " ";
		cout << endl;
	}*/
	
	vector<vector<bool> > adyacencias(n, vector<bool>(n, false));
	
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
		adyacencias[x][y] = true;
		adyacencias[y][x] = true;
	}
	
	/*cout << endl << "Adyacencias:" << endl;
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++)
			cout << adyacencias[i][j] << " ";
		cout << endl;
	}*/
	
	pair <vector<vector<bool> >, vector<vector<bool> > > result;
	result.first = adyacencias;
	result.second = particion;
	
    file.close();
	
	return result;
	
}
