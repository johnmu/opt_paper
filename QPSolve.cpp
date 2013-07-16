#include "QPSolve.h"
#include <ilcplex/ilocplex.h>
#include <numeric>
#include <cmath>
#include <cstdio>
using namespace std;

ILOSTLBEGIN

extern "C" {
	void dgetrf_(const int* m, const int* n, double* A, const int* lda, int* ipiv, int* info);
}

double Determinant(vector<vector<double> >& mat) {
	int m = mat.size();
	int n = mat.size();
	int* ipiv = new int[mat.size()];
	int lda = mat.size();
	double* A = new double[mat.size() * mat.size()];
	for(int i = 0; i < mat.size(); ++i) {
		for(int j = 0; j < mat.size(); ++j) {
			A[i + mat.size() * j] = mat[i][j];
		}
	}
	int info;
	dgetrf_(&m, &n, A, &lda, ipiv, &info);
	/*
	for(int i = 0; i < mat.size(); ++i) {
		for(int j = 0; j < mat.size(); ++j) {
			printf("%15.8lf", A[i + mat.size() * j]);
		}
		printf("\n");
	}
	for(int i = 0; i < mat.size(); ++i) {
		printf("%10d", ipiv[i]);
	}
	printf("\n");
	*/
	int numOfChange = 0;
	for(int i = 0; i < mat.size(); ++i) {
		if(ipiv[i] != i + 1) {
			++numOfChange;
		}
	}
	
	//printf("%d\n", numOfChange % 2);
	double det = 1.;
	for(int i = 0; i < mat.size(); ++i) {
		det *= A[i + mat.size() * i];
	}
	if(numOfChange % 2 == 1) {
		det = -det;
	}
	//printf("%15.8lf\n", det);
	delete [] ipiv;
	delete [] A;
	return det;
}

void QPSolve(vector<vector<double> >& p, vector<vector<int> >& simpleces, vector<double>& volumes, vector<double>& simplexArea, vector<double>& vals, double lambda, vector<double>& ratios) {
	vals.resize(p.size());
	double sum;
	sum = accumulate(simplexArea.begin(), simplexArea.end(), 0.0);
	for(int i = 0; i < simplexArea.size(); ++i) {
		simplexArea[i] /= sum;
	}
	IloEnv env;
	try {
		IloModel model(env);
		IloNumVarArray var(env);
		IloRangeArray con(env);
		for(int i = 0; i < vals.size(); ++i) {
			var.add(IloNumVar(env));
		}

		IloExpr objective(env);
		IloExpr constraints(env);
		//least square part
		for(int i = 0; i < simpleces.size(); ++i) {
			IloExpr tmp(env);
			for(int j = 0; j < simpleces[i].size(); ++j) {
				tmp += var[simpleces[i][j]];
			}
			constraints += simplexArea[i] * tmp / simpleces[i].size();
			objective += (simplexArea[i] * tmp / simpleces[i].size() - volumes[i] * simplexArea[i]) * (simplexArea[i] * tmp / simpleces[i].size() - volumes[i] * simplexArea[i]) / simplexArea[i];
			//objective += (tmp / simpleces[i].size() - volumes[i]) * (tmp / simpleces[i].size() - volumes[i]);
		}
		//penalty part
		vector<vector<double> > matDeterminant;
		vector<vector<double> > mat;
		matDeterminant.resize(p[0].size());
		mat.resize(p[0].size());
		for(int i = 0; i < mat.size(); ++i) {
			mat[i].resize(mat.size());
		}
		for(int i = 0; i < matDeterminant.size(); ++i) {
			matDeterminant[i].resize(simpleces[0].size());
		}
		int factor = 1;
		for(int i = 0; i < p[0].size(); ++i) {
			factor *= i + 1;
		}
		for(int i = 0; i < simpleces.size(); ++i) {
			for(int j = 0; j < matDeterminant.size(); ++j) {
				for(int k = 0; k < matDeterminant[j].size(); ++k) {
					int row = 0; 
					for(int ii = 0; ii < simpleces[i].size(); ++ii) {
						if(ii != k) {
							int col = 0;
							for(int jj = 0; jj < p[simpleces[i][ii]].size(); ++jj) {
								if(jj != j) {
									mat[row][col++] = p[simpleces[i][ii]][jj];
								}
							}
							mat[row++][col] = 1.;
						}
					}
					matDeterminant[j][k] = Determinant(mat);
				}
			}
			IloExpr Area(env);
			for(int j = 0; j < matDeterminant.size(); ++j) {
				IloExpr tmp(env);
				for(int k = 0; k < matDeterminant[j].size(); ++k) {
					int sign = k % 2 == 0 ? 1 : -1;
					tmp += sign * var[simpleces[i][k]] * matDeterminant[j][k];
				}
				Area += tmp * tmp;
			}
			objective += Area / factor / factor / simplexArea[i] / simplexArea[i] * lambda;
		}
		model.add(IloMinimize(env, objective));
		con.add(constraints == 1);
		model.add(con);
		
		IloCplex cplex(model);
		//cplex.setParam(IloCplex::TiLim, 3600);
		//cplex.setParam(IloCplex::EpGap, 0.05);
		//cplex.setParam(IloCplex::BarAlg, 1);
		//cplex.setParam(IloCplex::SolutionTarget, 2);
		//env.out() << "SolutionTarget: " << cplex.getParam(IloCplex::SolutionTarget) << endl;
		cplex.solve();
		
      	//env.out() << "Solution status = " << cplex.getStatus() << endl;
      	//env.out() << "Solution value  = " << cplex.getObjValue() << endl;

      	IloNumArray optimal(env);
      	cplex.getValues(optimal, var);
      	//env.out() << "optimal        = " << optimal << endl;
      	for(int i = 0; i < vals.size(); ++i) {
      		vals[i] = optimal[i];
      	}
      	//cplex.getSlacks(vals, con);
      	//env.out() << "Slacks        = " << vals << endl;

     	//cplex.exportModel("miqpex1.lp");
   	}
   	catch (IloException& e) {
    	  cerr << "Concert exception caught: " << e << endl;
   	}
   	catch (...) {
    	  cerr << "Unknown exception caught" << endl;
   	}

   	env.end();
	
}
