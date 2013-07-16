#include <vector>
#include <cstdio>
#include <string>
#include <algorithm>
#include "mtrand.h"
#include <cmath>
#include <cstring>
const double PI = 3.14159265358979323846;
using namespace std;
extern "C" {
void dgesv_(int* n, int* nrhs, double* A, int* lda, int* ipiv, double* B, int* ldb, int* info);
}

bool inSimplex(vector<vector<double> >& simplex, vector<double>& vals, vector<double>& p, double& fval) {
    int dim = p.size();
    double* A = new double[(dim + 1) * (dim + 1)];
    for(int i = 0; i < dim + 1; ++i) {
        for(int j = 0; j < dim; ++j) {
            A[i * (dim + 1) + j] = simplex[i][j];
        }
        A[i * (dim + 1) + dim] = 1;
    }
    /*
    for(int i = 0; i < (dim + 1) * (dim + 1); ++i) {
        printf("%lf\n", A[i]);
    }
    */
    double* B = new double[dim + 1];
    for(int i = 0; i < dim; ++i) {
        B[i] = p[i];
    }
    B[dim] = 1;
    /*
    for(int i = 0; i < dim + 1; ++i) {
        printf("%lf\n", B[i]);
    }
    */
    char trans = 'N';
    int n = dim + 1;
    int lda = dim + 1;
    int nrhs = 1;
    int* ipiv = new int[dim + 1];
    int ldb = dim + 1;
    int info;
    dgesv_(&n, &nrhs, A, &lda, ipiv, B, &ldb, &info);
    //printf("%d\n", info);
    bool flag = true;
    for(int i = 0; i < dim + 1; ++i) {
        //printf("%lf\n", B[i]);
        if(B[i] < 0) {
            flag = false;
            break;
        }
    }
    if(flag) {
        fval = 0;
        for(int i = 0; i < dim + 1; ++i) {
            fval += B[i] * vals[i];
        }
    }
    delete [] A;
    delete [] B;
    delete [] ipiv;
    return flag;
}

bool isWithin(const vector<double>& p, const vector<vector<double> >& bounds) {
	bool flag = true;
	for(int i = 0; i < p.size(); ++i) {
		if(p[i] - bounds[i][0] < 0 || p[i] - bounds[i][1] > 0) {
			flag = false;
			break;
		}
	}
	return flag;
}

void ReadFile(string PartitionFile, string ValFile, string SimplexFile, string PointFile, string CountFile, string SimLocFile, vector<vector<vector<double> > >& bounds, vector<vector<int> >& simpleces, vector<vector<double> >& p, vector<double>& vals, vector<double>& VolInCubes, vector<vector<int> >& SimLoc, vector<int>& SimCubeMap) {
    FILE* pParFile = fopen(PartitionFile.c_str(), "r");
    FILE* pValFile = fopen(ValFile.c_str(), "r");
    FILE* pSimFile = fopen(SimplexFile.c_str(), "r");
    FILE* pPntFile = fopen(PointFile.c_str(), "r");
    FILE* pCntFile = fopen(CountFile.c_str(), "r");
    FILE* pSimLocFile = fopen(SimLocFile.c_str(), "r");
    int dim, num;
    fscanf(pParFile, "%d%d", &dim, &num);
    bounds.resize(num);
    VolInCubes.resize(num);
    for(int i = 0; i < num; ++i) {
        bounds[i].resize(dim);
		for(int j = 0; j < dim; ++j) {
			bounds[i][j].resize(2);
			fscanf(pParFile, "%lf", &bounds[i][j][0]);
			fscanf(pParFile, "%lf", &bounds[i][j][1]);
		}
		fscanf(pParFile, "%lf", &VolInCubes[i]);     
    }
    int simNum, valNum, PointNum;
    fscanf(pCntFile, "%d%d%d", &simNum, &valNum, &PointNum);
    //printf("%d%d%d\n", simNum, valNum, PointNum);
    simpleces.resize(simNum);
    SimCubeMap.resize(simNum);
    SimLoc.resize(num);
    for(int i = 0; i < simNum; ++i) {
        simpleces[i].resize(dim + 1);
        int tmpInd;
        fscanf(pSimLocFile, "%d", &tmpInd);
        SimLoc[tmpInd].push_back(i);
        SimCubeMap[i] = tmpInd;
        for(int j = 0; j < dim + 1; ++j) {
            fscanf(pSimFile, "%d", &simpleces[i][j]);
            //printf("%10d", simpleces[i][j]);
        }
        //printf("\n");
    }
    /*
    for(int i = 0; i < SimLoc.size(); ++i) {
    	for(int j = 0; j < SimLoc[i].size(); ++j) {
    		printf("%d, %d\n", i, SimLoc[i][j]);
    	}
    }
    */
    printf("enterval");
    vals.resize(valNum);
    for(int i = 0; i < valNum; ++i) {
        fscanf(pValFile, "%lf", &vals[i]);
        //printf("%20.10lf\n", vals[i]);
    }
    printf("enterp");
    p.resize(PointNum);
    for(int i = 0; i < PointNum; ++i) {
        p[i].resize(dim);
        for(int j = 0; j < dim; ++j) {
            fscanf(pPntFile, "%lf", &p[i][j]);
            //printf("%20.10lf", p[i][j]);
        }
        //printf("\n");
    }
    fclose(pParFile);
    fclose(pValFile);
    fclose(pSimFile);
    fclose(pPntFile);
    fclose(pCntFile);
    fclose(pSimLocFile);
}

double Evaluate(int SimNum, double(*Density)(vector<double>&), vector<vector<vector<double> > >& bounds, vector<vector<int> >& simpleces, vector<vector<double> >& p, vector<double>& vals, vector<double>& VolInCubes, vector<vector<int> >& SimLoc, vector<int>& SimCubeMap) {
	int dim = simpleces[0].size() - 1;
	vector<double> TrueFval(SimNum, 0);
	vector<double> OPTFval(SimNum, 0);
	vector<double> FITFval(SimNum, 0);
	vector<vector<double> > simplex;
	vector<double> simVal;
	double fVal;
	simplex.resize(dim + 1);
	simVal.resize(dim + 1);
	//MTRand drand;
	FILE* fp = fopen("DENSITY", "r");
	for(int i = 0; i < SimNum; ++i) {
		vector<double> RandP(dim, 0);
		/*
		for(int j = 0; j < dim; ++j) {
			RandP[j] = drand();
		}
		TrueFval[i] = Density(RandP);
		*/
		for(int j = 0; j < dim; ++j) {
			fscanf(fp, "%lf", &RandP[j]);
		}
		fscanf(fp, "%lf", &TrueFval[i]);
		/*
		for(int j = 0; j < simpleces.size(); ++j) {
			for(int k = 0; k < simpleces[j].size(); ++k) {
				simplex[k] = p[simpleces[j][k]];
				simVal[k] = vals[simpleces[j][k]];
			}
			if(inSimplex(simplex, simVal, RandP, fVal)) {
				FITFval[i] = fVal;
				break;
			}
		}
		*/
		int cubeLoc;
		for(int j = 0; j < bounds.size(); ++j) {
			if(isWithin(RandP, bounds[j])) {
				OPTFval[i] = VolInCubes[j];
				cubeLoc = j;
				break;
			}
		}
		
		for(int k = 0; k < SimLoc[cubeLoc].size(); ++k) {
			for(int l = 0; l < simpleces[SimLoc[cubeLoc][k]].size(); ++l) {
				simplex[l] = p[simpleces[SimLoc[cubeLoc][k]][l]];
				simVal[l] = vals[simpleces[SimLoc[cubeLoc][k]][l]];
			}
			if(inSimplex(simplex, simVal, RandP, fVal)) {
				FITFval[i] = fVal;
				break;
			}
		}
		
		//printf("%10.5lf%10.5lf%10.5lf%10.5lf%10.5lf\n", RandP[0], RandP[1], TrueFval[i], FITFval[i], OPTFval[i]);
	}
	double Hellinger[2] = {0, 0};
	for(int i = 0; i < SimNum; ++i) {
		Hellinger[0] += pow(sqrt(TrueFval[i]) - sqrt(OPTFval[i]), 2.) / 2 / SimNum;
		Hellinger[1] += pow(sqrt(TrueFval[i]) - sqrt(FITFval[i]), 2.) / 2 / SimNum;
	}
	fclose(fp);
	
	//FILE* fpp = fopen((string("Hellinger") + file).c_str(), "w");
	printf("Hellinger Distance:\n OPT: %20.10lf\nFIT: %20.10lf\n", Hellinger[0], Hellinger[1]);
	//fclose(fpp);
	
}

double Normal_One(double p, double mean, double sd) {
	return 1. / sqrt(2 * PI) / sd * exp(-pow(p - mean, 2) / 2 / sd / sd);
}

double Mix_Normal(vector<double>& p) {
	double mean[2] = {0.2, 0.7};
	double sd[2] = {0.1, 0.05};
	double res[2] = {1, 1};
	for(int i = 0; i < p.size(); ++i) {
		res[0] *= Normal_One(p[i], mean[0], sd[0]);
	}
	for(int i = 0; i < p.size(); ++i) {
		res[1] *= Normal_One(p[i], mean[1], sd[1]);
	}
	return (res[0] + res[1]) / 2;
}

double pnorm(double x) {
	return exp(-x*x/2)/sqrt(2*M_PI);
}

double pnorm(double x, double mu, double sigma = 1) {
	return pnorm((x-mu)/sigma)/sigma;
}
//resultnormal
double sample_normal_density(vector<double> &x) {
   double density = 1.0;
   for (int i = 0; i < (int)x.size(); i++) {
       density *= pnorm(x[i], 0.35, 0.1);
   }
   return density;
}

//resultmixnormal
double sample_mix2normal_density(vector<double> &x) {
   	double mu1=0.3, mu2 =0.6, mu3=0.7, mu4=0.4, sd1=0.05, sd2=0.05, ratio=0.4;
   	double density = ratio*pnorm(x[0],mu1,sd1)*pnorm(x[1],mu2,sd2)+(1-ratio)*pnorm(x[0], mu3,sd1)*pnorm(x[1], mu4,sd2);

   	return density;
}

int main(int argc, char** argv) {
	string PartitionFile = argv[1];
	string ValFile = argv[2];
	string SimplexFile = argv[3];
	string PointFile = argv[4];
	string CountFile = argv[5];
	string SimLocFile = argv[6];
	int SimNum = atoi(argv[7]);
	vector<vector<vector<double> > > bounds;
	vector<vector<int> > simpleces;
	vector<vector<double> > p;
	vector<double> vals;
	vector<double> VolInCubes;
	vector<vector<int> > SimLoc;
	vector<int> SimCubeMap;
	ReadFile(PartitionFile, ValFile, SimplexFile, PointFile, CountFile, SimLocFile, bounds, simpleces, p, vals, VolInCubes, SimLoc, SimCubeMap);
	printf("enter");
	Evaluate(SimNum, sample_mix2normal_density, bounds, simpleces, p, vals, VolInCubes, SimLoc, SimCubeMap);
	/*
    vector<vector<double> > simplex;
    simplex.resize(4);
    for(int i = 0; i < simplex.size(); ++i) {
        simplex[i].resize(3);
    }    
    simplex[0][0] = 0;
    simplex[0][1] = 0;
    simplex[0][2] = 0;
    simplex[1][0] = 1;
    simplex[1][1] = 0;
    simplex[1][2] = 0;
    simplex[2][0] = 0;
    simplex[2][1] = 1;
    simplex[2][2] = 0;
    simplex[3][0] = 0;
    simplex[3][1] = 0;
    simplex[3][2] = 1;
    vector<double> vals(4, 0);
    vals[0] = 1;
    vals[1] = 2;
    vals[2] = 4;
    vals[3] = 5;
    vector<double> points(3, 0);
    points[0] = 1./4;
    points[1] = 1./4;
    points[2] = 1./4;
    double fval;
    bool yes = inSimplex(simplex, vals, points, fval);
    printf("%d, %lf\n", yes, fval);
    */
    return 0;
}
