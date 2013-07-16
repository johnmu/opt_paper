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

void Evaluate(vector<vector<vector<double> > >& bounds, vector<vector<int> >& simpleces, vector<vector<double> >& p, vector<double>& vals, vector<double>& VolInCubes, vector<vector<int> >& SimLoc, vector<int>& SimCubeMap, vector<vector<double> >& samples) {
	int dim = simpleces[0].size() - 1;
	vector<vector<double> > simplex;
	vector<double> simVal;
	double fVal;
	double max = *max_element(vals.begin(), vals.end());
	simplex.resize(dim + 1);
	simVal.resize(dim + 1);
	int sim = 10000;
	samples.resize(sim);
	MTRand drand;
	for(int i = 0; i < sim; ++i) {
		bool flag = false;
		while(!flag) {
			vector<double> RandP(dim, 0);
			for(int j = 0; j < dim; ++j) {
				RandP[j] = drand();
			}
		
			int cubeLoc;
			for(int j = 0; j < bounds.size(); ++j) {
				if(isWithin(RandP, bounds[j])) {
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
					break;
				}
			}
			if(drand() <= fVal / max) {
				samples[i] = RandP;
				flag = true;
			}
		}
	}
}

int main(int argc, char** argv) {
	string PartitionFile = argv[1];
	string ValFile = argv[2];
	string SimplexFile = argv[3];
	string PointFile = argv[4];
	string CountFile = argv[5];
	string SimLocFile = argv[6];
	vector<vector<vector<double> > > bounds;
	vector<vector<int> > simpleces;
	vector<vector<double> > p;
	vector<double> vals;
	vector<double> VolInCubes;
	vector<vector<int> > SimLoc;
	vector<int> SimCubeMap;
	vector<vector<double> > samples;
	printf("enter");
	ReadFile(PartitionFile, ValFile, SimplexFile, PointFile, CountFile, SimLocFile, bounds, simpleces, p, vals, VolInCubes, SimLoc, SimCubeMap);
	Evaluate(bounds, simpleces, p, vals, VolInCubes, SimLoc, SimCubeMap, samples);
	FILE* fp = fopen("Samples", "w");
	for(int i = 0; i < samples.size(); ++i) {
		for(int j = 0; j < samples[i].size(); ++j) {
			fprintf(fp, "%20.15lf", samples[i][j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
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
