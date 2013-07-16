extern "C" {
#define qh_QHimport
#include "qhull_a.h"
void dgesvd_(const char* jobu, const char* jobvt, const int* M, const int* N,
        double* A, const int* lda, double* S, double* U, const int* ldu,
        double* VT, const int* ldvt, double* work,const int* lwork, const
        int* info);
}
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <climits>
#include <string>
#include <iostream>
#include "QPSolve.h"
using namespace std;

void delaunayn(const vector<vector<double> >& p, vector<vector<int> >& triangles) {
	unsigned dim, n;
	boolT ismalloc;
	char flags[250];
	
	double* pt_array;
	
	FILE* outfile;
	FILE* errfile = stderr;
	
	dim = p[0].size();
	n = p.size();
	
	if(n > dim + 1) {
		int exitcode;
		int curlong, totlong;
		
		pt_array = new double[n * dim];
		for(int i = 0; i < n; ++i) {
			for(int j = 0; j < dim; ++j) {
				pt_array[dim * i + j] = p[i][j];
			}
		}
		printf("%d, %d\n", dim, n);
		ismalloc = False;
		if(dim < 5) {
			sprintf(flags, "qhull d Qz Qt Qc Qbb Pg PF1e-20 i p s");// A-0.99999999 A0.99999999 C-0.000000001 C0.000000001 i p s");
		}
		else {
			sprintf(flags, "qhull d Qx Qz Qt Qbb Qc Pg PF1e-20 i p s");// A-0.99999999 A0.99999999 C-0.000000001 C0.000000001 i s");
		}
		outfile = fopen("tri.txt", "w");
		exitcode = qh_new_qhull(dim, n, pt_array, ismalloc, flags, outfile, errfile);
		qh_freeqhull(!qh_ALL);
		fclose(outfile);
		qh_memfreeshort (&curlong, &totlong);
		FILE* infile = fopen("tri.txt", "r");
		int num;
		int index;
		fscanf(infile, "%d", &num);
		triangles.resize(num);
		for(int i = 0; i < num; ++i) {
			triangles[i].resize(dim + 1);
			for(int j = 0; j < dim + 1; ++j) {
				fscanf(infile, "%d", &triangles[i][j]);
			}
		}
		fclose(infile);
		remove("tri.txt");
	}
}

double SimplexArea(vector<vector<double> > p) {
	int dim = p.size() - 1;
	double vol = 1.;
	char jobu = 'N';
	char jobvt = 'N';
	int M = dim + 1;
	int N = dim + 1;
	int lda = dim + 1;
	double* S = new double[dim + 1];
	double* U = new double[(dim + 1) * (dim + 1)];
	double* A = new double[(dim + 1) * (dim + 1)];
	int ldu = dim + 1;
	double* VT = new double[(dim + 1) * (dim + 1)];
	int ldvt = dim + 1;
	double* work = new double[5 * (dim + 1)];
	int lwork = 5 * (dim + 1);
	int info;
	for(int i = 0; i < p.size(); ++i) {
		for(int j = 0; j < dim; ++j) {
			A[i + j * (dim + 1)] = p[i][j];
		}
	}
	for(int i = 0; i < dim + 1; ++i) {
		A[i + dim * (dim + 1)] = 1.;
	}
	dgesvd_(&jobu, &jobvt, &M, &N, A, &lda, S, U, &ldu, VT, &ldvt, work, &lwork, &info);
	for(int i = 0; i < dim + 1; ++i) {
		vol *= S[i];
	}
	double coeff = 1.;
	for(int i = 1; i <= dim; ++i) {
		coeff *= i;
	}
	delete [] S;
	delete [] U;
	delete [] A;
	delete [] VT;
	delete [] work;
	return vol / coeff;
}

double Dist(const vector<double>& p1, const vector<double>& p2) {
	double SquareSum = 0;
	for(int i = 0; i < p1.size(); ++i) {
		SquareSum += pow(p1[i] - p2[i], 2);
	}
	return sqrt(SquareSum);
}

bool isWithin(const vector<double> p, const vector<vector<double> >& bounds) {
	bool flag = true;
	for(int i = 0; i < p.size(); ++i) {
		if(p[i] - bounds[i][0] < -1E-5 || p[i] - bounds[i][1] > 1E-5) {
			flag = false;
			break;
		}
	}
	return flag;
}

void ReadFile(string FileName, vector<vector<double> >& p, vector<vector<int> >& simpleces, vector<double>& volumes, vector<double>& simplexAreaList, vector<double>& ratios, string filename, string lamstr) {
	FILE* pf = fopen(FileName.c_str(), "r");
	int dim, num;
	fscanf(pf, "%d%d", &dim, &num);
	printf("%d, %d\n", dim, num);
	vector<vector<vector<double> > > PointsInCubes;
	vector<vector<vector<double> > > bounds;
	vector<double> VolInCubes;
	vector<double> RatioInCubes;
	PointsInCubes.resize(num);
	VolInCubes.resize(num);
	bounds.resize(num);
	RatioInCubes.resize(num);
	for(int i = 0; i < num; ++i) {
		bounds[i].resize(dim);
		for(int j = 0; j < dim; ++j) {
			bounds[i][j].resize(2);
			fscanf(pf, "%lf", &bounds[i][j][0]);
			fscanf(pf, "%lf", &bounds[i][j][1]);
		}
		fscanf(pf, "%lf", &VolInCubes[i]);
		unsigned x = UINT_MAX;
		x <<= dim;
		x = ~x;
		vector<double> tmpPoint;
		tmpPoint.resize(dim);
		for(unsigned ii = 0; ii <= x; ++ii) {
			unsigned tmp = ii;
			for(int jj = 0; jj < dim; ++jj) {
				tmpPoint[jj] = bounds[i][jj][tmp % 2];
				tmp >>= 1;
			}
			/*
			for(int jj = 0; jj < dim; ++jj) {
				printf("%20.10lf", tmpPoint[jj]);
			}
			printf("\n");
			*/
			PointsInCubes[i].push_back(tmpPoint);
		}
		/*
		double CubeVol = 1;
		vector<double> EdgeLen;
		for(int ii = 0; ii < bounds[i].size(); ++ii) {
			CubeVol *= (bounds[i][ii][1] - bounds[i][ii][0]);
			EdgeLen.push_back(bounds[i][ii][1] - bounds[i][ii][0]);
		}
		double Max = *max_element(EdgeLen.begin(), EdgeLen.end());
		double Min = *min_element(EdgeLen.begin(), EdgeLen.end());
		RatioInCubes[i] = Max / Min;
		if(CubeVol < pow(10, -1 * dim) && RatioInCubes[i] < 2) {
			vector<double> tmpPoint;
			tmpPoint.resize(dim);
			for(int ii = 0; ii < bounds[i].size(); ++ii) {
				tmpPoint[ii] = (bounds[i][ii][0] + bounds[i][ii][1]) / 2;
			}
			PointsInCubes[i].push_back(tmpPoint);
		}
		*/
	}
	/*
	double maxVol = *max_element(VolInCubes.begin(), VolInCubes.end());
	for(int i = 0; i < num; ++i) {
		double CubeVol = 1;
		vector<double> EdgeLen;
		for(int ii = 0; ii < bounds[i].size(); ++ii) {
			CubeVol *= (bounds[i][ii][1] - bounds[i][ii][0]);
			EdgeLen.push_back(bounds[i][ii][1] - bounds[i][ii][0]);
		}
		double Max = *max_element(EdgeLen.begin(), EdgeLen.end());
		double Min = *min_element(EdgeLen.begin(), EdgeLen.end());
		RatioInCubes[i] = Max / Min;
		if(VolInCubes[i] > 0.99 * maxVol && RatioInCubes[i] < 2.5) {
			vector<double> tmpPoint;
			tmpPoint.resize(dim);
			for(int ii = 0; ii < bounds[i].size(); ++ii) {
				tmpPoint[ii] = (bounds[i][ii][0] + bounds[i][ii][1]) / 2;
			}
			PointsInCubes[i].push_back(tmpPoint);
		}	
	}
	*/
	vector<vector<int> > LabelsInCubes;
	vector<vector<int> > SimplexLoc;
	LabelsInCubes.resize(num);
	SimplexLoc.resize(num);
	//printf("%d\n", PointsInCubes.size());
	for(int i = 0; i < PointsInCubes.size(); ++i) {
		//printf("%d\n", PointsInCubes[i].size());
		for(int j = 0; j < PointsInCubes[i].size(); ++j) {
			bool flag = true;
			//printf("%d\n", p.size());
			for(int k = 0; k < p.size(); ++k) {
				if(Dist(PointsInCubes[i][j], p[k]) < 1E-5) {
					flag = false;
					break;
				}
			}
			if(flag) {
				p.push_back(PointsInCubes[i][j]);
			}
		}
	}

	for(int i = 0; i < p.size(); ++i) {
		for(int j = 0; j < bounds.size(); ++j) {
			if(isWithin(p[i], bounds[j])) {
				LabelsInCubes[j].push_back(i);
			}
		}
	}
	/*
	for(int i = 0; i < LabelsInCubes.size(); ++i) {
		for(int j = 0; j < LabelsInCubes[i].size(); ++j) {
			printf("%d, ", LabelsInCubes[i][j]);
		}
		printf("\n");
	}
	*/
	int SimLocInd = 0;
	for(int i = 0; i < LabelsInCubes.size(); ++i) {
		vector<vector<double> > points;
		vector<vector<int> > tri;
		for(int j = 0; j < LabelsInCubes[i].size(); ++j) {
			points.push_back(p[LabelsInCubes[i][j]]);
		}
		delaunayn(points, tri);
		/*
		for(int j = 0; j < tri.size(); ++j) {
			for(int k = 0; k < tri[j].size(); ++k) {
				printf("%10d", tri[j][k]);
			}
			printf("\n");
		}
		*/
		//printf("%d\n", tri.size());
		for(int j = 0; j < tri.size(); ++j) {
			vector<int> simplexInd;
			simplexInd.resize(tri[j].size());
			for(int k = 0; k < tri[j].size(); ++k) {
				simplexInd[k] = LabelsInCubes[i][tri[j][k]];
			}
			simpleces.push_back(simplexInd);
			SimplexLoc[i].push_back(SimLocInd++);
			//printf("%d\n", points.size());
			vector<vector<double> > simplexNode;
			for(int k = 0; k < tri[j].size(); ++k) {
				//printf("%d\n", tri[j][k]);
				simplexNode.push_back(points[tri[j][k]]);
			}
			simplexAreaList.push_back(SimplexArea(simplexNode));
			volumes.push_back(VolInCubes[i]);
			ratios.push_back(RatioInCubes[i]);
		}
	}
	FILE* fp = fopen((string("output/simplecesLoc") + lamstr + filename).c_str(), "w");
	for(int i = 0; i < SimplexLoc.size(); ++i) {
		for(int j = 0; j < SimplexLoc[i].size(); ++j) {
			printf("%d, %d\n", i, SimplexLoc[i][j]);
			fprintf(fp, "%d\n", i);
		}
	}
	fclose(fp);
	/*
	for(int i = 0; i < p.size(); ++i) {
		printf("%d = ", i);
		for(int j = 0; j < p[i].size(); ++j) {
			printf("%20.10lf", p[i][j]);
		}
		printf("\n");
	}
	for(int i = 0; i < simpleces.size(); ++i) {
		for(int j = 0; j < simpleces[i].size(); ++j) {
			printf("%10d", simpleces[i][j]);
		} 
		printf("\n");
	}
	double s1 = 0, s2 = 0;
	for(int i = 0; i < simplexAreaList.size(); ++i) {
		s1 += simplexAreaList[i];
		s2 += simplexAreaList[i] * volumes[i];
	}
	
	for(int i = 0; i < ratios.size(); ++i) {
		printf("%lf\n", ratios[i]);
	}
	printf("volumes = \n");
	for(int i = 0; i < volumes.size(); ++i) {
		printf("%lf, ", volumes[i]);
	}
	*/
	//printf("%lf, %lf", s1, s2);
}

void test(vector<vector<double> >&p, vector<vector<int> >& tri) {
	vector<double> vol(tri.size(), 0.);
	for(int i = 0; i < tri.size(); ++i) {
		vector<vector<double> > Simplex;
		for(int j = 0; j < tri[i].size(); ++j) {
			Simplex.push_back(p[tri[i][j]]);
		}
		vol[i] = SimplexArea(Simplex);
	}
	double min = 1;
	for(int i = 0; i < vol.size(); ++i) {
		if(min > vol[i]) min = vol[i];
		printf("%20.18f", vol[i]);
	}
	printf("\n%20.18f%20.18f\n", accumulate(vol.begin(), vol.end(), 0.), min);
}

int main(int argc, char** argv) {
	vector<vector<double> > p;
	vector<vector<int> > simpleces;
	vector<double> volumes;
	vector<double> simplexAreaList;
	vector<double> ratios;
	string filename = string(argv[1]);
	double lambda;
	lambda = atof(argv[2]);
	char buffer[500];
	vector<double> vals;
	ReadFile(filename, p, simpleces, volumes, simplexAreaList, ratios, argv[1], argv[2]);
	QPSolve(p, simpleces, volumes, simplexAreaList, vals, lambda, ratios);
	sprintf(buffer, (string("output/points") + argv[2] + argv[1]).c_str());
	FILE* pf = fopen(buffer, "w");
	for(int i = 0; i < p.size(); ++i) {
		for(int j = 0; j < p[i].size(); ++j) {
			fprintf(pf, "%20.10lf", p[i][j]);
		}
		fprintf(pf, "\n");
	}
	fclose(pf);
	sprintf(buffer, (string("output/simpleces") + argv[2] + argv[1]).c_str());
	pf = fopen(buffer, "w");
	for(int i = 0; i < simpleces.size(); ++i) {
		for(int j = 0; j < simpleces[i].size(); ++j) {
			fprintf(pf, "%10d", simpleces[i][j]);
		}
		fprintf(pf, "\n");
	}
	fclose(pf);
	sprintf(buffer, (string("output/val") + argv[2] + argv[1]).c_str());
	pf = fopen(buffer, "w");
	for(int i = 0; i < vals.size(); ++i) {
		fprintf(pf, "%20.10lf\n", vals[i]);
	}
	fclose(pf);
	pf = fopen((string("output/count") + argv[2] + argv[1]).c_str(), "w");
	fprintf(pf, "%10d%10d%10d\n", simpleces.size(), vals.size(), p.size());
	fclose(pf);
	/*
	for(int ind = 0; ind < 15; ++ind) {
		vector<double> vals;
		lambda = pow(10, -1 * ind * 9. / 15. - 1);
		QPSolve(p, simpleces, volumes, simplexAreaList, vals, lambda, ratios);
		sprintf(buffer, "output/points%d.txt", ind);
		FILE* pf = fopen(buffer, "w");
		for(int i = 0; i < p.size(); ++i) {
			for(int j = 0; j < p[i].size(); ++j) {
				fprintf(pf, "%20.10lf", p[i][j]);
			}
			fprintf(pf, "\n");
		}
		fclose(pf);
		sprintf(buffer, "output/simpleces%d.txt", ind);
		pf = fopen(buffer, "w");
		for(int i = 0; i < simpleces.size(); ++i) {
			for(int j = 0; j < simpleces[i].size(); ++j) {
				fprintf(pf, "%10d", simpleces[i][j]);
			}
			fprintf(pf, "\n");
		}
		fclose(pf);
		sprintf(buffer, "output/val%d.txt", ind);
		pf = fopen(buffer, "w");
		for(int i = 0; i < vals.size(); ++i) {
			fprintf(pf, "%20.10lf\n", vals[i]);
		}
		fclose(pf);
	}
	*/
	/*
	vector<vector<double> > mat;
	mat.resize(1);
	for(int i = 0; i < mat.size(); ++i) {
		mat[i].resize(mat.size());
	}
	mat[0][0] = 5;
	
	FILE* fp = fopen("Mymatrix.txt", "r");
	for(int i = 0; i < mat.size(); ++i) {
		for(int j = 0; j < mat.size(); ++j) {
			fscanf(fp, "%lf", &mat[i][j]);
			printf("%5.2lf", mat[i][j]);
		}
		printf("\n");
	}
	
	Determinant(mat);
	*/
	return 0;
}
