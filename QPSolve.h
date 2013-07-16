#ifndef _QPSOLVE_
#define _QPSOLVE_
#include <vector>
double Determinant(std::vector<std::vector<double> >&);
void QPSolve(std::vector<std::vector<double> >&, std::vector<std::vector<int> >&, std::vector<double>&, std::vector<double>&, std::vector<double>&, double, std::vector<double>&);
#endif
