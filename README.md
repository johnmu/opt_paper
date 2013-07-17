opt_paper
=========
The code needs QHULL Cplex(IBM) and Lapack(Fortran Lib) to compile.
------------------
make
------------------
./TRI partition-file lambda<br />
------------------
./HELL partition-file output/val-lambda-partition-file output/simpleces-lambda-partition-file output/points-lambda-partition-file output/count-lambda-partition-file output/simplecesLoc-lambda-partition-file N<br />
------------------
./SAMPLE partition-file output/val-lambda-partition-file output/simpleces-lambda-partition-file output/points-lambda-partition-file output/count-lambda-partition-file output/simplecesLoc-lambda-partition-file<br />
------------------

TRI will produce the files needed by HELL and SAMPLE in folder output/<br />
N is the number of samples used to calculate the helliger distance.
