CPLEX_INCLUDE = /home/michaelyoung/ILOG/CPLEX_Studio_Academic123/cplex/include /home/michaelyoung/ILOG/CPLEX_Studio_Academic123/concert/include
CPLEX_LIBS = /home/michaelyoung/ILOG/CPLEX_Studio_Academic123/concert/lib/x86_sles10_4.1/static_pic /home/michaelyoung/ILOG/CPLEX_Studio_Academic123/cplex/lib/x86_sles10_4.1/static_pic
CPLEXLIBS = ilocplex cplex concert
LDFLAGS = -pthread $(addprefix -L, $(CPLEX_LIBS)) $(addprefix -l, $(CPLEXLIBS))

CXX = g++
DESTDIR = /home/michaelyoung/qhull-qhull
INCDIR = $(DESTDIR)/src
LIBDIR = $(DESTDIR)/lib
OBJ = /home/michaelyoung/qhull-qhull/src/libqhullstatic

TRI: TRIAN QPSOLVE HELLINGER SAMPLING
	$(CXX) Triangulation.o QPSolve.o -I$(INCDIR)/libqhull/ -L$(LIBDIR) -lqhull6 -llapack -lblas -lm $(addprefix -I, $(CPLEX_INCLUDE)) $(LDFLAGS) -o TRI
TRIAN: Triangulation.cpp
	$(CXX) -c Triangulation.cpp -I$(INCDIR)/libqhull/ -L$(LIBDIR) -lqhull6 -llapack -lblas -lm
QPSOLVE: QPSolve.cpp
	$(CXX) -c QPSolve.cpp $(addprefix -I, $(CPLEX_INCLUDE)) $(LDFLAGS) -llapack -lblas -lm
HELLINGER: Hellinger.cpp
	$(CXX) mtrand.cpp Hellinger.cpp -llapack -lblas -lm -o HELL
SAMPLING: Sampling.cpp
	$(CXX) mtrand.cpp Sampling.cpp -llapack -lblas -lm -o SAMPLE
.PHONY: clean
clean: 
	rm -f *.o *~
