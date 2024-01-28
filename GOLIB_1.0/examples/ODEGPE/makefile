# Executable name
EXE = ODEGPE
GPEXE = ODEGPE_GP

# Compilation options
PATH_GO = $(shell cd ../../ ; pwd)
include $(PATH_GO)/src/makeoptions.mk

LIB_ALL = $(LIB_IPOPT) $(LIB_CVODES) $(LIB_GSL) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_GUROBI) $(LIB_CPLEX) $(LIB_LAPACK)
INC_ALL = -I$(PATH_GO)/include $(INC_MCPP) $(INC_IPOPT) $(INC_CVODES) $(INC_GSL) $(INC_PROFIL) $(INC_FILIB) $(INC_GUROBI) $(INC_CPLEX) $(INC_ODEBND) $(INC_LAPACK)
#OBJS = test_tmodel2.o
OBJS = main.o
GPOBJS = gnuplot.o

#####

$(EXE) : $(OBJS)
	$(LINK) $(PROF) $(FLAGS_LINK) -o $(EXE) $(OBJS) $(LIB_ALL)
       
main.o: main.cpp main.hpp
	$(CPP) -c $(PROF) $(FLAGS_CPP) $(INC_ALL) -o main.o main.cpp
       
test_tmodel2.o: test_tmodel2.cpp main.hpp
	$(CPP) -c $(PROF) $(FLAGS_CPP) $(INC_ALL) -o test_tmodel2.o test_tmodel2.cpp

#####

GNUPLOT : $(GPOBJS)
	$(LINK) $(PROF) $(FLAGS_LINK) -o $(GPEXE) $(GPOBJS)
       
gnuplot.o: gnuplot.cpp
	$(CPP) -c $(PROF) $(FLAGS_CPP) -o gnuplot.o gnuplot.cpp

#####

clean :
	rm -f $(EXE) $(OBJS) $(GPEXE) $(GPOBJS)
