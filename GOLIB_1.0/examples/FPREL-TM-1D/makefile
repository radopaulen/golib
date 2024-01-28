# Executable name
EXE = FPREL-1D

PATH_MC = $(shell cd ../../ ; pwd)
LIB_MC = -L$(PATH_MC)/lib $(LIB_PROFIL) $(LIB_FILIB) $(LIB_GUROBI) $(LIB_CPLEX) -llapack -lblas
INC_MC = -I$(PATH_MC)/include $(INC_PROFIL) $(INC_FILIB) $(INC_GUROBI) $(INC_CPLEX)
OBJS = main.o

# Compilation options
include $(PATH_MC)/src/makeoptions.mk

#####

$(EXE) : $(OBJS)
	$(LINK) $(PROF) $(FLAGS_LINK) -o $(EXE) $(OBJS) $(LIB_MC)
       
main.o: main.cpp
	$(CPP) -c $(PROF) $(FLAGS_CPP) $(INC_MC) -o main.o main.cpp

#####

MC :
	(cd $(PATH_MC)/src; make; make install)

#####

clean :
	rm -f $(EXE) $(OBJS)
