# Executable name
EXE = MCIPOPT

# Compilation options
PATH_GO = $(shell cd ../../ ; pwd)
include $(PATH_GO)/src/makeoptions.mk

LIB_ALL = $(LIB_IPOPT)
INC_ALL = -I$(PATH_GO)/include $(INC_MCPP) $(INC_IPOPT)
OBJS = main.o

#####

$(EXE) : $(OBJS)
	$(LINK) $(PROF) $(FLAGS_LINK) -o $(EXE) $(OBJS) $(LIB_ALL)
       
main.o: main.cpp
	$(CPP) -c $(PROF) $(FLAGS_CPP) $(INC_ALL) -o main.o main.cpp

#####

clean :
	rm -f $(EXE) $(OBJS)
