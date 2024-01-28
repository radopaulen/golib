# Executable name
EXE = fftest

PATH_MC  = $(shell cd ../../ ; pwd)
LIB_MC   = -L$(PATH_MC)/lib
INC_MC   = -I$(PATH_MC)/include
FILE_DOT = $(shell ls *.dot)
FILE_SVG = $(shell ls *.svg)

OBJS     = main.o
#OBJS     = test_EA.o

# Compilation options
include $(PATH_MC)/src/makeoptions.mk

#####

$(EXE) : $(OBJS)
	$(LINK) $(PROF) $(FLAGS_LINK) -o $(EXE) $(OBJS)

main.o: main.cpp ffdep.hpp ffunc.hpp
	$(CPP) -c $(PROF) $(FLAGS_CPP) $(INC_MC) -o main.o main.cpp

test_EA.o: test_EA.cpp ffdep.hpp ffunc.hpp EA.hpp
	$(CPP) -c $(PROF) $(FLAGS_CPP) $(INC_MC) -o test_EA.o test_EA.cpp

#####

graph :
	dot -Tsvg -O $(FILE_DOT)
	@for FILE in $(FILE_SVG); do \
	  echo $$FILE; \
	  (display $$FILE &); \
	done

#####

clean :
	rm -f $(EXE) $(OBJS) *.dot *.svg
