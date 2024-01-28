# Executable name
target=SBB

rootdir     = $(shell cd ../../ ; pwd)
libpath     = $(rootdir)/lib
includepath = $(rootdir)/include
profilpath  = /home/bchachua/Programs/other/Profil-2.0.8
filibpath   = /home/bchachua/Programs/other/filibsrc
ipoptpath   = /home/bchachua/Programs/other/Ipopt-3.7.0

SYSTEM     = x86_debian4.0_4.1
LIBFORMAT  = static_pic
CPLEXLIBDIR   = $(CPLEX_HOME)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERT_HOME)/lib/$(SYSTEM)/$(LIBFORMAT)

# Compilation options
include $(rootdir)/src/makeoptions.mk

libobjs = main.o
liblink = -L$(profilpath)/lib -lProfilPackages -lProfil -lBias -llr \
          -L$(filibpath)/lib -lprim -L$(ipoptpath)/lib -lipopt \
          -L$(GUROBI_HOME)/lib -lgurobi_c++ -lgurobi40 \
          -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert \
          -m32 -lm -pthread -llapack -lblas
#-L$(libpath) -lglobal \
          
inclink = -I$(includepath) -I$(profilpath)/include -I$(filibpath)/interval \
          -I$(ipoptpath)/include -I$(GUROBI_HOME)/include \
          -I$(CPLEX_HOME)/include -I$(CONCERT_HOME)/include

#####

$(target) : $(libobjs)
	$(link) $(prof) $(linkflags) -o $(target) $(libobjs) $(liblink)

main.o: main.cpp
	$(cpp) -c $(prof) $(cppflags) $(inclink) -o main.o main.cpp

#####

MC :
	(cd $(rootdir)/src; make; make install)

#####

clean :
	rm -f $(target) $(libobjs)
