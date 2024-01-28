# THIRD-PARTY LIBRARIES -- CHANGE AS APPROPRIATE

PATH_PROFIL  = /home/bchachua/Programs/ThirdParty/Profil-2.0.8
#PATH_PROFIL  = /opt/Profil-2.0.8
LIB_PROFIL = -L$(PATH_PROFIL)/lib -lProfilPackages -lProfil -lBias -llr
INC_PROFIL = -I$(PATH_PROFIL)/include

PATH_FILIB = /opt/filib++
LIB_FILIB = -L$(PATH_FILIB)/lib -lprim
INC_FILIB = -I$(PATH_FILIB)/include/ -I$(PATH_FILIB)/include/interval
FLAGS_FILIB = -frounding-math -ffloat-store

PATH_MCPP = /home/bchachua/Programs/Devel/MC++/MC++_1.0
INC_MCPP  = -I$(PATH_MCPP)/include

PATH_LAPACK = /opt/cpplapack-2010.03.27/
LIB_LAPACK = -llapack -lblas
INC_LAPACK = -I$(PATH_LAPACK)/include

PATH_ODEBND = /home/bchachua/Programs/Devel/MC++/ODEBND_1.1
INC_ODEBND  = -I$(PATH_ODEBND)/src

PATH_CPLEX = /opt/ibm/ILOG/CPLEX_Studio1251/cplex
PATH_CONCERT = /opt/ibm/ILOG/CPLEX_Studio1251/concert
LIB_CPLEX = -L$(PATH_CPLEX)/lib/x86-64_sles10_4.1/static_pic -lilocplex -lcplex \
            -L$(PATH_CONCERT)/lib/x86-64_sles10_4.1/static_pic -lconcert \
	    -lm -pthread
INC_CPLEX = -I$(PATH_CPLEX)/include -I$(PATH_CONCERT)/include
FLAGS_CPLEX = -m64 -fPIC -fexceptions -DIL_STD
#FLAGS_CPLEX = -m64 -fPIC -fexceptions -DNDEBUG -DIL_STD

PATH_GUROBI = $(GUROBI_HOME)
LIB_GUROBI = -L$(PATH_GUROBI)/lib -lgurobi_c++ -lgurobi56 -pthread
INC_GUROBI = -I$(PATH_GUROBI)/include

PATH_IPOPT = /opt/Ipopt-3.10.2
LIB_IPOPT = -L$(PATH_IPOPT)/lib -lipopt -lcoinhsl -llapack -lblas -ldl
INC_IPOPT = -I$(PATH_IPOPT)/include

PATH_SUNDIALS = /opt/sundials-2.5.0
LIB_CVODES = -L$(PATH_SUNDIALS)/lib -lsundials_cvodes -lsundials_nvecserial \
             -llapack -lblas
INC_CVODES = -I$(PATH_SUNDIALS)/include -I$(PATH_SUNDIALS)/include/sundials \
             -I$(PATH_SUNDIALS)/include/cvodes -I$(PATH_SUNDIALS)/include/nvector

LIB_GSL = -lgsl -lgslcblas -lm
INC_GSL = 

# COMPILATION -- CHANGE AS APPROPRIATE

# PROF = -pg
OPTIM = -O2
#OPTIM = 
DEBUG = -g
WARN  = -Wall

CC = gcc
CPP = g++
FLAGS_CPP = $(DEBUG) $(OPTIM) $(WARN) $(FLAGS_FILIB) $(FLAGS_CPLEX)

archdyn = ld
archdflags = -r

archive  = ar
archflags = -rs

LINK = $(CPP)
FLAGS_LINK = 