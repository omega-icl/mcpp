# THIRD-PARTY LIBRARIES <<-- CHANGE AS APPROPRIATE -->>

PATH_3RD_PARTY = $(HOME)/Programs/bitbucket/mcpp23/src/3rdparty
PATH_EXTERN = $(HOME)/Programs/bitbucket/mcpp23/extern

PATH_LAPACK = $(PATH_3RD_PARTY)/cpplapack-2015.05.11-1
LIB_LAPACK  = -llapack -lblas
INC_LAPACK  = -I$(PATH_LAPACK)/include
FLAG_LAPACK =

PATH_FADBAD = $(PATH_3RD_PARTY)/fadbad++
LIB_FADBAD  = 
INC_FADBAD  = -I$(PATH_FADBAD)
FLAG_FADBAD =

PATH_PROFIL = /opt/Profil-2.0.8
LIB_PROFIL  = -L$(PATH_PROFIL)/lib -lProfilPackages -lProfil -lBias -llr
INC_PROFIL  = -I$(PATH_PROFIL)/include
FLAG_PROFIL = #-DMC__USE_PROFIL

PATH_FILIB  = /opt/filib++
LIB_FILIB   = -L$(PATH_FILIB)/lib -lprim
INC_FILIB   = -I$(PATH_FILIB)/include -I$(PATH_FILIB)/include/interval
FLAG_FILIB = -frounding-math #-DMC__USE_FILIB

PATH_BOOST = $(PATH_3RD_PARTY)/boost
LIB_BOOST  = 
INC_BOOST  = -I/usr/include -I$(PATH_BOOST) 
FLAG_BOOST = -DBOOST_UBLAS_NO_STD_CERR -DMC__USE_BOOST

PATH_HSL = 
LIB_HSL  = -lmc13 -lmc21 -lmc33 -lgfortran
INC_HSL  = 
FLAG_HSL = -DMC__USE_HSL

PATH_BOOSTPYTHON = /opt/boost_1_83_0
INC_BOOSTPYTHON = -I/usr/include/python3.10 -I$(PATH_BOOSTPYTHON)/include
LIB_BOOSTPYTHON = -lpython3.10 -L$(PATH_BOOSTPYTHON)/lib -lboost_python310

INC_PYBIND11 = -I/usr/include/python3.10 -I$(PATH_EXTERN)/pybind11/include

LIB_CPPUNIT = -lcppunit

# COMPILATION <<-- CHANGE AS APPROPRIATE -->>

# PROF = -pg
OPTIM = #-O2
DEBUG = -g
WARN  = -Wall -Wno-misleading-indentation -Wno-unknown-pragmas -Wno-parentheses
CPP17 = -std=c++17
CC    = gcc-12
CPP   = g++-12
# CPP   = icpc

# <<-- NO CHANGE BEYOND THIS POINT -->>

FLAG_CPP  = $(DEBUG) $(OPTIM) $(CPP17) $(WARN) $(PROF)
LINK      = $(CPP)
FLAG_LINK = $(PROF)

FLAG_MC  = $(FLAG_FADBAD) $(FLAG_LAPACK) $(FLAG_BOOST)
LIB_MC   = $(LIB_FADBAD) $(LIB_LAPACK) $(LIB_BOOST)
INC_MC   = -I$(HOME)/Programs/bitbucket/mcpp23/src/mc $(INC_FADBAD) $(INC_LAPACK) $(INC_BOOST)

ifneq (,$(findstring -DMC__USE_HSL, $(FLAG_HSL)))
 FLAG_MC += $(FLAG_HSL)
 INC_MC  += $(INC_HSL)
 LIB_MC  += $(LIB_HSL)
endif

ifneq (,$(findstring -DMC__USE_PROFIL, $(FLAG_PROFIL)))
 FLAG_MC += $(FLAG_PROFIL)
 INC_MC  += $(INC_PROFIL)
 LIB_MC  += $(LIB_PROFIL)
endif

ifneq (,$(findstring -DMC__USE_FILIB, $(FLAG_FILIB)))
 FLAG_MC += $(FLAG_FILIB)
 INC_MC  += $(INC_FILIB)
 LIB_MC  += $(LIB_FILIB)
endif


