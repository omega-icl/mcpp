# THIRD-PARTY LIBRARIES <<-- CHANGE AS APPROPRIATE -->>

PATH_MC_MK := $(dir $(lastword $(MAKEFILE_LIST)))

PATH_MC        = $(abspath $(PATH_MC_MK)../)
PATH_3RD_PARTY = $(PATH_MC)/src/3rdparty
PATH_EXTERN    = $(PATH_MC)/extern

PATH_LAPACK = #$(PATH_3RD_PARTY)/cpplapack-2015.05.11-1
LIB_LAPACK  = -larmadillo -llapack -lblas
INC_LAPACK  = #-I$(PATH_LAPACK)/include
FLAG_LAPACK = -DMC__USE_ARMADILLO

PATH_FADBAD = #$(PATH_3RD_PARTY)/fadbad++
LIB_FADBAD  =
INC_FADBAD  = #-I$(PATH_FADBAD)
FLAG_FADBAD = #-DMC__USE_FADBAD #-DMC__DAG_FADIFF -DMC__DAG_BADIFF #-DMC__USE_FADBAD 

PATH_PROFIL = $(PROFIL_HOME)
LIB_PROFIL  = -L$(PATH_PROFIL)/lib -lProfilPackages -lProfil -lBias -llr
INC_PROFIL  = -I$(PATH_PROFIL)/include
FLAG_PROFIL = #-DMC__USE_PROFIL

PATH_FILIB  = $(FILIB_HOME)
LIB_FILIB   = -L$(PATH_FILIB)/lib -lprim
INC_FILIB   = -I$(PATH_FILIB)/include -I$(PATH_FILIB)/include/interval
FLAG_FILIB = -frounding-math #-DMC__USE_FILIB

PATH_BOOST = #$(PATH_3RD_PARTY)/boost
LIB_BOOST  =
INC_BOOST  = -I/usr/include #-I$(PATH_BOOST)
#FLAG_BOOST = -DBOOST_UBLAS_NO_STD_CERR -DMC__USE_BOOST
FLAG_BOOST = -DMC__USE_BOOST

PATH_HSL =
LIB_HSL  = -lmc13 -lmc21 -lmc33 -lgfortran
INC_HSL  =
FLAG_HSL = -DMC__USE_HSL

PATH_TORCH = $(LIBTORCH_HOME)
LIB_TORCH  = -L$(PATH_TORCH)/lib -ltorch_cpu -lc10
#-Wl,-rpath,$(LIBTORCH_HOME)/lib
INC_TORCH  = -I$(PATH_TORCH)/include -I$(PATH_TORCH)/include/torch/csrc/api/include
FLAG_TORCH = -DMC__USE_TORCH

INC_PYTHON = $(shell python3 -c "from sysconfig import get_paths; print(get_paths()['include'])")
INC_PYBIND11 = -I$(INC_PYTHON) -I$(PATH_EXTERN)/pybind11/include
#LIB_CPPUNIT = -lcppunit

# COMPILATION <<-- CHANGE AS APPROPRIATE -->>

PROF = #-pg
OPTIM = #-O2
DEBUG = -g
WARN  = -Wall -Wno-misleading-indentation -Wno-unknown-pragmas -Wno-parentheses -Wno-return-type
CPP17 = -std=c++17
CC    = gcc
CPP   = g++
# CPP   = icpc

# <<-- NO CHANGE BEYOND THIS POINT -->>

FLAG_CPP  = $(DEBUG) $(OPTIM) $(CPP17) $(WARN) $(PROF)
LINK      = $(CPP)
FLAG_LINK = $(PROF)

FLAG_MC  = -DMC__USE_THREAD $(FLAG_FADBAD) $(FLAG_LAPACK) $(FLAG_BOOST)
LIB_MC   = $(LIB_FADBAD) $(LIB_LAPACK) $(LIB_BOOST)
INC_MC   = -I$(PATH_MC)/src/mc $(INC_FADBAD) $(INC_LAPACK) $(INC_BOOST)

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

ifneq (,$(findstring -DMC__WITH_SPQR, $(FLAG_SUITESPARSE)))
 FLAG_MC += $(FLAG_SUITESPARSE)
 INC_MC  += $(INC_SUITESPARSE)
 LIB_MC  += $(LIB_SUITESPARSE)
endif

ifneq (,$(findstring -DMC__USE_TORCH, $(FLAG_TORCH)))
 FLAG_MC += $(FLAG_TORCH)
 INC_MC  += $(INC_TORCH)
 LIB_MC  += $(LIB_TORCH)
endif
