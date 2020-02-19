# THIRD-PARTY LIBRARIES <<-- CHANGE AS APPROPRIATE -->>

PATH_3RD_PARTY = $(shell cd ../../src/ ; pwd )/3rdparty

PATH_PROFIL  = $(HOME)/Programs/ThirdParty/Profil-2.0.8
LIB_PROFIL = -L$(PATH_PROFIL)/lib -lProfilPackages -lProfil -lBias -llr
INC_PROFIL = -I$(PATH_PROFIL)/include

PATH_BOOST = 
LIB_BOOST =
INC_BOOST = 

#PATH_FILIB = /opt/filib++
#LIB_FILIB = -L$(PATH_FILIB)/lib -lprim
#INC_FILIB = -I$(PATH_FILIB)/include/
#FLAGS_FILIB = -frounding-math -ffloat-store

PATH_LAPACK = $(PATH_3RD_PARTY)/cpplapack-2015.05.11-1
LIB_LAPACK = -llapack -lblas
INC_LAPACK = -I$(PATH_LAPACK)/include

PATH_FADBAD = $(PATH_3RD_PARTY)/fadbad++
LIB_FADBAD = 
INC_FADBAD = -I$(PATH_FADBAD)

#LIB_SDPA   = -lsdpa -ldmumps_seq
#INC_SDPA   = #-I/usr/include/

PATH_HSL = 
LIB_HSL = -lmc13 -lmc21 -lmc33 -lgfortran
INC_HSL = 
FLAGS_HSL = -DMC__USE_HSL

LIB_CPPUNIT = -lcppunit

# COMPILATION <<-- CHANGE AS APPROPRIATE -->>

# PROF = -pg
#OPTIM = -Ofast
DEBUG = -g
WARN  = -Wall -Wno-unknown-pragmas
CPP17 = -std=c++1z

CC = gcc-8
CPP = g++-8
# CPP = icpc
FLAGS_CPP = $(DEBUG) $(OPTIM) $(CPP17) $(WARN) $(FLAGS_FILIB)

LINK = $(CPP)
FLAGS_LINK = 
