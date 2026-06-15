PATH_MC = $(shell cd ../../ ; pwd)

# Compilation options
include $(PATH_MC)/src/makeoptions.mk

#####

all: SPOLY-TEST

SPOLY-TEST : SPOLY-TEST.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

%.o: %.cpp 
	$(CPP) -c $(FLAG_CPP) $(FLAG_MC) $(INC_MC) $< -o $@

#####

clean :
	rm -f SPOLY-TEST
	-(rm -f *.o *.out *.eps *.dot *.svg *.lp *~)

