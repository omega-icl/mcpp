PATH_MC = $(shell cd ../../ ; pwd)

# Compilation options
include $(PATH_MC)/src/makeoptions.mk

#####

all: test_fadbad test_fdiff test_bdiff test_tdiff test_vs_fadbad test_perf test_forms test_mc_dep

test_fadbad : test_fadbad.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

test_fdiff : test_fdiff.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

test_bdiff : test_bdiff.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

test_tdiff : test_tdiff.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

test_vs_fadbad : test_vs_fadbad.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

test_perf : test_perf.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

test_forms : test_forms.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

test_mc_dep : test_mc_dep.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

%.o: %.cpp
	$(CPP) -c $(FLAG_CPP) $(FLAG_MC) $(INC_MC) $< -o $@

#####

clean :
	rm -f test_fadbad test_fdiff test_bdiff test_tdiff test_vs_fadbad test_perf test_forms test_mc_dep
	rm -f *.o
