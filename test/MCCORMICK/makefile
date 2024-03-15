PATH_MC = $(shell cd ../../ ; pwd)

# Compilation options
include $(PATH_MC)/src/makeoptions.mk

#####

all: MC-1D MC-2D

MC-1D : MC-1D.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

MC-2D : MC-2D.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

%.o: %.cpp 
	$(CPP) -c $(FLAG_CPP) $(FLAG_MC) $(INC_MC) $< -o $@

#####

graph : gendot gensvg

gendot :
	dot -Tsvg -O $(shell ls *.dot)

gensvg :
	@for FILE in $(shell ls *.svg); do \
	  echo $$FILE; \
	  (display $$FILE &); \
	done

#####

clean :
	rm -f MC-1D MC-2D *.o *.out *.eps *.dot *.svg *~

