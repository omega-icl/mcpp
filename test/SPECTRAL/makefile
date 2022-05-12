PATH_MC = $(shell cd ../../ ; pwd)

# Compilation options
include $(PATH_MC)/src/makeoptions.mk

#####

all: SB-2D SB-C

SB-2D : SB-2D.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

SB-C : SB-C.o
	$(LINK) $(FLAGS_LINK) -o $@ $< $(LIB_MC)

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
	rm -f SB-2D SB-C *.o *.out *.eps *.dot *.svg *~

