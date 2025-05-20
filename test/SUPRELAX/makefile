PATH_MC = $(shell cd ../../ ; pwd)

# Compilation options
include $(PATH_MC)/src/makeoptions.mk

#####

SM-2D : SM-2D.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

testpwlu : testpwlu.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

testpwcu : testpwcu.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

%.o: %.cpp 
	$(CPP) -c $(FLAG_CPP) $(FLAG_MC) $(INC_MC) $(INC_GUROBI) $< -o $@

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
	rm -f SM-2D testpwlu testpwcu *.o *.out *.eps *.dot *.svg *~

