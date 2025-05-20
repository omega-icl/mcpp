PATH_MC = $(shell cd ../../ ; pwd)

# Compilation options
include $(PATH_MC)/src/makeoptions.mk

#####

test_MLP : test_MLP.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

test_MLPREG : test_MLPREG.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

%.o: %.cpp 
	$(CPP) -c $(FLAG_CPP) $(FLAG_MC) $(INC_MC) $< -o $@

#####

graph : gendot gensvg

gendot :
	dot -Tpng -O $(shell ls *.dot)

gensvg :
	@for FILE in $(shell ls *.png); do \
	  echo $$FILE; \
	  convert -resize 70% $$FILE $$FILE; \
	  (display $$FILE &); \
	done

#####

clean :
	rm -f test_MLP test_MLPREG
	-(rm -f *.o *.out *.eps *.dot *.png *.lp *~)

