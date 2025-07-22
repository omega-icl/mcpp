PATH_MC = $(shell cd ../../ ; pwd)

# Compilation options
include $(PATH_MC)/src/makeoptions.mk

PATH_CANON = /home/bchachua/Programs/bitbucket/canon40/src
INC_CANON  = -I$(PATH_CANON)

PATH_GAMS  = /opt/gams/gams50.1_linux_x64_64_sfx
LIB_GAMS   =
INC_GAMS   = -I$(PATH_GAMS)/apifiles/C/api
FLAG_GAMS  = -DMC__WITH_GAMS=\"$(PATH_GAMS)\"

#####

SM-2D : SM-2D.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

SM-HD : SM-HD.o gmomcc.o gevmcc.o optcc.o palmcc.o
	$(LINK) $(FLAG_LINK) -o $@ $^ $(LIB_MC)

testpwlu : testpwlu.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

testpwcu : testpwcu.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

%.o: %.cpp 
	$(CPP) -c $(FLAG_CPP) $(FLAG_MC) $(FLAG_GAMS) $(INC_MC) $(INC_GUROBI) $(INC_GAMS) $(INC_CANON) $< -o $@

%.o : %.c
	$(CPP) -c $(FLAG_CPP) $(INC_GAMS) $< -o $@

%.c : $(PATH_GAMS)/apifiles/C/api/%.c
	cp $< $@

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
	rm -f SM-2D SM-HD testpwlu testpwcu *.o *.out *.eps *.dot *.svg *~

