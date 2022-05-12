###############################################################################
##                                 Makefile                                  ##
###############################################################################

TARGET = A.OUT

###############################################################################

include $(HOME)/local/cpplapack/makefiles/Makefile.g++
#include $(HOME)/local/cpplapack/makefiles/Makefile.icpc

###############################################################################

THIS_FILE:= $(firstword $(MAKEFILE_LIST))
OPTIONS:= $(MAKEFLAGS)

HEADERS:= $(shell find ./ -iname '*.hpp')
SOURCES:= $(shell find ./ -iname '*.cpp')
SOURCES:= $(sort $(SOURCES))

OBJECTS:= $(SOURCES:%.cpp=%.o)

ifeq ($(OMP_NUM_THREADS),)
	OMP_NUM_THREADS:= $(shell cat /proc/cpuinfo | grep processor | wc -l)
endif

###############################################################################
###############################################################################
###############################################################################

debug: depend
	@echo -e "################################"
	@echo -e "#### building in debug mode ####"
	@echo -e "################################\n"
	@$(MAKE) -j $(OMP_NUM_THREADS) -f $(THIS_FILE) --no-print-directory _build DEBUG=1
	@echo -e "#############################"
	@echo -e "#### built in debug mode ####"
	@echo -e "#############################\n"

release: depend
	@echo -e "##################################"
	@echo -e "#### building in release mode ####"
	@echo -e "##################################\n"
	@$(MAKE) -j $(OMP_NUM_THREADS) -f $(THIS_FILE) --no-print-directory _build RELEASE=1
	@echo -e "###############################"
	@echo -e "#### built in release mode ####"
	@echo -e "###############################\n"

profile: depend
	@echo -e "##################################"
	@echo -e "#### building in profile mode ####"
	@echo -e "##################################\n"
	@$(MAKE) -j $(OMP_NUM_THREADS) -f $(THIS_FILE) --no-print-directory _build PROFILE=1
	@echo -e "###############################"
	@echo -e "#### built in profile mode ####"
	@echo -e "###############################\n"

verbose: depend
	@echo -e "##################################"
	@echo -e "#### building in verbose mode ####"
	@echo -e "##################################\n"
	@$(MAKE) -j $(OMP_NUM_THREADS) -f $(THIS_FILE) --no-print-directory _build DEBUG=1 VERBOSE=1
	@echo -e "###############################"
	@echo -e "#### built in verbose mode ####"
	@echo -e "###############################\n"

###############################################################################

.SUFFIXES: .cpp .o
.cpp.o:
	$(CXX) -c $< $(FLAGS) $(CFLAGS) $(INCLUDE_DIRS) $(MACROS) -o $@
	@echo

_build: $(OBJECTS)
	$(CXX) $(OBJECTS) $(LIBS) $(FLAGS) $(LFLAGS) $(LIB_DIRS) -o $(TARGET)
	@echo

###############################################################################

depend:
	makedepend -f- -Y $(SOURCES) > Makefile.depend 2> /dev/null
#	gccmakedep -- -I./ -MM -- $(SOURCES)
#	$(CXX) -MM -I./ $(SOURCES) > Makefile.depend
#	$(CXX) -MM $(INCLUDE_DIRS) $(SOURCES) > Makefile.depend
	@echo

clean:
	rm -f $(OBJECTS)

fullclean:
	rm -f $(shell find -name '*.o') Makefile.depend std err *~ $(TARGET)

###############################################################################
-include Makefile.depend
