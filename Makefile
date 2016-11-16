#------------------------------------------------------------------------
#   Copyright (C) Inria 2016
#
#------------------------------------------------------------------------
#   Project:    Neotrope
#   Created on: Feb 29, 2016
#   Author:     Vincent KUBICKI <vincent.kubicki@inria.fr>
#------------------------------------------------------------------------

# generic build variables
UNITY_ROOT = Unity-master
CFLAGS += -g -pthread -MD -MP -DTEST
CFLAGS += -Isrc -I$(UNITY_ROOT)/src
# CFLAGS += -DLINDEBUG

# library files
LIB_LINALG = lib/libLinAlg.a
SRC_LINALG += $(wildcard src/*.c)
O_LINALG = $(SRC_LINALG:%.c=%.o)

# unit-test files
EXE_UTEST = utest/utest

SRC_UTEST += $(wildcard utest/libLinAlg/*.c)
SRC_UTEST += utest/UTest.c
O_UTEST = $(SRC_UTEST:%.c=%.o)

SRC_UYTEST = $(UNITY_ROOT)/src/unity.c
O_UYTEST = $(SRC_UYTEST:%.c=%.o)

# make recipes
all: lib utest

clean: cleanBin cleanData cleanSrc cleanUTest

cleanData:
	rm -f data/smallOutput.csv

cleanBin:
	rm -f $(EXE_UTEST)
	rm -f $(LIB_LINALG)

cleanSrc: # only cleans sources
	find -L src -type f -name "*.o" -exec rm -f {} \;
	find -L src -type f -name "*.d" -exec rm -f {} \;
	find -L src -type f -name "*.a" -exec rm -f {} \;
	find -L src -type f -name "*.so" -exec rm -f {} \;
	find -L src -type f -name "*.exe" -exec rm -f {} \;
	find -L src -type f -name "*.dll" -exec rm -f {} \;

cleanUTest: # only cleans unit test
	find -L utest -type f -name "*.o" -exec rm -f {} \;
	find -L utest -type f -name "*.d" -exec rm -f {} \;
	find -L utest -type f -name "*.a" -exec rm -f {} \;
	find -L utest -type f -name "*.so" -exec rm -f {} \;
	find -L utest -type f -name "*.exe" -exec rm -f {} \;
	find -L utest -type f -name "*.dll" -exec rm -f {} \;

lib: $(LIB_LINALG)

$(LIB_LINALG): $(O_LINALG)
	$(AR) $(ARFLAGS) $@ $^

utest: $(EXE_UTEST)

$(EXE_UTEST): $(LIB_LINALG) $(O_UYTEST) $(O_UTEST)
	$(CXX) $^ -lpthread -o $@

-include $(SRC_LINALG:%.c=%.d)
-include $(SRC_UTEST:%.c=%.d)
-include $(SRC_UYTEST:%.c=%.d)
