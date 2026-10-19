#-------------------------------------------------------------------------------
# SPEX/SPEX_SDIPPS/Makefile: build and test SDIPPS
#-------------------------------------------------------------------------------

# The SDIPPS library code (Source/*.c) is compiled into libspex by SPEX's own
# build, unchanged: SPEX's CMakeLists.txt compiles every SPEX*/Source/*.c.
# This Makefile builds that library, then the test program, and runs it.
#
#   make            build libspex (with SDIPPS) and the test, run the tests
#   make library    build libspex only (same as "make local" in SPEX)
#   make test       build the test program only
#   make run        run the tests (every option combination)
#   make lp         stress test on a simplex run: SDIPPS on the last N
#                   updates of the pivot log LOG (metrics: <name>_sdipps.csv)
#   make lp_spex    the same with SPEX's update, SPEX_update_lu_colrep
#                   (metrics: <name>_spex.csv)
#   make lp_refactor  the baseline: SPEX LU from scratch every Y updates
#                   (metrics: <name>_refactor.csv)
#   make clean      remove the test programs
#
# The LP targets take LOG, N and Y, e.g.
#   make lp LOG=Test/Data/cycle_bases.txt N=100
# Defaults: LOG=Test/Data/co9_bases.txt, N=50, Y=5.  <name> is LOG's file
# name without "_bases.txt" (co9, cycle).

SPEX = ..
ROOT = ../..

CC     ?= gcc
CFLAGS ?= -O2 -Wall -Wextra -Wno-unused-parameter
INC     = -IInclude -I$(SPEX)/Include -I$(SPEX)/Demo \
          -I$(ROOT)/SuiteSparse_config
LIBS    = $(SPEX)/build/libspex.a $(ROOT)/COLAMD/build/libcolamd.a \
          $(ROOT)/AMD/build/libamd.a \
          $(ROOT)/SuiteSparse_config/build/libsuitesparseconfig.a \
          -lmpfr -lgmp -lm -fopenmp

# SPEX demo utilities reused by the test: reading a triplet matrix, checking
# a solution exactly, and printing error messages
DEMO = $(SPEX)/Demo/Utilities/spex_demo_tripread.c \
       $(SPEX)/Demo/Utilities/spex_demo_check_solution.c \
       $(SPEX)/Demo/Utilities/spex_demo_determine_error.c

TEST = Test/sdipps_test

all: library test run

library:
	( cd $(SPEX) && $(MAKE) local )

test: $(TEST)

$(TEST): Test/sdipps_test.c Include/SDIPPS.h $(SPEX)/build/libspex.a
	$(CC) $(CFLAGS) $(INC) Test/sdipps_test.c $(DEMO) $(LIBS) -o $@

# each data set under every combination of the options
OPTIONS = "" "-g" "-l" "-e" "-g -e" "-l -e" \
          "-c" "-c -g" "-c -l" "-c -e" "-c -g -e" "-c -l -e"

run: $(TEST)
	@for opts in $(OPTIONS) ; do \
	    ./$(TEST) $(SPEX)/ExampleMats/example.mat.txt Test/Data/example.upd.txt $$opts || exit 1 ; \
	    ./$(TEST) Test/Data/pm1_12.mat.txt Test/Data/pm1_12.upd.txt $$opts || exit 1 ; \
	    ./$(TEST) $(SPEX)/ExampleMats/10teams.mat.txt Test/Data/10teams.upd.txt $$opts || exit 1 ; \
	done
	@echo "SDIPPS: all tests passed"

# stress tests on the bases of a simplex run (a QSopt_ex pivot log)
LP_TEST     = Test/sdipps_lp_test
LP_REFACTOR = Test/sdipps_lp_refactor
LP_COMMON   = Test/sdipps_pivot_log.c Test/sdipps_pivot_log.h

LOG  ?= Test/Data/co9_bases.txt
N    ?= 50
Y    ?= 5
NAME  = $(subst _bases,,$(basename $(notdir $(LOG))))

$(LP_TEST): Test/sdipps_lp_test.c $(LP_COMMON) Include/SDIPPS.h \
        $(SPEX)/build/libspex.a
	$(CC) $(CFLAGS) $(INC) Test/sdipps_lp_test.c Test/sdipps_pivot_log.c \
	    $(DEMO) $(LIBS) -o $@

$(LP_REFACTOR): Test/sdipps_lp_refactor.c $(LP_COMMON) $(SPEX)/build/libspex.a
	$(CC) $(CFLAGS) $(INC) Test/sdipps_lp_refactor.c Test/sdipps_pivot_log.c \
	    $(DEMO) $(LIBS) -o $@

lp: $(LP_TEST)
	./$(LP_TEST) $(LOG) -n $(N) -o $(NAME)_sdipps.csv

lp_spex: $(LP_TEST)
	./$(LP_TEST) $(LOG) -n $(N) -s -o $(NAME)_spex.csv

lp_refactor: $(LP_REFACTOR)
	./$(LP_REFACTOR) $(LOG) -n $(N) -y $(Y) -o $(NAME)_refactor.csv

clean:
	rm -f $(TEST) $(LP_TEST) $(LP_REFACTOR)

.PHONY: all library test run lp lp_spex lp_refactor clean
