CC      ?= gcc
CFLAGS  ?= -O2 -std=c99 -Wall -Wextra -pedantic
LDLIBS  = -lm

# make OMP=1 : run the windows of east_reus in parallel with OpenMP
ifeq ($(OMP),1)
CFLAGS += -fopenmp
LDLIBS += -fopenmp
endif

# exact diagonalization (continuous time)
ED_SRC  = src/model.c src/basis.c src/tgen.c src/eig_dense.c src/lanczos.c \
          src/cg.c src/ed.c src/cli_util.c
# discrete-time simulation and analysis
SIM_SRC = src/rng.c src/input.c src/dyn.c src/stats.c src/tps.c src/clone.c \
          src/tm.c src/mbar.c src/peaks.c src/reus.c
HEADERS = $(wildcard src/*.h)

PROGS = bin/east_ed bin/east_test bin/east_tps bin/east_clone bin/east_tm \
        bin/east_mbar bin/east_reus bin/east_validate bin/east_edan

all: $(PROGS)

bin:
	mkdir -p bin

bin/east_ed: $(ED_SRC) src/east_ed.c $(HEADERS) | bin
	$(CC) $(CFLAGS) -o $@ $(ED_SRC) src/east_ed.c $(LDLIBS)

bin/east_test: $(ED_SRC) $(SIM_SRC) src/east_test.c $(HEADERS) | bin
	$(CC) $(CFLAGS) -o $@ $(ED_SRC) $(SIM_SRC) src/east_test.c $(LDLIBS)

bin/east_tps: $(ED_SRC) $(SIM_SRC) src/east_tps.c $(HEADERS) | bin
	$(CC) $(CFLAGS) -o $@ $(ED_SRC) $(SIM_SRC) src/east_tps.c $(LDLIBS)

bin/east_clone: $(ED_SRC) $(SIM_SRC) src/east_clone.c $(HEADERS) | bin
	$(CC) $(CFLAGS) -o $@ $(ED_SRC) $(SIM_SRC) src/east_clone.c $(LDLIBS)

bin/east_tm: $(ED_SRC) $(SIM_SRC) src/east_tm.c $(HEADERS) | bin
	$(CC) $(CFLAGS) -o $@ $(ED_SRC) $(SIM_SRC) src/east_tm.c $(LDLIBS)

bin/east_mbar: $(ED_SRC) $(SIM_SRC) src/east_mbar.c $(HEADERS) | bin
	$(CC) $(CFLAGS) -o $@ $(ED_SRC) $(SIM_SRC) src/east_mbar.c $(LDLIBS)

bin/east_reus: $(ED_SRC) $(SIM_SRC) src/east_reus.c $(HEADERS) | bin
	$(CC) $(CFLAGS) -o $@ $(ED_SRC) $(SIM_SRC) src/east_reus.c $(LDLIBS)

bin/east_validate: $(ED_SRC) $(SIM_SRC) src/east_validate.c $(HEADERS) | bin
	$(CC) $(CFLAGS) -o $@ $(ED_SRC) $(SIM_SRC) src/east_validate.c $(LDLIBS)

bin/east_edan: src/east_edan.c | bin
	$(CC) $(CFLAGS) -o $@ src/east_edan.c $(LDLIBS)

test: bin/east_test
	./bin/east_test

validate: bin/east_validate
	./bin/east_validate

clean:
	rm -rf bin

.PHONY: all test validate clean
