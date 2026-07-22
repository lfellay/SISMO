FC = gfortran
FFLAGS ?= -O2 -Wall -Wextra -fimplicit-none -ffree-form -ffree-line-length-none -fopenmp
MODEL ?= tests/data/release_1M5.osc.mod
REFERENCE ?= doc/release_reference_1M5.sismo
CONFIG ?= sismo.conf

SRC = \
	src/sismo_precision.f90 \
	src/sismo_types.f90 \
	src/sismo_config.f90 \
	src/sismo_io.f90 \
	src/sismo_poisson.f90 \
	src/sismo_asymptotic.f90 \
	src/sismo_core2.f90 \
	src/main.f90

OBJ = $(SRC:.f90=.o)

.PHONY: all clean test regression

all: sismo

sismo: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -f sismo src/*.o src/*.mod *.mod

test: sismo
	./tests/test_cli.sh ./sismo

regression: sismo
	@test -n "$(MODEL)" || { echo "MODEL=/path/to/model.osc.mod is required"; exit 2; }
	@test -n "$(REFERENCE)" || { echo "REFERENCE=/path/to/reference.sismo is required"; exit 2; }
	./tests/test_regression.sh ./sismo "$(MODEL)" "$(REFERENCE)" "$(CONFIG)"
