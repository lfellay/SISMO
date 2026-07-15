FC = gfortran
FFLAGS ?= -O2 -Wall -Wextra -fimplicit-none -ffree-form -ffree-line-length-none -fopenmp

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

.PHONY: all clean

all: sismo

sismo: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -f sismo src/*.o src/*.mod *.mod
