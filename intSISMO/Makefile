FC = gfortran

debug = yes

ifeq ($(debug),yes)
  FFLAGS = -ffree-form -ffpe-trap=invalid,zero,overflow -fimplicit-none \
    -Wall -g -fbounds-check -fbacktrace -cpp -fdollar-ok -fcheck=all \
    -fno-align-commons -frecursive
else
  FFLAGS = -ffree-form -fimplicit-none -O3 -cpp -fdollar-ok -fno-align-commons
endif

TARGET = intSISMO
OBJECTS = physctes.o madmodlib.o intSISMO_lib.o intSISMO_main.o

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $(OBJECTS)

physctes.o: physctes.f95
	$(FC) $(FFLAGS) -c $< -o $@

madmodlib.o: madmodlib.f95
	$(FC) $(FFLAGS) -c $< -o $@

intSISMO_lib.o: intSISMO_lib.f95 physctes.o madmodlib.o
	$(FC) $(FFLAGS) -c $< -o $@

intSISMO_main.o: intSISMO_main.f95 intSISMO_lib.o madmodlib.o
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -f $(TARGET) *.o *.mod
