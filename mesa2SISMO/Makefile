FC = gfortran
FFLAGS ?= -O2 -Wall -Wextra -fimplicit-none -ffree-form

TARGET = mesa2SISMO
SOURCE = mesa2SISMO.f90

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SOURCE)
	$(FC) $(FFLAGS) -o $@ $<

clean:
	rm -f $(TARGET) *.o *.mod
