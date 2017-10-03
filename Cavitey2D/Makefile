SHELL=/bin/sh

FC = f95

#                              for pi-VizStudio
#FC = ifort -CB -O3 -xAVX
#FC = ifort -debug

#                              for pi-computer
#FC = frtpx -Kfast -Kparallel
#FC = frtpx -Kfast -Kopenmp

LD = $(FC)

.SUFFIXES: .f90

TARGET = cavity

SRCS = main.f90 set_params.f90 alloc_variables.f90 \
       make_coefmatrix.f90 set_boundary.f90 calc_rhs.f90 \
       solve_poisson.f90 linsor.f90 uv_updates.f90 output.f90


OBJS = main.o set_params.o alloc_variables.o \
       make_coefmatrix.o set_boundary.o calc_rhs.o \
       solve_poisson.o linsor.o uv_updates.o output.o

all: $(TARGET) 

$(TARGET): $(OBJS)
	$(LD) -o $@ $(OBJS)

.f90.o: 
	$(FC) -c $<

clean:
	@rm -f $(TARGET) variables.mod *.o
