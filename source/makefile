# Definition of variables 
SHELL := /bin/bash

objects = input.o parameters.o gammamatrices.o propagators.o treelevel_gX.o \
          treelevel_g1.o treelevel_gVt.o treelevel_lY.o treelevel_l1.o \
          treelevel_lVt.o treelevel_4fermions.o vertices.o oneloop_gX.o \
          oneloop_g1.o oneloop_gVt.o oneloop_lY.o oneloop_l1.o oneloop_lVt.o \
          oneloop_4f_G1.o oneloop_4f_G2.o oneloop_4f_Gk1.o \
		  oneloop_4f_Gk2.o compute.o main.o
obs_mods= treelevel_4fermions.mod treelevel_g1.mod treelevel_gVt.mod \
          treelevel_gX.mod treelevel_gX_p.mod treelevel_l1.mod treelevel_lVt.mod \
          treelevel_lY.mod oneloop_4f_G1.mod oneloop_4f_G2.mod oneloop_4f_Gk1.mod\
		  oneloop_4f_Gk2.mod oneloop_g1.mod oneloop_gVt.mod \
          oneloop_gX.mod oneloop_gX_p.mod oneloop_l1.mod oneloop_lVt.mod oneloop_lY.mod 

obs_source= treelevel_4fermions.f90 treelevel_g1.f90 treelevel_gVt.f90 \
          treelevel_gX.f90 treelevel_gX_p.f90 treelevel_l1.f90 treelevel_lVt.f90 \
          treelevel_lY.f90 oneloop_4f_G1.f90 oneloop_4f_G2.f90 oneloop_4f_Gk1.f90\
		  oneloop_4f_Gk2.f90 oneloop_g1.f90 oneloop_gVt.f90 \
          oneloop_gX.f90 oneloop_gX_p.f90 oneloop_l1.f90 oneloop_lVt.f90 oneloop_lY.f90 

f90comp = ifort
switch = -O2 -openmp
#libs = -L$(LAPACK_LIB) -llapack -L$(BLAS_LIB) -lblas -lgfortran
libs = -llapack -lblas



pt_simple.out: $(objects)
	$(f90comp) -o pt_simple.out $(switch) $(objects) $(libs)

compute.mod: compute.o compute.f90
	$(f90comp) -c $(switch) compute.f90
compute.o: $(obs_mods) $(obs_source)
	$(f90comp) -c $(switch) compute.f90

main.mod: main.o main.f90
	$(f90comp) -c $(switch) main.f90
main.o:  gammamatrices.mod input.mod parameters.mod propagators.mod compute.mod main.f90
	$(f90comp) -c $(switch) main.f90


input.mod: input.o input.f90
	$(f90comp) -c $(switch) input.f90
input.o: input.f90
	$(f90comp) -c $(switch) input.f90

parameters.mod: parameters.o parameters.f90
	$(f90comp) -c $(switch) parameters.f90
parameters.o: parameters.f90
	$(f90comp) -c $(switch) parameters.f90

gammamatrices.mod: gammamatrices.o gammamatrices.f90
	$(f90comp) -c $(switch) gammamatrices.f90
gammamatrices.o: parameters.mod gammamatrices.f90
	$(f90comp) -c $(switch) gammamatrices.f90

propagators.mod: propagators.o propagators.f90
	$(f90comp) -c $(switch) propagators.f90
propagators.o: gammamatrices.mod input.mod parameters.mod propagators.f90
	$(f90comp) -c $(switch) propagators.f90

treelevel_4fermions.mod: treelevel_4fermions.o treelevel_4fermions.f90
	$(f90comp) -c $(switch) treelevel_4fermions.f90
treelevel_4fermions.o: gammamatrices.mod input.mod parameters.mod propagators.mod treelevel_4fermions.f90
	$(f90comp) -c $(switch) treelevel_4fermions.f90

treelevel_g1.mod: treelevel_g1.o treelevel_g1.f90
	$(f90comp) -c $(switch) treelevel_g1.f90
treelevel_g1.o: gammamatrices.mod input.mod parameters.mod propagators.mod treelevel_g1.f90
	$(f90comp) -c $(switch) treelevel_g1.f90

treelevel_gVt.mod: treelevel_gVt.o treelevel_gVt.f90
	$(f90comp) -c $(switch) treelevel_gVt.f90
treelevel_gVt.o: gammamatrices.mod input.mod parameters.mod propagators.mod treelevel_gVt.f90
	$(f90comp) -c $(switch) treelevel_gVt.f90

treelevel_gX.mod: treelevel_gX.o treelevel_gX.f90
	$(f90comp) -c $(switch) treelevel_gX.f90
treelevel_gX.o: gammamatrices.mod input.mod parameters.mod propagators.mod treelevel_gX.f90
	$(f90comp) -c $(switch) treelevel_gX.f90

treelevel_gX_p.mod: treelevel_gX_p.o treelevel_gX_p.f90
	$(f90comp) -c $(switch) treelevel_gX_p.f90
treelevel_gX_p.o: gammamatrices.mod input.mod parameters.mod propagators.mod treelevel_gX_p.f90
	$(f90comp) -c $(switch) treelevel_gX_p.f90

treelevel_l1.mod: treelevel_l1.o treelevel_l1.f90
	$(f90comp) -c $(switch) treelevel_l1.f90
treelevel_l1.o: gammamatrices.mod input.mod parameters.mod propagators.mod treelevel_l1.f90
	$(f90comp) -c $(switch) treelevel_l1.f90

treelevel_lVt.mod: treelevel_lVt.o treelevel_lVt.f90
	$(f90comp) -c $(switch) treelevel_lVt.f90
treelevel_lVt.o: gammamatrices.mod input.mod parameters.mod propagators.mod treelevel_lVt.f90
	$(f90comp) -c $(switch) treelevel_lVt.f90

treelevel_lY.mod: treelevel_lY.o treelevel_lY.f90
	$(f90comp) -c $(switch) treelevel_lY.f90
treelevel_lY.o: gammamatrices.mod input.mod parameters.mod propagators.mod treelevel_lY.f90
	$(f90comp) -c $(switch) treelevel_lY.f90


vertices.mod: vertices.o vertices.f90
	$(f90comp) -c $(switch) vertices.f90
vertices.o: gammamatrices.mod input.mod parameters.mod vertices.f90
	$(f90comp) -c $(switch) vertices.f90

oneloop_4f_G1.mod: oneloop_4f_G1.o oneloop_4f_G1.f90
	$(f90comp) -c $(switch) oneloop_4f_G1.f90
oneloop_4f_G1.o: gammamatrices.mod input.mod parameters.mod propagators.mod vertices.mod treelevel_4fermions.f90
	$(f90comp) -c $(switch) oneloop_4f_G1.f90

oneloop_4f_G2.mod: oneloop_4f_G2.o oneloop_4f_G2.f90
	$(f90comp) -c $(switch) oneloop_4f_G2.f90
oneloop_4f_G2.o: gammamatrices.mod input.mod parameters.mod propagators.mod vertices.mod treelevel_4fermions.f90
	$(f90comp) -c $(switch) oneloop_4f_G2.f90

oneloop_4f_Gk1.mod: oneloop_4f_Gk1.o oneloop_4f_Gk1.f90
	$(f90comp) -c $(switch) oneloop_4f_Gk1.f90
oneloop_4f_Gk1.o: gammamatrices.mod input.mod parameters.mod propagators.mod vertices.mod treelevel_4fermions.f90
	$(f90comp) -c $(switch) oneloop_4f_Gk1.f90

oneloop_4f_Gk2.mod: oneloop_4f_Gk2.o oneloop_4f_Gk2.f90
	$(f90comp) -c $(switch) oneloop_4f_Gk2.f90
oneloop_4f_Gk2.o: gammamatrices.mod input.mod parameters.mod propagators.mod vertices.mod treelevel_4fermions.f90
	$(f90comp) -c $(switch) oneloop_4f_Gk2.f90


oneloop_g1.mod: oneloop_g1.o oneloop_g1.f90
	$(f90comp) -c $(switch) oneloop_g1.f90
oneloop_g1.o: gammamatrices.mod input.mod parameters.mod propagators.mod vertices.mod treelevel_g1.f90
	$(f90comp) -c $(switch) oneloop_g1.f90

oneloop_gVt.mod: oneloop_gVt.o oneloop_gVt.f90
	$(f90comp) -c $(switch) oneloop_gVt.f90
oneloop_gVt.o: gammamatrices.mod input.mod parameters.mod propagators.mod vertices.mod treelevel_gVt.f90
	$(f90comp) -c $(switch) oneloop_gVt.f90

oneloop_gX.mod: oneloop_gX.o oneloop_gX.f90
	$(f90comp) -c $(switch) oneloop_gX.f90
oneloop_gX.o: gammamatrices.mod input.mod parameters.mod propagators.mod vertices.mod treelevel_gX.f90
	$(f90comp) -c $(switch) oneloop_gX.f90

oneloop_gX_p.mod: oneloop_gX_p.o oneloop_gX_p.f90
	$(f90comp) -c $(switch) oneloop_gX_p.f90
oneloop_gX_p.o: gammamatrices.mod input.mod parameters.mod propagators.mod vertices.mod treelevel_gX_p.f90
	$(f90comp) -c $(switch) oneloop_gX_p.f90

oneloop_l1.mod: oneloop_l1.o oneloop_l1.f90
	$(f90comp) -c $(switch) oneloop_l1.f90
oneloop_l1.o: gammamatrices.mod input.mod parameters.mod propagators.mod vertices.mod treelevel_l1.f90
	$(f90comp) -c $(switch) oneloop_l1.f90

oneloop_lVt.mod: oneloop_lVt.o oneloop_lVt.f90
	$(f90comp) -c $(switch) oneloop_lVt.f90
oneloop_lVt.o: gammamatrices.mod input.mod parameters.mod propagators.mod vertices.mod treelevel_lVt.f90
	$(f90comp) -c $(switch) oneloop_lVt.f90

oneloop_lY.mod: oneloop_lY.o oneloop_lY.f90
	$(f90comp) -c $(switch) oneloop_lY.f90
oneloop_lY.o: gammamatrices.mod input.mod parameters.mod propagators.mod vertices.mod treelevel_lY.f90
	$(f90comp) -c $(switch) oneloop_lY.f90







clean:
	rm *.mod 
	rm *.o
