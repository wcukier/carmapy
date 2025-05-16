# Sub makefile for test files

# This is intended to be included by another makefile to actually
# build the system. It has all the dependency information for the
# files in the test tree.

# NOTE: In the future (or in parnet models) this could be replaced
# by automatic dependency generation.

atmosphere_mod.o : atmosphere_mod.F90
	$(FORTRAN) $(FFLAGS) -c $<

carma_testutils.o : carma_testutils.F90
	$(FORTRAN) $(FFLAGS) -c $<

carma_carmapy.o : carma_carmapy.F90 atmosphere_mod.mod carma_mod.mod carma_planet_mod.mod carma_condensate_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<