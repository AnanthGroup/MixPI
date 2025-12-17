FC = mpif90

#libraries need for compilation. First line are CP2K general flags. Second line
#are flags that match specific compilation version of CP2K.
LIB_CP2K = -lcp2k -ldbcsr
LIB = $(LIB_CP2K) -lplumed -lxcf03 -lxc -lGXCommon -lgx_ac -lgx_minimax -lint2 -lxsmmf -lxsmmf.1 -lxsmmext -lxsmm -lhdf5_fortran -lhdf5_f90cstub -lhdf5 -lz -lsz  -lcosma_prefixed_pxgemm -lcosma -lcosta -lfftw3_mpi -lfftw3_omp -lfftw3 -lscalapack -lgsl -lopenblas
FLAGS = -O3 -ffree-form -fopenmp -fbacktrace -fallow-argument-mismatch

HDF5_HOME = $(shell brew --prefix hdf5)
SZ_HOME = $(shell brew --prefix libaec)
#-lhdf5_fortran -lhdf5_f90cstub -lhdf5 -lsz -ltblite -ldftd4 -ls-dftd3 -lmulticharge -lmctc-lib -ltoml-f -ltrexio 

#Global paths for the cp2k and dbcsr libraries. This is the location CP2K 
#puts the library directories
CP2K_LIB = /Users/john708/Documents/Software/cp2k/lib/Darwin-gnu-arm64/psmp 
DBCSR_LIB = /Users/john708/Documents/Software/cp2k/lib/Darwin-gnu-arm64/psmp/exts/dbcsr
BREW_LIB += -L$(HDF5_HOME)/lib
BREW_LIB += -L$(SZ_HOME)/lib

#Global paths for the cp2k and dbcsr include/object directories. They contain
#the .F, .mod, and .o files for ALL of the subroutines/modules used by CP2K.
#Must move the libcp2k.h file from the src/start folder into the CP2K_INC folder
#or have another path that points towards this directory. 
CP2K_INC = /Users/john708/Documents/Software/cp2k/obj/Darwin-gnu-arm64/psmp
DBCSR_INC = /Users/john708/Documents/Software/cp2k/lib/Darwin-gnu-arm64/psmp/exts/dbcsr
BREW_INC += -I$(HDF5_HOME)/include
BREW_INC += -I$(SZ_HOME)/include

#Path from the compile statement of CP2K. It is the base part of the path for the FFTW,
#SCALAPACK, MKL, and BLACS libraries.
#MKL_PATH = /share/apps/intel/2020u4/mkl



# --------------------------
#
#IF have additional .MOD files you want to compile, include additional
#definitions for MO2, MO3, etc.
MAIN = cp2k_pimd
MO4 = variable_types
MO3 = input_destroy
MO2 = input_setup
MO1 = md_run
.DEFAULT_GOAL := $(MAIN).out


#Makes .o files from all of the .F files. Include additional 
#routines for additional files (it's own object generation as
#well as inclusion OBJ defintion).
$(MO4).o: $(MO4).f90
	$(FC) -c $(MO4).f90 $(FLAGS)
$(MO3).o: $(MO3).f90 $(MO4).f90
	$(FC) -c $(MO3).f90 -I$(CP2K_INC) -I$(DBCSR_INC) $(BREW_INC)  $(LIB) $(FLAGS)
$(MO2).o: $(MO2).f90 $(MO3).f90 $(MO4).f90
	$(FC) -c $(MO2).f90 -I$(CP2K_INC) -I$(DBCSR_INC) $(BREW_INC)  $(LIB) $(FLAGS)
$(MO1).o: $(MO1).f90 $(MO4).f90
	$(FC) -c $(MO1).f90 -I$(CP2K_INC) -I$(DBCSR_INC) $(BREW_INC)  $(LIB) $(FLAGS)
$(MAIN).o: $(MAIN).f90 $(MO1).f90 $(MO2).f90 $(MO3).f90 $(MO4).f90 
	$(FC) -c $(MAIN).f90 -I$(CP2K_INC) -I$(DBCSR_INC) $(BREW_INC) $(LIB) $(FLAGS)
 
OBJ = $(MAIN).o $(MO1).o $(MO2).o $(MO3).o $(MO4).o

#Compile statement. Makes an executable (MAIN).out that is dependent on (MAIN).o and $(MO1).o (must list mod files so that
#they are generated before they are used. Add any additional .mod files to the list.
#The 3/4 library statements point to the locations of liblfftw3xf.a and  libmkl_scalapack_lp64, libmkl_gf_lp64, libmkl_core, libmkl_sequential, 
#libmkl_blacs_intelmpi_lp64. (See the lib flags.)
#For dependencies, library flags come AFTER main.o
$(MAIN).out	: $(MO4).o $(MO3).o $(MO2).o $(MO1).o $(MAIN).o
	$(FC) -I$(CP2K_INC) -I$(DBCSR_INC) $(BREW_INC) -L$(CP2K_LIB) -L$(DBCSR_LIB) $(BREW_LIB) $(OBJ) $(LIB) $(FLAGS) -lstdc++ -lz -ldl  -o $@ 


#Removes the .mod, .o, and .out files
clean:
	rm *.o
	rm *.mod
	rm $(MAIN).out
