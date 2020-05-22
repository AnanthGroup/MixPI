FC = mpif90

#libraries need for compilation. First line are CP2K general flags. Second line
#are flags that match specific compilation version of CP2K.
LIB_CP2K = -lcp2k -ldbcsr 
LIB = $(LIB_CP2K) -lfftw3xf_gnu -lmkl_scalapack_lp64 -lmkl_gf_lp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_lp64
FLAGS = -O3 -ffree-form -C

#Global paths for the cp2k and dbcsr libraries. This is the location CP2K 
#puts the library directories
CP2K_LIB = /people/john708/cp2k_editlib/lib/constance-gf-mkl-mvapich/sopt/
DBCSR_LIB = /people/john708/cp2k_editlib/lib/constance-gf-mkl-mvapich/sopt/exts/dbcsr/

#Global paths for the cp2k and dbcsr include/object directories. They contain
#the .F, .mod, and .o files for ALL of the subroutines/modules used by CP2K.
#Must move the libcp2k.h file from the src/start folder into the CP2K_INC folder
#or have another path that points towards this directory. 
CP2K_INC = /people/john708/cp2k_editlib/obj/constance-gf-mkl-mvapich/sopt/
DBCSR_INC = /people/john708/cp2k_editlib/obj/constance-gf-mkl-mvapich/sopt/exts/dbcsr

#Path from the compile statement of CP2K. It is the base part of the path for the FFTW,
#SCALAPACK, MKL, and BLACS libraries.
MKL_PATH = /share/apps/intel/2015u2/mkl


# --------------------------
#
#IF have additional .MOD files you want to compile, include additional
#definitions for MO2, MO3, etc.
MAIN = cp2k_pimd
MO3 = variable_types
MO2 = cp2k_run_forces
MO1 = pimd_mod
.DEFAULT_GOAL := $(MAIN).out


#Makes .o files from all of the .F files. Include additional 
#routines for additional files (it's own object generation as
#well as inclusion OBJ defintion).
$(MO3).o: $(MO3).F
	$(FC) -c $(MO3).F $(FLAGS)
$(MO1).o: $(MO1).F
	$(FC) -c $(MO1).F $(MO3).F -I$(CP2K_INC) -I$(DBCSR_INC) -L$(CP2K_LIB) -L$(DBCSR_LIB) -L$(MKL_PATH)/interfaces/fftw3xf/ -L$(MKL_PATH)/lib/intel64/ $(LIB) $(FLAGS)
$(MO2).o: $(MO2).F 
	$(FC) -c $(MO2).F -I$(CP2K_INC) -I$(DBCSR_INC) -L$(CP2K_LIB) -L$(DBCSR_LIB) -L$(MKL_PATH)/interfaces/fftw3xf/ -L$(MKL_PATH)/lib/intel64/ $(LIB) $(FLAGS)
$(MAIN).o: $(MAIN).F 
	$(FC) -c $(MAIN).F -I$(CP2K_INC) -I$(DBCSR_INC) -L$(CP2K_LIB) -L$(DBCSR_LIB) -L$(MKL_PATH)/interfaces/fftw3xf/ -L$(MKL_PATH)/lib/intel64/ $(MO1).F $(MO2).F $(MO3).F $(LIB) $(FLAGS)
 
OBJ = $(MAIN).o $(MO1).o $(MO2).o $(MO3).o

#Compile statement. Makes an executable (MAIN).out that is dependent on (MAIN).o and $(MO1).o (must list mod files so that
#they are generated before they are used. Add any additional .mod files to the list.
#The 3/4 library statements point to the locations of liblfftw3xf.a and  libmkl_scalapack_lp64, libmkl_gf_lp64, libmkl_core, libmkl_sequential, 
#libmkl_blacs_intelmpi_lp64. (See the lib flags.)
#For dependencies, library flags come AFTER main.o
$(MAIN).out	: $(MO3).o $(MO2).o $(MO1).o $(MAIN).o
	$(FC) -I$(CP2K_INC) -I$(DBCSR_INC) -L$(CP2K_LIB) -L$(DBCSR_LIB) -L$(MKL_PATH)/interfaces/fftw3xf/ -L$(MKL_PATH)/lib/intel64/ $(OBJ) $(LIB) $(FLAGS) -o $@  


#Removes the .mod, .o, and .out files
clean:
	rm *.o
	rm *.mod
	rm $(MAIN).out
