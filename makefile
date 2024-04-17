FC = mpif90

#libraries need for compilation. First line are CP2K general flags. Second line
#are flags that match specific compilation version of CP2K.
LIB_CP2K = -lcp2k -ldbcsr
LIB = $(LIB_CP2K) -lgfortran 
FLAGS = -O3 -ffree-form -fopenmp

#Global paths for the cp2k and dbcsr libraries. This is the location CP2K 
#puts the library directories
CP2K_LIB = /people/john708/bin/cp2k-v2023/cp2k-2023.2-tip4p/lib/deception-gf-mkl-openmpi-libint/psmp
DBCSR_LIB = /people/john708/bin/cp2k-v2023/cp2k-2023.2-tip4p/lib/deception-gf-mkl-openmpi-libint/psmp/exts/dbcsr/

#Global paths for the cp2k and dbcsr include/object directories. They contain
#the .F, .mod, and .o files for ALL of the subroutines/modules used by CP2K.
#Must move the libcp2k.h file from the src/start folder into the CP2K_INC folder
#or have another path that points towards this directory. 
CP2K_INC = /people/john708/bin/cp2k-v2023/cp2k-2023.2-tip4p/obj/deception-gf-mkl-openmpi-libint/psmp
DBCSR_INC = /people/john708/bin/cp2k-v2023/cp2k-2023.2-tip4p/obj/deception-gf-mkl-openmpi-libint/psmp/exts/dbcsr

#Path from the compile statement of CP2K. It is the base part of the path for the FFTW,
#SCALAPACK, MKL, and BLACS libraries.
MKL_PATH = /share/apps/intel/2020u4/mkl


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
$(MO4).o: $(MO4).F
	$(FC) -c $(MO4).F $(FLAGS)
$(MO3).o: $(MO3).F $(MO4).F
	$(FC) -c $(MO3).F -I$(CP2K_INC) -I$(DBCSR_INC)  -I$(MKL_PATH)/include -I$(MKL_PATH)/include/fftw $(LIB) $(FLAGS)
$(MO2).o: $(MO2).F $(MO3).F $(MO4).F
	$(FC) -c $(MO2).F -I$(CP2K_INC) -I$(DBCSR_INC)  -I$(MKL_PATH)/include -I$(MKL_PATH)/include/fftw $(LIB) $(FLAGS)
$(MO1).o: $(MO1).F $(MO4).F
	$(FC) -c $(MO1).F -I$(CP2K_INC) -I$(DBCSR_INC)  -I$(MKL_PATH)/include -I$(MKL_PATH)/include/fftw $(LIB) $(FLAGS)
$(MAIN).o: $(MAIN).F $(MO1).F $(MO2).F $(MO3).F $(MO4).F 
	$(FC) -c $(MAIN).F -I$(CP2K_INC) -I$(DBCSR_INC) $(LIB) $(FLAGS)
 
OBJ = $(MAIN).o $(MO1).o $(MO2).o $(MO3).o $(MO4).o

#Compile statement. Makes an executable (MAIN).out that is dependent on (MAIN).o and $(MO1).o (must list mod files so that
#they are generated before they are used. Add any additional .mod files to the list.
#The 3/4 library statements point to the locations of liblfftw3xf.a and  libmkl_scalapack_lp64, libmkl_gf_lp64, libmkl_core, libmkl_sequential, 
#libmkl_blacs_intelmpi_lp64. (See the lib flags.)
#For dependencies, library flags come AFTER main.o
$(MAIN).out	: $(MO4).o $(MO3).o $(MO2).o $(MO1).o $(MAIN).o
	$(FC) -I$(CP2K_INC) -I$(DBCSR_INC) -I$(MKL_PATH)/include -I$(MKL_PATH)/include/fftw -L$(CP2K_LIB) -L$(DBCSR_LIB) $(OBJ) $(LIB) -Wl,--start-group $(MKL_PATH)/lib/intel64/libmkl_scalapack_lp64.a $(MKL_PATH)/lib/intel64/libmkl_blacs_openmpi_lp64.a $(MKL_PATH)/lib/intel64/libmkl_gf_lp64.a $(MKL_PATH)/lib/intel64/libmkl_gnu_thread.a $(MKL_PATH)/lib/intel64/libmkl_core.a -Wl,--end-group $(FLAGS) -lstdc++ -lz -ldl -lgomp -lpthread -o $@ 


#Removes the .mod, .o, and .out files
clean:
	rm *.o
	rm *.mod
	rm $(MAIN).out
