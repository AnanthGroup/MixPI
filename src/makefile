FC = mpif90

#libraries need for compilation. First line are CP2K general flags. Second line
#are flags that match specific compilation version of CP2K.
LIB_CP2K = -lcp2k -ldbcsr
#libraries from cp2k arch file, example below
LIB = $(LIB_CP2K) -lgfortran -lspla -lspfft -lsymspg -lvdwxc -lhdf5 -lhdf5_hl -lz -lgsl -lelpa_openmp -lcosma_prefixed_pxgemm -lcosma -lcosta  -lxsmmf -lxsmm -ldl -lpthread -lxcf03 -lxc -lint2 -lfftw3_mpi -lfftw3 -lfftw3_omp   -lmpi  -L/share/apps/intel/2020u4/mkl/lib/intel64 -Wl,-rpath,/share/apps/intel/2020u4/mkl/lib/intel64 -lmkl_scalapack_lp64 -Wl,--start-group -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_openmpi_lp64 -Wl,--end-group -lpthread -lm -ldl -lvori -lstdc++ -lstdc++
FLAGS = -O3 -ffree-form -fopenmp

#Global paths for the cp2k and dbcsr libraries. This is the location CP2K 
#puts the library directories (example below)
CP2K_LIB = /bin/cp2k-2023.2/lib/local/psmp 
DBCSR_LIB = /bin/cp2k-2023.2/lib/local/psmp/exts/dbcsr

#Global paths for the cp2k and dbcsr include/object directories. They contain
#the .F, .mod, and .o files for ALL of the subroutines/modules used by CP2K.
#Must move the libcp2k.h file from the src/start folder into the CP2K_INC folder
#or have another path that points towards this directory. Example below 
CP2K_INC = /bin/cp2k-2023.2/obj/local/psmp 
DBCSR_INC = /bin/cp2k-2023.2/obj/local/psmp/exts/dbcsr



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
	$(FC) -c $(MO3).F -I$(CP2K_INC) -I$(DBCSR_INC)  $(LIB) $(FLAGS)
$(MO2).o: $(MO2).F $(MO3).F $(MO4).F
	$(FC) -c $(MO2).F -I$(CP2K_INC) -I$(DBCSR_INC)  $(LIB) $(FLAGS)
$(MO1).o: $(MO1).F $(MO4).F
	$(FC) -c $(MO1).F -I$(CP2K_INC) -I$(DBCSR_INC)   $(LIB) $(FLAGS)
$(MAIN).o: $(MAIN).F $(MO1).F $(MO2).F $(MO3).F $(MO4).F 
	$(FC) -c $(MAIN).F -I$(CP2K_INC) -I$(DBCSR_INC) $(LIB) $(FLAGS)
 
OBJ = $(MAIN).o $(MO1).o $(MO2).o $(MO3).o $(MO4).o

#Compile statement. Makes an executable (MAIN).out that is dependent on (MAIN).o and $(MO1).o (must list mod files so that
#they are generated before they are used. Add any additional .mod files to the list.
#For dependencies, library flags come AFTER main.o
$(MAIN).out	: $(MO4).o $(MO3).o $(MO2).o $(MO1).o $(MAIN).o
	$(FC) -I$(CP2K_INC) -I$(DBCSR_INC) -L$(CP2K_LIB) -L$(DBCSR_LIB) $(OBJ) $(LIB) $(FLAGS) -lstdc++ -lz -ldl -lgomp -lpthread -o $@ 


#Removes the .mod, .o, and .out files
clean:
	rm *.o
	rm *.mod
	rm $(MAIN).out
