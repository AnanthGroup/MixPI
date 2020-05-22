! ------------------------------------------------ !
! --------CP2K_RUN_FORCES.F--------------------- !
! -This program is a test program that will call-- !
! -various cp2k subroutines and external modulues- !
! -to determine how to print out the atomistic --- !
! -forces at each step in the MD run and then ---- !
! -update the forces from external peturbation --- !
! -specified by external function. --------------- !
! B.Johnson: created on 12/12/19 ------------------!
! ------------------------------------------------ ! 

MODULE cp2k_run_forces
! ------------------------------------------------ !
! must call libcp2k.mod (and access libcp2k.h) --- !
! in order to call the CP2K subroutines. The ----- !
! location is specified by the library and include !
! paths in the compile statement. ---------------- !
! ------------------------------------------------ !
use libcp2k
use pimd_mod
!use force_env_types
implicit none

character(len=17)     :: version_cp2k

character(len=128)    :: cp2k_input_file, cp2k_output_file

!TYPE(force_env_type)  :: cp2k_force_env

integer                :: force_env_flag, istep
integer                :: num_classical_atoms
real(kind=8),ALLOCATABLE,DIMENSION(:)       :: atom_forces, atom_positions, atom_velocities


contains

SUBROUTINE set_up_cp2k
implicit none

! ------------------------------------------------ !
! CALL cp2k_get_version(*,*) from the CP2K library !
! to confirm correct libary linking as well as --- !
! the version of CP2K. Current code was made with- !
! VERSION 7.0 ------------------------------------ !
! ------------------------------------------------ !
CALL cp2k_get_version(version_cp2k,LEN(version_cp2k))  
write(*,*)' USING ', TRIM(version_cp2k)

! ------------------------------------------------ !
! CALL cp2k_init() from the CP2K library --------- !
! This subroutine initializes CP2K with MPI ------ !
! Must be paired with cp2k_finalize() ------------ !
! ???????????????????????????????????????????????? !
! NOT completely sure what this does, but it must- !
! Set up arrays and environments? ---------------- !
! ------------------------------------------------ !
CALL cp2k_init()

END SUBROUTINE set_up_cp2k


SUBROUTINE destroy_cp2k
! ------------------------------------------------ !
! CALL cp2 k_finalize() from the CP2K library ----- !
! This subroutine deallocates/nullifies the ------ !
! structures to prevent leakage. Pairs with ------ !
! cp2k_init() ------------------------------------ !
CALL cp2k_finalize()
! ------------------------------------------------ !

END SUBROUTINE destroy_cp2k

SUBROUTINE classical_atom_initialize
! ------------------------------------------------ !
! Subroutine that reads in the atom positions and  !
! information from the input file (.xyz or .pdb) - !
! ------------------------------------------------ !
implicit none
INTEGER          ::i, hold_i
REAL(KIND=8),ALLOCATABLE,DIMENSION(:)    :: position_array
real(kind=8)     ::pot_energy
character(len=6)  ::hold_c

! ------------------------------------------------ !
! READS IN INPUT FROM CLASSICAL.INP. THIS FILE --- !
! CONTAINS ANY INPUT RELEVANT TO THE CLASSICAL --- !
! VARIABLES. ------------------------------------- !
! ------------------------------------------------ !
write(*,*)' read in classical.inp '
open(101,file='classical.inp')
read(101,*) classical%cp2k_input_file
read(101,*) classical%cp2k_output_file
close(101)


! ------------------------------------------------ !
! CREATES THE FORCE ENVIRONEMENT IN CP2K THAT ---- !
! IS BASED UPONN THE ABOVE INPUT FILE AND THE ---- !
! FORCE FIELD SPECIFIED IN THAT FILE. ------------ !
! ------------------------------------------------ !
write(*,*)' create_force_env '
write(*,*) classical%cp2k_input_file
write(*,*) classical%cp2k_output_file
CALL cp2k_create_force_env(classical%cp2k_env_flag,trim(classical%cp2k_input_file),trim(classical%cp2k_output_file))

! ------------------------------------------------ !
! USES THE CP2K FORCE_ENV TO READ IN THE NUMBER OF !
! CLASSICAL ATOMS (NATOMS) AND THEIR RESPECTIVE -- !
! POSITIONS.
! ALSO READ IN THE CLASSICAL-CLASSICAL FORCES. --- !
! ------------------------------------------------ !
CALL cp2k_get_natom(classical%cp2k_env_flag, classical%natoms)
write(*,*)'natoms is ', classical%natoms

ALLOCATE(classical%i_position(classical%natoms,md%space_dim))
ALLOCATE(position_array(classical%natoms*md%space_dim))
CALL cp2k_get_positions(classical%cp2k_env_flag,position_array,classical%natoms*md%space_dim)

OPEN(201,file='initial_classical_position.out')
do i = 1, classical%natoms
  classical%i_position(i,1:md%space_dim) = position_array( (i-1)*md%space_dim+1:i*md%space_dim )
  write(201,'(6f12.4)') position_array( (i-1)*md%space_dim+1: i*md%space_dim), &
    position_array( (i-1)*md%space_dim+1:i*md%space_dim)/AtoBohr
enddo
CLOSE(201)

CALL cp2k_calc_energy_force(classical%cp2k_env_flag)
CALL cp2k_get_potential_energy(classical%cp2k_env_flag, pot_energy)
write(*,*)' initial potential energy is ', pot_energy

DEALLOCATE(position_array)

! ----------------------------------------------- !
! EDIT NEEDED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WRITE SUBROUTINE/DETERMINE CP2K INTERNAL ------ !
! SUBROUTINE TO OUTPUT THE PARTIAL ATOMIC CHARGES !
! FROM THE FORCE FIELDS. ------------------------ !
ALLOCATE(classical%ff_atom_charges(classical%natoms))
!classical%ff_atom_charges = 0.d0
ALLOCATE(classical%ff_atom_lj(classical%natoms,2))
!classical%ff_atom_lj(:,1) = 119.8d0 * KtoH  
!classical%ff_atom_lj(:,2) = 3.405d0 * AtoBohr 
! ----------------------------------------------- !

ALLOCATE(classical%force(classical%natoms,md%space_dim))
classical%force = 0.d0



! ----------------------------------------------- !
! EDIT NEEDED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WRITE OR FIND SUBROUTINE TO GET ATOMIC KINDS -- !
! ----------------------------------------------- !
ALLOCATE(classical%atom_kinds(classical%natoms))
ALLOCATE(classical%mass(classical%natoms))
!classical%atom_kinds = 'AR'
!classical%mass = 39.95d0 * amutome
write(*,*)' amu to me ', amutome
write(*,*)'the classical mass is ', classical%mass(1)




! ------------------------------------------------ !
! READS IN THE DATA FROM THE PDB FILE ------------ !
! ------------------------------------------------ !
open(333,file='dft_c2h4.pdb')
do i = 1, 5
  read(333,*)
enddo
do i = 1, classical%natoms
  read(333,*) hold_c, hold_i, classical%atom_kinds(i)
enddo
close(333)
do i = 1, classical%natoms
  if (trim(classical%atom_kinds(i)) == 'O') then
    classical%mass(i) = 15.998d0 * amutome
    classical%ff_atom_charges(i) = -0.75d0
    classical%ff_atom_lj(i,1) = 85.0 * KtoH
    classical%ff_atom_lj(i,2) = 3.30 * AtoBohr 
  else if (trim(classical%atom_kinds(i)) == 'Si') then
    classical%mass(i) = 28.08d0 * amutome
    classical%ff_atom_charges(i) = 1.5d0
    classical%ff_atom_lj(i,1) = 22.0 * KtoH
    classical%ff_atom_lj(i,2) = 2.30 * AtoBohr
  else
    write(*,*)' do not recognize the atom types '
  endif
enddo
 





! ----------------------------------------------- !
! INITIALIZE THE CLASSICAL VELOCITIES USING THE - !
! RANDOM NUMBER GENERATOR ON A MAXWELL-BOLTZMANN  !
! DISTRIBUTION. 
! CALCULATE THE INITIAL MOMENTA ----------------- !
! ----------------------------------------------- !
classical%beta = (1.0d0/ (md%temp * kb) )
write(*,*) 'classical beta is ', classical%beta
CALL initialize_classical_velocities

open(222,file='initial_atom_velocities.dat')
do i = 1, classical%natoms
  write(222,'(3d20.6)') classical%i_velocity(i,:)
enddo
close(222)
!open(222,file='initial_atom_velocities.dat')
!do i = 1, classical%natoms
!  read(222,*) classical%i_velocity(i,:)
!enddo
!close(222)

ALLOCATE(classical%x(classical%natoms,md%space_dim))
ALLOCATE(classical%p(classical%natoms,md%space_dim))

classical%x = classical%i_position
do i = 1, classical%natoms
  classical%p(i,:) = classical%i_velocity(i,:) * classical%mass(i)
enddo


ALLOCATE(classical%traj(0:md%nsteps/md%freq_print,classical%natoms,2*md%space_dim))
classical%traj = 0.0d0

do i = 1, classical%natoms
  classical%traj(0,i,1:md%space_dim) = classical%x(i,:)
  classical%traj(0,i,md%space_dim+1:2*md%space_dim) = classical%p(i,:)/classical%mass(i)
enddo




ALLOCATE(classical%potential(0:md%nsteps), classical%kinetic(0:md%nsteps))
classical%potential = 0.0d0
classical%kinetic = 0.0d0

END SUBROUTINE classical_atom_initialize





SUBROUTINE classical_atom_destroy

CALL cp2k_destroy_force_env(classical%cp2k_env_flag)
DEALLOCATE(classical%i_position)

END SUBROUTINE









END MODULE cp2k_run_forces
