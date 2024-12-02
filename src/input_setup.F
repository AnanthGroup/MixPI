MODULE input_setup
use variable_types
use libcp2k
implicit none


!TYPE(path_integral_env)                ::path_integral
!TYPE(classical_type)                   ::classical
!TYPE(md_info)                           ::md
!TYPE(force_field)                       ::ff


contains


! -------------------------------------------- !
! Reads in the variables from the md_run.inp - !
! input file. Generally concerned with the --- !
! details of the MD run.                       !
! nsteps = number of MD steps                  !
! freq_print = how often to save the positions !
!       and momenta to the traj array--------- !
! dt = time step in fs ----------------------- !
! space_dim = dimension of space (1D, 2D, 3D)- !
! obtain_force = method for calculcating       !
!       force. Options include cp2k, model     !
! model_no = ?? Don't think used anymore       !
! temp = temperature (K) of the system         !
! -------------------------------------------- !
SUBROUTINE read_in_md(md)
implicit none
TYPE(md_info)     :: md
INTEGER           :: icount, REASON, tcount, i
CHARACTER(LEN=1)  ::skip_check
CHARACTER(LEN=3)  ::degree


open(101,file='md_run.inp')
read(101,*) md%nsteps
read(101,*) md%freq_print
read(101,*) md%dt
read(101,*) md%method 
read(101,*) md%ensemble
if (md%method /= 'pimd' .and. md%method /= 'rpmd') write(*,*)' this method is not supported'
if (md%ensemble == 'nvt') read(101,*) md%thermostat
read(101,*) md%space_dim
read(101,*) md%obtain_force
read(101,*) md%model_no
read(101,*) md%temp
read(101,*) md%propagation
read(101,*) md%init_vel
read(101,*) md%constraint
read(101,*) md%classical_lang
read(101,*) md%rp_lang
read(101,*) md%print_level
if (trim(adjustl(md%print_level)) /= 'low' .and. trim(adjustl(md%print_level)) /= 'high') &
    write(*,*)' print_levels are "high" and "low" '
close(101)

! -------------------------------------------- !
! WE ASSUME THAT TIME WAS READ IN IN FS. HERE  !
! WE CONVERT THIS TIME INTO ATOMIC UNITS. ---- !
! -------------------------------------------- !
!convert time into atomic units
md%dt = md%dt * fstota

md%kinetic_primitive = 0.0d0
md%kinetic_virial = 0.0d0
md%potential_virial = 0.0d0

! --------------------------------------------- !
! READS IN FILE FF_NONBONDED.DAT WHICH CONTAINS !
! ALL THE FORCE-FIELD PARAMETERS (SUCH AS LJ    !
! CONSTANTS AND ATOMIC CHARGES). THESE          !
! PARAMETERS ARE THEN STORED IN FF%NONBOND_* -- !
! VECTORS. ------------------------------------ !
! FF%NONBOND_ATOM = ATOM LABEL
! FF%NONBON_PARAM = (mass (amu), charge (au), LJ epsilon (K), LJ sigma (A)), 
! lang(thermostat)
! --------------------------------------------- ! 
if (md%constraint) then
open(101, file='md_constraint.inp')
read(101,*) md%fixed
if (md%fixed) then
  read(101,*) md%n_fixed_atoms                 !NEW EDIT (08/18)
  allocate(md%fixed_list(md%n_fixed_atoms,md%space_dim+1))
  do i = 1, md%n_fixed_atoms   
    read( 101, *) md%fixed_list(i,1), degree
    if (trim(adjustl(degree)) == 'x') then
      md%fixed_list(i,2:md%space_dim+1) = (/1,0,0/)
    else if (trim(adjustl(degree)) == 'y') then
      md%fixed_list(i,2:md%space_dim+1) = (/0,1,0/)
    else if (trim(adjustl(degree)) == 'z') then
      md%fixed_list(i,2:md%space_dim+1) = (/0,0,1/)
    else if (trim(adjustl(degree)) == 'xy') then
      md%fixed_list(i,2:md%space_dim+1) = (/1,1,0/)
    else if (trim(adjustl(degree)) == 'xz') then
      md%fixed_list(i,2:md%space_dim+1) = (/1,0,1/)
    else if (trim(adjustl(degree)) == 'yz') then
      md%fixed_list(i,2:md%space_dim+1) = (/0,1,1/)
    else if (trim(adjustl(degree)) == 'xyz') then
      md%fixed_list(i,2:md%space_dim+1) = (/1,1,1/)
    endif
    write(*,*)' checking fixed_list'
    write(*,*) md%fixed_list(i,1:md%space_dim+1)
  enddo
endif
read(101,*) md%centroid
if (md%centroid) then
  read(101,*) md%n_centroids 
  allocate(md%centroid_list(md%n_centroids,md%space_dim+1))
  do i = 1, md%n_centroids
    read(101,*)md%centroid_list(i,1),degree
    if (trim(adjustl(degree)) == 'x') then
      md%centroid_list(i,2:md%space_dim+1) = (/1,0,0/)
    else if (trim(adjustl(degree)) == 'y') then
      md%centroid_list(i,2:md%space_dim+1) = (/0,1,0/)
    else if (trim(adjustl(degree)) == 'z') then
      md%centroid_list(i,2:md%space_dim+1) = (/0,0,1/)
    else if (trim(adjustl(degree)) == 'xy') then
      md%centroid_list(i,2:md%space_dim+1) = (/1,1,0/)
    else if (trim(adjustl(degree)) == 'xz') then
      md%centroid_list(i,2:md%space_dim+1) = (/1,0,1/)
    else if (trim(adjustl(degree)) == 'yz') then
      md%centroid_list(i,2:md%space_dim+1) = (/0,1,1/)
    else if (trim(adjustl(degree)) == 'xyz') then
      md%centroid_list(i,2:md%space_dim+1) = (/1,1,1/)
    endif
    write(*,*)' checking centroid_list'
    write(*,*) md%centroid_list(i,1:md%space_dim+1)
  enddo
endif
read(101,*) md%virtual_site
if (md%virtual_site) then
  read(101,*) md%vs_type
  if (trim(adjustl(md%vs_type)) /= 'tip4p') write(*,*)' THESES VIRTUAL SITES ARE NOT SUPPORTED '
  read(101,*) md%n_virtual_sites
  allocate(md%virtual_list(md%n_virtual_sites))
  md%virtual_list(:) = 0
endif
else
  md%centroid = .false.
  md%fixed = .false.
  md%virtual_site = .false.
endif

write(*,*)' end of read in md'


END SUBROUTINE read_in_md
! --------------------------------------------- !

SUBROUTINE set_up_cp2k(para_var)
implicit none
TYPE(parallel_env)    :: para_var
character(len=25)     :: version_cp2k
integer               :: ierror
include 'mpif.h'

! ------------------------------------------------ !
! CALL cp2k_get_version(*,*) from the CP2K library !
! to confirm correct libary linking as well as --- !
! the version of CP2K. Current code was made with- !
! VERSION 7.0 ------------------------------------ !
! ------------------------------------------------ !
CALL cp2k_get_version(version_cp2k,LEN(version_cp2k))  
!write(*,*)' USING ', TRIM(version_cp2k)

! ------------------------------------------------ !
! CALL cp2k_init() from the CP2K library --------- !
! This subroutine initializes CP2K with MPI ------ !
! Must be paired with cp2k_finalize() ------------ !
! ???????????????????????????????????????????????? !
! NOT completely sure what this does, but it must- !
! Set up arrays and environments? ---------------- !
! ------------------------------------------------ !
!CALL MPI_INIT(ierror)
CALL cp2k_init()
call MPI_BARRIER(MPI_COMM_WORLD,ierror)

call MPI_COMM_SIZE(MPI_COMM_WORLD, para_var%size_rank, ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD, para_var%my_rank, ierror)

END SUBROUTINE set_up_cp2k


! -------------------------------------------------- !
! A subroutine that reads in the path-integral ----- !
! related variables from pimd.inp and creates the -- !
! path_integral structures necessary for each ring - !
! polymer. ----------------------------------------- !
! From pimd.inp                                      !
!       num_rp= number of individual ring polymers   !
!       num_beads = number of beads for that RP      !
!       particle_type = not quite sure what info this!
!               carries anymore                      !
!       kinetic_mass = mass of RP in kinetic term    !
!       harm_mass = mass of RP used in spring terms  !
!       charge = charge of the particle being        !
!               replaced by RP                       !
!       coord = how to read in the initial coord     !
!               of the RP                            !
!From system_parameters.inp                          !
!       initialize = how to initalize the RP coord.  !
!       debug_initialize = whether to print info     !
!               used to debug the mc_initalize       !
!       init_steps = number of mc steps to initalize !
!       init_print_steps = frequency to steps in mc  !
!               initalize scheme to print            !
!       mc_step = step size of MC step               !
!       propagation = method for propagation (single !
!               bead moves, normal, Cayle, etc.      !  
! -------------------------------------------------- !
SUBROUTINE pimd_initialize(path_integral, md)
implicit none
TYPE atom_type
  integer                                     ::nbead
  real(kind=8)                                ::update_mass
END TYPE atom_type
TYPE molecule_type
  character(len=160)                          ::psf_input     
  integer                                     ::natoms, total_beads, num_repeat
  TYPE(atom_type),allocatable,dimension(:)    ::atoms
  logical                                     ::first=.true.
END TYPE molecule_type

TYPE(path_integral_env) ::path_integral
TYPE(md_info)           ::md
integer                 ::i,j,icount, iatom, irp
integer                 ::ibead
character(len=6)        ::char_h

integer                 ::num_par
real(kind=8),dimension(3)   ::com_vel
real(kind=8)            ::kin_energy
TYPE(molecule_type),allocatable,dimension(:)     ::molecules
character(len=160)      ::psf_file 
character(len=20)        :: label
character(len=16)        :: label_long
integer                 ::psf_repeat, imolc


!For now, will not try any "fancy" reads. Modify this later.
! ------------------------------------------------------- !
! READ IN THE GENERAL INFORMATION ABOUT EACH OF THE RING- !
! POLYMERS IN THE SYSTEM. SEE NOTES IN PIMD.INP FOR MORE- !
! INFO ON EACH VARIABLE. -------------------------------- !
! ------------------------------------------------------- ! 
OPEN(101,FILE='rp_molecule_info.inp')
read(101,*) path_integral%num_molecules
allocate(molecules(path_integral%num_molecules))
molecules(:)%natoms = 0
molecules(:)%num_repeat = 0
molecules(:)%psf_input = ' '
irp = 0
path_integral%num_rp = 0
do while (irp < path_integral%num_molecules)
  irp = irp+1
  read(101,*) molecules(irp)%natoms, molecules(irp)%num_repeat, molecules(irp)%psf_input 
  !write(*,*)' natoms for molecule ', irp, ' is ', molecules(irp)%natoms
  allocate(molecules(irp)%atoms(molecules(irp)%natoms))
  molecules(irp)%atoms(:)%update_mass = -1.0
  do iatom = 1, molecules(irp)%natoms
    if (trim(adjustl(md%method)) == 'pimd') then
      read(101,*) molecules(irp)%atoms(iatom)%nbead, molecules(irp)%atoms(iatom)%update_mass
    else
      read(101,*) molecules(irp)%atoms(iatom)%nbead
    endif
  enddo
  molecules(irp)%total_beads = sum(molecules(irp)%atoms(:)%nbead)
  path_integral%num_rp = path_integral%num_rp + molecules(irp)%natoms*molecules(irp)%num_repeat
  do i = 2, molecules(irp)%num_repeat
    molecules(irp+i-1)%natoms = molecules(irp)%natoms
    allocate(molecules(irp+i-1)%atoms(molecules(irp)%natoms))
    molecules(irp+i-1)%atoms(:)%nbead = molecules(irp)%atoms(:)%nbead
    if (trim(adjustl(md%method)) == 'pimd') then
      molecules(irp+i-1)%atoms(:)%update_mass = molecules(irp)%atoms(:)%update_mass
    endif
    molecules(irp+i-1)%first = .false.
  enddo
  irp = irp+molecules(irp)%num_repeat-1
enddo
close(101)



allocate(path_integral%rp(path_integral%num_rp))
irp = 0
do i = 1, path_integral%num_molecules
  do j = 1, molecules(i)%natoms
    irp = irp + 1
    path_integral%rp(irp)%num_beads = molecules(i)%atoms(j)%nbead
    if (trim(adjustl(md%method)) == 'pimd') then
      path_integral%rp(irp)%kinetic_mass = molecules(i)%atoms(j)%update_mass*amutome
    endif
  enddo
enddo

open(101,file='RP_molset.pot')
irp = 0
do i = 1, SIZE(molecules)
  if (molecules(i)%first) then
    !write(*,*)' molecule ', i
    irp = irp+1
    if (irp/=1) read(101,*)
    !write(*,*)' irp is ', irp
    read(101,*)
    read(101,*) label, psf_repeat
    if (psf_repeat /= molecules(i)%num_repeat)  then
      write(*,*)' error in match RP_molset and rp_molecule_info.inp '
      STOP
    endif
    read(101,*)
    read(101,*) label_long, psf_file 
    path_integral%rp(irp)%psf_file = trim(adjustl(psf_file))
    path_integral%rp(irp)%molecule_atom_number = 1
    do j = 2, molecules(i)%natoms
      irp = irp+1
      !write(*,*)' irp is ', irp, ' as same molecule '
      path_integral%rp(irp)%psf_file = trim(adjustl(psf_file))
      path_integral%rp(irp)%molecule_atom_number = j
    enddo
    do j = 2, molecules(i)%num_repeat
      do iatom = 1, molecules(i)%natoms
        irp = irp+1
      !write(*,*)' irp is ', irp, ' as repeat molecules'
        path_integral%rp(irp)%psf_file = trim(adjustl(psf_file))
        path_integral%rp(irp)%molecule_atom_number = iatom
      enddo
    enddo
  endif
enddo
close(101)

if (irp /= path_integral%num_rp) then
  write(*,*)' matching of psf file/atom numbers does not match calcuclated num_rp'
  STOP
endif



path_integral%total_beads = 0
do i = 1, path_integral%num_rp
  !write(*,*) 'num beads for ', i, 'is ', path_integral%rp(i)%num_beads
  path_integral%total_beads = path_integral%total_beads + path_integral%rp(i)%num_beads
enddo
!write(*,*)' total number of beads are ', path_integral%total_beads
  
! -------------------------------------------------------- !
! ALLOCATES THE INITIAL POSITION ARRAYS FOR ALL BEAD ----- !
! AND CENTROID POSITIONS --------------------------------- !
! READS IN EITHER THE INITIAL COORDINATES OR INITIAL ----- !
! CENTROID POSITIONS FROM FILES RP_COORD.XYZ OR ---------- !
! -------------------------------------------------------- !
! -------------------------------------------------------- !
do i = 1, path_integral%num_rp
  ALLOCATE(path_integral%rp(i)%i_position(path_integral%rp(i)%num_beads,md%space_dim))
  path_integral%rp(i)%i_position = 0.0d0
  ALLOCATE(path_integral%rp(i)%i_centroid(md%space_dim))
  path_integral%rp(i)%i_centroid = 0.0d0
enddo



END SUBROUTINE pimd_initialize
! ------------------------------------------------------ !

! ---------------------------------------------- !
! Subroutine to initalize the velocoties of the- !
! beads for each RP. They are random numbers---- !
! generaed by a normal distribution to mimic the !
! Maxwell-Boltzmann with mean = 0 and ---------- !
! sigma = sqrt(1/kinetic_mass*betaN) in each --- !
! direction.------------------------------------ !
SUBROUTINE initialize_bead_velocities(path_integral, md, para_var)
implicit none

TYPE(path_integral_env)  ::path_integral
TYPE(md_info)            ::md
TYPE(parallel_env)       ::para_var
integer                  ::irp, ibead, i, j
real(kind=8)             ::rnd_vel, com
integer                  ::ierror, icount
include 'mpif.h'

if (para_var%my_rank == 0) then
do irp = 1, path_integral%num_rp

! ---------------------------------------------- !
! ALLOCATES THE INITALIZE_VELOCITY ARRAYS ------ !
! ---------------------------------------------- !
!  ALLOCATE(path_integral%rp(irp)%initialized_velocities(path_integral%rp(irp)%num_beads, &
!md%space_dim))
!  path_integral%rp(irp)%initialized_velocities = 0.0d0

    !write(*,*) 'initalizing bead velocities ', irp
  do ibead = 1, path_integral%rp(irp)%num_beads
    do i = 1, md%space_dim
! ---------------------------------------------- !
! CALL RANDOM GUASSIAN NUMBER GENERATOR FOR EACH !
! BEAD IN EACH DIMENSION OF SPACE -------------- !
! ---------------------------------------------- !
      CALL rnd_gauss(0.0d0, dsqrt(1.0d0/(path_integral%rp(irp)%beta * &
         path_integral%rp(irp)%kinetic_mass)), rnd_vel)
      path_integral%rp(irp)%initialized_velocities(ibead,i) = rnd_vel
    enddo
  enddo
enddo
if (md%fixed) then
  icount = 1
  do while (md%fixed_list(icount,1) <= path_integral%num_rp)
    write(*,*)' fixing rp velocities, atom ', md%fixed_list(icount,1)
    if (md%fixed_list(icount,2) == 1) then
      path_integral%rp(md%fixed_list(icount,1))%initialized_velocities(:,1) = 0.0d0
    endif
    if (md%fixed_list(icount,3) == 1) then
      path_integral%rp(md%fixed_list(icount,1))%initialized_velocities(:,2) = 0.0d0
    endif
    if (md%fixed_list(icount,4) == 1) then
      path_integral%rp(md%fixed_list(icount,1))%initialized_velocities(:,3) = 0.0d0
    endif
    icount = icount+1
    if (icount > md%n_fixed_atoms ) EXIT
  enddo
endif
if (md%centroid) then
  icount = 1
  do while (md%centroid_list(icount,1) <= path_integral%num_rp)
    !write(*,*)' fixing rp centroid velocities, atom ', md%centroid_list(icount,1)
    if (md%centroid_list(icount,2) == 1) then
      com = sum(path_integral%rp(md%centroid_list(icount,1))%initialized_velocities(:,1))/&
        dble(real(path_integral%rp(md%centroid_list(icount,1))%num_beads))
      path_integral%rp(md%centroid_list(icount,1))%initialized_velocities(:,1) = &
        path_integral%rp(md%centroid_list(icount,1))%initialized_velocities(:,1) - com
    endif
    if (md%centroid_list(icount,3) == 1) then
      com = sum(path_integral%rp(md%centroid_list(icount,1))%initialized_velocities(:,2))/&
        dble(real(path_integral%rp(md%centroid_list(icount,1))%num_beads))
      path_integral%rp(md%centroid_list(icount,1))%initialized_velocities(:,2) = &
        path_integral%rp(md%centroid_list(icount,1))%initialized_velocities(:,2) - com
    endif
    if (md%centroid_list(icount,4) == 1) then
      com = sum(path_integral%rp(md%centroid_list(icount,1))%initialized_velocities(:,3))/&
        dble(real(path_integral%rp(md%centroid_list(icount,1))%num_beads))
      path_integral%rp(md%centroid_list(icount,1))%initialized_velocities(:,3) = &
        path_integral%rp(md%centroid_list(icount,1))%initialized_velocities(:,3) - com
    endif
    icount = icount+1
    if (icount > md%n_centroids ) EXIT
  enddo
endif
  !open(222,file='rp-vel-2.dat')
  !do i = 1, path_integral%num_rp
  !  do j = 1, path_integral%rp(i)%num_beads
  !    write(222,'(3f10.6)') path_integral%rp(i)%initialized_velocities(j,:)
  !  enddo
  !enddo
  !close(222)
else
  do irp = 1, path_integral%num_rp
    path_integral%rp(irp)%initialized_velocities(:,:) = 0.0d0
  enddo
endif
call MPI_BARRIER(MPI_COMM_WORLD,ierror)
!write(*,*)' on rank ', para_var%my_rank, ' during initialization '
!write(*,*)' the velocities are ', path_integral%rp(1)%initialized_velocities(1,:)

END SUBROUTINE initialize_bead_velocities
! --------------------------------------------- !


! -------------------------------------------------- ! 
! Subroutine uses the Box-Mueller method to generate ! 
! a psedo-random number from the normal distribution ! 
! with a mean=x and a sigma/standard deviation=sig.- !
! -------------------------------------------------- !
SUBROUTINE rnd_gauss(x,sig,grnd)
implicit none

real(kind=8)          ::x,sig
real(kind=8)          ::v1,v2,w
real(kind=8)          ::grnd,z0
real(kind=8)          ::rnd1,rnd2

w=2.0d0

do
  call random_number(rnd1)
  call random_number(rnd2)

  v1 = 2.0d0*rnd1-1.0d0
  v2 = 2.0d0*rnd2-1.0d0
  
  w = v1*v1 + v2*v2

  if (w .lt. 1.d0) exit
enddo



z0 = v1*dsqrt(-2.0d0*log(w)/w)

grnd = z0*sig + x

END SUBROUTINE rnd_gauss
! ------------------------------------------------- !






SUBROUTINE cp2k_pimd_initialize(path_integral, classical, md, para_var)
! ------------------------------------------------ !
! Subroutine that reads in the atom positions and  !
! information from the input file (.xyz or .pdb) - !
! ------------------------------------------------ !
use iso_c_binding
implicit none
TYPE(path_integral_env)                ::path_integral
TYPE(classical_type)                   ::classical
TYPE(md_info)                           ::md
TYPE(parallel_env)                      ::para_var
INTEGER          ::i, j, ibead, irp, icount, i_vs
CHARACTER(LEN=10) ::hold_c
CHARACTER(LEN=9)  :: input_file, output_file, structure_file
INTEGER           ::hold_i, length
INTEGER           ::length_i,length_o,length_s
INTEGER           ::total_natoms
REAL(KIND=8),ALLOCATABLE,DIMENSION(:)    :: position_array
real(kind=8)     ::pot_energy
real(kind=8)     ::omega_n, omega_k

integer,allocatable,dimension(:)      ::repeat_mols
character(len=60),allocatable,dimension(:)    ::psf_mols
character(len=160)                   :: psf_format
integer                              :: file_check, natom_psf, iatom
character(len=5)                    :: label
character(len=16)                    :: label_long, char_check
integer                              :: i_dum, mol_num, nskip, imolc, num_mols
real(kind=8)                         :: atom_charge, atom_mass
character(len=160)                   :: mol_name, res_id, atom_id, atom_ff_type
integer                              :: ierror

include 'mpif.h'

! ------------------------------------------------ !
! READS IN INPUT FROM CLASSICAL.INP. THIS FILE --- !
! CONTAINS ANY INPUT RELEVANT TO THE CLASSICAL --- !
! VARIABLES. ------------------------------------- !
! ------------------------------------------------ !
write(*,*)' read in classical.inp '
open(101,file='file_info.inp')
read(101,*) classical%cp2k_input_file
read(101,*) classical%cp2k_output_file
read(101,*) classical%cp2k_structure_file
read(101,*) classical%atom_present
read(101,*) classical%cp2k_conn_file
close(101)


! ------------------------------------------------ !
! CREATES THE FORCE ENVIRONEMENT IN CP2K THAT ---- !
! IS BASED UPONN THE ABOVE INPUT FILE AND THE ---- !
! FORCE FIELD SPECIFIED IN THAT FILE. ------------ !
! ------------------------------------------------ !
!write(*,*) classical%cp2k_input_file
!write(*,*) classical%cp2k_structure_file
call MPI_BARRIER(MPI_COMM_WORLD,ierror)
CALL cp2k_create_force_env(classical%cp2k_env_flag,trim(classical%cp2k_input_file)//c_null_char, &
    trim(classical%cp2k_output_file)//c_null_char)
!CALL cp2k_create_force_env(classical%cp2k_env_flag,classical%cp2k_input_file,classical%cp2k_output_file)

! ------------------------------------------------ !
! USES THE CP2K FORCE_ENV TO READ IN THE NUMBER OF !
! TOTAL ATOMS (THE RP BEADS AND THEN THE CLASSICAL !
! ATOMS. THESE ARE THEN PARTITIONED INTO THE ----- !
! RELEVANT ARRAYS FOR STORAGE. ------------------- !
! ------------------------------------------------ !
!if (para_var%my_rank == 0) then
  CALL cp2k_get_natom(classical%cp2k_env_flag, total_natoms)
  write(*,*)'natoms is ', total_natoms

  CALL pimd_initialize(path_integral,md)
  

  classical%natoms = total_natoms - path_integral%total_beads
  !write(*,*)' total number of classical atoms is ', classical%natoms

  ALLOCATE(classical%i_position(classical%natoms,md%space_dim))
  ALLOCATE(position_array(total_natoms*md%space_dim))
  !write(*,*)' call cp2k get positions ', total_natoms*md%space_dim
CALL cp2k_get_positions(classical%cp2k_env_flag,position_array,total_natoms*md%space_dim)
!write(*,*)' done calling cp2k get positions '

icount = 0
do irp = 1, path_integral%num_rp
  !write(*,*)' path integral number ', irp
  do ibead = 1, path_integral%rp(irp)%num_beads
    icount = icount + 1
    !write(*,*)' icout is ', icount
    !write(*,*)' position index is ', (icount-1)*md%space_dim+1, icount*md%space_dim
    path_integral%rp(irp)%i_position(ibead,1:md%space_dim) = position_array( (icount-1)*md%space_dim+1:icount*md%space_dim)
  enddo
  do j = 1, md%space_dim
      !write(*,*)' path integral centroid dimension ', j
      path_integral%rp(irp)%i_centroid(j) = sum(path_integral%rp(irp)%i_position(:,j))/real(path_integral%rp(irp)%num_beads)
      !!write(*,*)' the initial calculcated i_centroid is ', path_integral%rp(irp)%i_centroid(j)/AtoBohr
  enddo
enddo
if (icount /= path_integral%total_beads) WRITE(*,*)" ERROR WITH PATH INTEGRAL POSITIONS "
do i = 1, classical%natoms
  classical%i_position(i,1:md%space_dim) = position_array( (i-1+path_integral%total_beads)*md%space_dim+1:&
     (i+path_integral%total_beads)*md%space_dim)
enddo

!write(*,*)' call cp2k calc energy force '
call MPI_BARRIER(MPI_COMM_WORLD,ierror)
CALL cp2k_calc_energy_force(classical%cp2k_env_flag)
CALL cp2k_get_potential_energy(classical%cp2k_env_flag, pot_energy)
!write(*,*)' initial potential energy is ', pot_energy

DEALLOCATE(position_array)

! ----------------------------------------------- !
! EDIT NEEDED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! READ IN THE ATOM IDENTIFIES/BEAD LABELS FROM -- !
! THE STRUCTURE FILE (XYZ/PDB) ------------------ !
! ----------------------------------------------- !
ALLOCATE(classical%ff_atom_charges(classical%natoms))
ALLOCATE(classical%atom_kinds(classical%natoms))
ALLOCATE(classical%mass(classical%natoms))
ALLOCATE(classical%force(classical%natoms,md%space_dim))
if (md%ensemble == 'nvt') ALLOCATE(classical%lang_atom_coeff(classical%natoms))

do irp = 1, path_integral%num_rp
  ALLOCATE(path_integral%rp(irp)%particle_label(path_integral%rp(irp)%num_beads))
  if (md%ensemble == 'nvt') then 
    ALLOCATE(path_integral%rp(irp)%lang_coeff(path_integral%rp(irp)%num_beads))
    path_integral%rp(irp)%lang_coeff = 0.0d0
  endif
enddo
classical%force = 0.d0

  psf_format = '(I10,T12,A7,T21,I8,T30,A7,T39,A6,T47,A6,T53,F10.6,T69,F12.8,T88,I1)'

if (classical%atom_present) then
open(101,file=trim(adjustl(classical%cp2k_conn_file)))
icount = 0
do
  read(101,*,iostat=file_check) label
  if (file_check > 0) then
    write(*,*)' something in classical read in went wrong '
    STOP
  else if (file_check < 0) then
    !write(*,*)' end of file '
    EXIT
  else
    icount = icount+1
  endif
enddo

REWIND(101)
!write(*,*)' icount is ', icount
if (mod(icount,5)/=0) then
  write(*,*)' the classical_molset file does not have the correct number of lines '
  STOP
endif

num_mols = icount/5
allocate(repeat_mols(num_mols))
allocate(psf_mols(num_mols))
!write(*,*)' num_mols is ', num_mols
do i = 1, num_mols
    read(101,*)
    read(101,*) label, repeat_mols(i)
    read(101,*)
    read(101,*) label_long, psf_mols(i)
    read(101,*)
enddo
close(101)

icount = 0
i_vs = 0
do i = 1, num_mols
  open(102, file=trim(adjustl(psf_mols(i))))
  read(102,*) char_check
  read(102,*) 
  read(102,*) nskip
  do j = 1, nskip
    read(102,*)
  enddo
  read(102,*)
  read(102,*) natom_psf
  do iatom = 1, natom_psf
    icount = icount + 1
    read(102,psf_format) i_dum, mol_name, mol_num, res_id, atom_id, atom_ff_type, atom_charge, atom_mass
    classical%atom_kinds(icount) = trim(adjustl(atom_ff_type))
    classical%mass(icount) = atom_mass * amutome
    classical%ff_atom_charges(icount) = atom_charge
    if (classical%atom_kinds(icount)(1:1) == 'M') then
      i_vs = i_vs + 1
      if (i_vs > md%n_virtual_sites) write(*,*)' VIRUTAL SITE INDEXING IS WRONG '
      md%virtual_list(i_vs) = icount
    endif
    !write(*,*) classical%atom_kinds(icount)
    !write(*,*) classical%mass(icount)
    !write(*,*) classical%ff_atom_charges(icount)
  enddo
  !write(*,*)' icount is ', icount
  !write(*,*)' mol number is ', i, ' out of ', num_mols
  !write(*,*)' completed natoms ', natom_psf
  !write(*,*)' now will repeat ', repeat_mols(i)
  do imolc = 2, repeat_mols(i)   !figure out this loop
    do iatom = 1, natom_psf
      icount = icount + 1
      classical%atom_kinds(icount) = classical%atom_kinds(icount-(imolc-1)*natom_psf)
      classical%mass(icount) = classical%mass(icount-(imolc-1)*natom_psf)
      classical%ff_atom_charges(icount) = classical%ff_atom_charges(icount-(imolc-1)*natom_psf)
      if (classical%atom_kinds(icount)(1:1) == 'M') then
        i_vs = i_vs + 1
        if (i_vs > md%n_virtual_sites) write(*,*)' VIRUTAL SITE INDEXING IS WRONG '
        md%virtual_list(i_vs) = icount
      endif
        !write(*,*) classical%atom_kinds(icount)
        !write(*,*) classical%mass(icount)
        !write(*,*) classical%ff_atom_charges(icount)
    enddo
  enddo
  close(102)
enddo


if (icount /= classical%natoms) then
  write(*,*)' issue with classical ff iteraction matching structure input '
  STOP
endif

endif


!open(111,file='classical_info_check.dat')
!do i = 1, classical%natoms
!  write(111,*) i, classical%atom_kinds(i), classical%mass(i)/amutome, classical%ff_atom_charges(i)
!enddo
!close(111)


if (md%ensemble == 'nvt') then
  classical%lang_atom_coeff(:) = 1.0d0/(md%classical_lang * fstota)
endif

! ------------------------------------------------------ !
! Path-integral force-fields --------------------------- !
! ------------------------------------------------------ !
do i = 1, path_integral%num_rp
  if (path_integral%rp(i)%molecule_atom_number == 1) then
    open(103,file=trim(adjustl(path_integral%rp(i)%psf_file)))
    !write(*,*)' opened ',trim(adjustl(path_integral%rp(i)%psf_file))
    read(103,*) char_check
    read(103,*) 
    read(103,*) nskip
    do j = 1, nskip
      read(103,*)
    enddo
    read(103,*)
    read(103,*) natom_psf
    !write(*,*)' natom_psf ', natom_psf
    read(103,psf_format) i_dum, mol_name, mol_num, res_id, atom_id, atom_ff_type, atom_charge, atom_mass
    if (i_dum /= 1) then
      write(*,*)' error in PI PSF read in for path integral ', i, ' with number ', &
	      path_integral%rp(i)%molecule_atom_number, 'and i_dum ', i_dum
      STOP
    endif

    path_integral%rp(i)%particle_label(:) = atom_ff_type

! ------------------------------------------------------- !
! READS IN VALUES FOR XX_CUTOFF IF AN ELECTRON RP IS ---- !
! PRESENT. ---------------------------------------------- !
! ------------------------------------------------------- !
    if (.NOT. (md%xx_present) .and. (trim(path_integral%rp(i)%particle_label(1)&
             (len(trim(path_integral%rp(i)%particle_label(1)))-1:)) == 'XX')) then
      md%xx_present = .true.
      open(111,file='xx_rcut.inp')
        read(111,*) md%xx_rcut
      close(111)
      md%xx_atom_num = i
      write(*,*)' electron is present ', md%xx_atom_num
    endif


    path_integral%rp(i)%harm_mass = atom_mass * amutome
    !path_integral%rp(i)%harm_mass = (atom_mass*real(path_integral%rp(i)%num_beads))*amutome
    path_integral%rp(i)%charge = atom_charge
    do ibead = 2, path_integral%rp(i)%num_beads
      read(103,*) i_dum
    enddo
    if (i_dum == natom_psf) then
      close(103)
      !write(*,*)' closed file '
    endif
    !write(*,*)' harm mass for ', i,' is ', path_integral%rp(i)%harm_mass/amutome
    !write(*,*)' kinetic mass for ', i,' is ', path_integral%rp(i)%kinetic_mass/amutome
  else
    read(103,psf_format) i_dum, mol_name, mol_num, res_id, atom_id, atom_ff_type, atom_charge, atom_mass
    !if ((i_dum-SUM(path_integral%rp(:i-1)%num_beads)+(i-1)) /= path_integral%rp(i)%molecule_atom_number) then
    !  write(*,*)' error in PI PSF read in for ', i_dum, SUM(path_integral%rp(:i-1)%num_beads), i
    !  STOP
    !endif

    path_integral%rp(i)%particle_label(:) = atom_ff_type
    path_integral%rp(i)%harm_mass = atom_mass * amutome
    !path_integral%rp(i)%harm_mass = (atom_mass*real(path_integral%rp(i)%num_beads))*amutome
    path_integral%rp(i)%charge = atom_charge
    do ibead = 2, path_integral%rp(i)%num_beads
      read(103,*) i_dum
    enddo

    if (i_dum == natom_psf) then
      close(103)
      !write(*,*)' closed file '
    endif
    !write(*,*)' harm mass for ', i,' is ', path_integral%rp(i)%harm_mass/amutome
    !write(*,*)' kinetic mass for ', i,' is ', path_integral%rp(i)%kinetic_mass/amutome

  endif

  if (trim(adjustl(md%method)) == 'rpmd') then
    path_integral%rp(i)%kinetic_mass = path_integral%rp(i)%harm_mass
  endif
    !write(*,*)' harm mass for ', i,' is ', path_integral%rp(i)%harm_mass/amutome
    !write(*,*)' kinetic mass for ', i,' is ', path_integral%rp(i)%kinetic_mass/amutome
enddo

!open(111,file='path_integral_info.dat')
!do i = 1, path_integral%num_rp
!    write(111,*) i, path_integral%rp(i)%particle_label, path_integral%rp(i)%kinetic_mass/amutome, &
!	    path_integral%rp(i)%harm_mass/amutome, path_integral%rp(i)%charge
!enddo
!close(111)

if (md%ensemble == 'nvt') then
  do irp = 1, path_integral%num_rp
    path_integral%rp(irp)%lang_coeff(1) = 1.0d0/(md%rp_lang * fstota)
    !write(*,*)' lang_pi for ', irp,' is ', path_integral%rp(irp)%lang_coeff(1)
  enddo
endif

! ------------------------------------------------------- !
! CALCULATES BETA_N VALUES FOR EACH OF THE RING-POLYMERS  !
! ------------------------------------------------------- !
classical%beta = (1.0d0/ (md%temp * kb) )
do i = 1, path_integral%num_rp
  path_integral%rp(i)%betaN = (1.0d0/(kb * md%temp)) / real(dble(path_integral%rp(i)%num_beads))
  path_integral%rp(i)%beta = (1.0d0/(kb * md%temp))
enddo

! ------------------------------------------------------- !
! INITALIZE THE RANDOM NUMBER GENERATOR. ---------------- !
! ------------------------------------------------------- !
!write(*,*)' init seed - BAJ error'
call init_rnd_seed()

! ------------------------------------------------------- !
! INITALIZE THE VELOCITIES OF THE BEADS AND ---- !
! ALLOCATES THE INITALIZE_VELOCITY ARRAYS ------ !
! ---------------------------------------------- !
do i = 1, path_integral%num_rp
  ALLOCATE(path_integral%rp(i)%initialized_velocities(path_integral%rp(i)%num_beads, &
md%space_dim))
  path_integral%rp(i)%initialized_velocities = 0.0d0
enddo

if (md%init_vel .eqv. .true.) then
  CALL initialize_bead_velocities(path_integral,md, para_var)

  !open(222,file='rp-vel.dat')
  !do i = 1, path_integral%num_rp
  !  do j = 1, path_integral%rp(i)%num_beads
  !    write(222,'(3f10.6)') path_integral%rp(i)%initialized_velocities(j,:)
  !  enddo
  !enddo
  !close(222)
else
  open(222,file='rp-vel.dat')
  do i = 1, path_integral%num_rp
    do j = 1, path_integral%rp(i)%num_beads
      read(222,*) path_integral%rp(i)%initialized_velocities(j,:)
    enddo
  enddo
  close(222)
endif

! ---------------------------------------------- !
! IF NVT SIMULATION, NEED TO DETERMINE THE ----- !
! LANG COEFFICIENT (BASED ON THE ONE INPUT ----- !
! PARAMETER ------------------------------------ !
! ---------------------------------------------- !
if (md%ensemble == 'nvt') then
  do irp = 1, path_integral%num_rp
    !path_intgral%rp(irp)%lang_coeff(1) = 1.0d0
    omega_n = real(path_integral%rp(irp)%num_beads)/path_integral%rp(irp)%beta
    !omega_n = (dble(real(path_integral%rp(irp)%num_beads))/path_integral%rp(irp)%beta)&
    !    *dsqrt(path_integral%rp(irp)%harm_mass/path_integral%rp(irp)%kinetic_mass)
    do ibead = 2, path_integral%rp(irp)%num_beads
      omega_k =2.0d0*omega_n*dsin( (ibead-1)*pi/real(path_integral%rp(irp)%num_beads) )
      path_integral%rp(irp)%lang_coeff(ibead) = 2.0d0 * omega_k
    enddo
  enddo
endif


! ----------------------------------------------- !
! INITIALIZE THE CLASSICAL AND RP VELOCITIES ---- !
! USING THE ------------------------------------- !
! RANDOM NUMBER GENERATOR ON A MAXWELL-BOLTZMANN  !
! DISTRIBUTION. 
! CALCULATE THE INITIAL MOMENTA ----------------- !
! ----------------------------------------------- !

ALLOCATE(classical%i_velocity(classical%natoms,md%space_dim))
classical%i_velocity = 0.0d0
CALL calc_constraints(classical, path_integral, md)
if (md%init_vel .eqv. .true. ) then 
  CALL initialize_classical_velocities(classical,path_integral,md, para_var)
  !open(222,file='atom-vel-1.dat')
  !do i = 1, classical%natoms
  !  write(222,'(3d20.6)') classical%i_velocity(i,:)
  !enddo
  !close(222)
  !write(*,*) 'rescale vel '
  CALL rescale_velocities(classical, path_integral, md, para_var)
  write(*,*) 'classical vel '
  open(222,file='atom-vel.dat')
  do i = 1, classical%natoms
    write(222,'(3d20.6)') classical%i_velocity(i,:)
  enddo
  close(222)

  !write(*,*)' rp vel '
  open(222,file='rp-vel.dat')
  do i = 1, path_integral%num_rp
    do j = 1, path_integral%rp(i)%num_beads
      write(222,'(3f20.15)') path_integral%rp(i)%initialized_velocities(j,:)
    enddo
  enddo
  close(222)


else
  open(222,file='atom-vel.dat')
  do i = 1, classical%natoms
    read(222,*) classical%i_velocity(i,:)
  enddo
  close(222)
endif
! ---------------------------------------------- !



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




classical%potential = 0.0d0
classical%kinetic = 0.0d0




do i = 1, path_integral%num_rp
  ALLOCATE(path_integral%rp(i)%initialized_positions(path_integral%rp(i)%num_beads,md%space_dim))
  do j = 1, path_integral%rp(i)%num_beads
    path_integral%rp(i)%initialized_positions(j,:) = path_integral%rp(i)%i_position(j,:)
  enddo
enddo



END SUBROUTINE cp2k_pimd_initialize



SUBROUTINE initialize_classical_velocities(classical,path_integral,md, para_var)
implicit none

TYPE(classical_type) ::classical
TYPE(md_info)        ::md
TYPE(path_integral_env)           ::path_integral
TYPE(parallel_env)   ::para_var
real(kind=8)         ::rnd_vel
integer              ::iatom,i, icount
include 'mpif.h'

!ALLOCATE(classical%i_velocity(classical%natoms,md%space_dim))
!classical%i_velocity = 0.0d0


if (para_var%my_rank == 0) then
do iatom = 1, classical%natoms
  !write(*,*)' initialize classical ', iatom
  do i = 1, md%space_dim
    CALL rnd_gauss(0.0d0,dsqrt(1.0d0/( classical%beta * classical%mass(iatom) )), rnd_vel)
    classical%i_velocity(iatom,i) = rnd_vel
    write(*,*)' rnd vel is ', rnd_vel
    write(*,*)' classical%beta is ', classical%beta
    write(*,*)' classical%mass is ', classical%mass(iatom)
  enddo
  open(222,file='atom-vel-2.dat')
  do i = 1, classical%natoms
    write(222,'(3d20.6)') classical%i_velocity(i,:)
  enddo
  close(222)
enddo
if (md%fixed) then
  icount = md%n_fixed_atoms
  do while (md%fixed_list(icount,1) > path_integral%num_rp)
    write(*,*)' fixing classical atom ', md%fixed_list(icount,1)
    if (md%fixed_list(icount,2) == 1) then
      classical%i_velocity(md%fixed_list(icount,1)-path_integral%num_rp,1) = 0.0d0
    endif
    if (md%fixed_list(icount,2) == 1) then
      classical%i_velocity(md%fixed_list(icount,1)-path_integral%num_rp,2) = 0.0d0
    endif
    if (md%fixed_list(icount,2) == 1) then
      classical%i_velocity(md%fixed_list(icount,1)-path_integral%num_rp,3) = 0.0d0
    endif
    if (icount == 1) EXIT
    icount = icount - 1
  enddo
endif
  open(222,file='atom-vel-1.dat')
  do i = 1, classical%natoms
    write(222,'(3d20.6)') classical%i_velocity(i,:)
  enddo
  close(222)
else
  classical%i_velocity(:,:) = 0.0d0
endif

END SUBROUTINE initialize_classical_velocities


SUBROUTINE calc_constraints(classical, path_integral, md)
implicit none
TYPE(classical_type)              ::classical
TYPE(path_integral_env)           ::path_integral
TYPE(md_info)                     ::md
integer                           ::i,j



!CALCULATE THE NUMBER OF CONTRAINTS IN THE SYSTEM!
write(*,*)' calc_contraints' 
md%n_constraints = 0 
md%x_cnst = 0
md%y_cnst = 0
md%z_cnst = 0
if (md%fixed) then
  do i = 1, path_integral%num_rp
    do j = 1, md%n_fixed_atoms
      if (md%fixed_list(j,1)==i) then 
        md%n_constraints = md%n_constraints + sum(md%fixed_list(j,2:md%space_dim+1))&
          *path_integral%rp(i)%num_beads
        write(*,*)' RP fixed constraint ', md%fixed_list(j,1)
        md%x_cnst = md%x_cnst + md%fixed_list(j,2)*path_integral%rp(i)%num_beads
        md%y_cnst = md%y_cnst + md%fixed_list(j,3)*path_integral%rp(i)%num_beads
        md%z_cnst = md%z_cnst + md%fixed_list(j,4)*path_integral%rp(i)%num_beads
      endif
    enddo
  enddo
  write(*,*) 'num path integral ', path_integral%num_rp
  do i = 1, classical%natoms
    do j = 1, md%n_fixed_atoms
      write(*,*) 'if state', md%fixed_list(j,1), path_integral%num_rp + i
      if (md%fixed_list(j,1)==path_integral%num_rp+i) then 
        md%n_constraints = md%n_constraints + &
          sum(md%fixed_list(j,2:md%space_dim+1))
        write(*,*)' classical constraint ', md%fixed_list(j,1)
        md%x_cnst = md%x_cnst + md%fixed_list(j,2)
        md%y_cnst = md%y_cnst + md%fixed_list(j,3)
        md%z_cnst = md%z_cnst + md%fixed_list(j,4)
      endif
    enddo
  enddo
endif
if (md%centroid) then
  do i = 1, md%n_centroids
    md%n_constraints = md%n_constraints + sum(md%centroid_list(i,2:md%space_dim+1)) 
      md%x_cnst = md%x_cnst + md%centroid_list(i,2)*path_integral%rp(i)%num_beads
      md%y_cnst = md%y_cnst + md%centroid_list(i,3)*path_integral%rp(i)%num_beads
      md%z_cnst = md%z_cnst + md%centroid_list(i,4)*path_integral%rp(i)%num_beads
  enddo
endif
if (md%virtual_site) then
  md%n_constraints = md%n_constraints + md%n_virtual_sites*real(dble(md%space_dim))
      md%x_cnst = md%x_cnst + md%n_virtual_sites
      md%y_cnst = md%y_cnst + md%n_virtual_sites
      md%z_cnst = md%z_cnst + md%n_virtual_sites
endif

write(*,*)' number of constraints in simulation ', md%n_constraints
if (classical%natoms+path_integral%total_beads > 2 ) then
  md%n_constraints = md%n_constraints + 6
else
  md%n_constraints = md%n_constraints + 5
endif
write(*,*)' number of constraints in simulation ', md%n_constraints

END SUBROUTINE calc_constraints



SUBROUTINE rescale_velocities(classical, path_integral, md, para_var)
implicit none

TYPE(classical_type)              ::classical
TYPE(path_integral_env)           ::path_integral
TYPE(md_info)                     ::md
TYPE(parallel_env)                ::para_var
integer                           ::i,j, k
real(kind=8),dimension(3)         ::com_vel
integer                           ::num_par
real(kind=8)                      ::kin_energy, lambda, instant_temp
integer                           ::icount, irp
real(kind=8),allocatable,dimension(:)    ::vel_array
real(kind=8)                      ::xmin, xmax, delx
integer                           ::xnum
integer, allocatable,dimension(:) ::his
integer                           ::ierror, root
include 'mpif.h'



!CALCULATE THE NUMBER OF CONTRAINTS IN THE SYSTEM!
write(*,*)' rescale velocities' 
!md%n_constraints = 0 
!x_cnst = 0
!y_cnst = 0
!z_cnst = 0
!if (md%fixed) then
!  do i = 1, path_integral%num_rp
 !   do j = 1, md%n_fixed_atoms
!      if (md%fixed_list(j,1)==i) then 
!        md%n_constraints = md%n_constraints + sum(md%fixed_list(j,2:md%space_dim+1))&
!          *path_integral%rp(i)%num_beads
!        write(*,*)' RP fixed constraint ', md%fixed_list(j,1)
!        x_cnst = x_cnst + md%fixed_list(j,2)*path_integral%rp(i)%num_beads
!        y_cnst = y_cnst + md%fixed_list(j,3)*path_integral%rp(i)%num_beads
!        z_cnst = z_cnst + md%fixed_list(j,4)*path_integral%rp(i)%num_beads
!      endif
!    enddo
!  enddo
!  write(*,*) 'num path integral ', path_integral%num_rp
!  do i = 1, classical%natoms
!    do j = 1, md%n_fixed_atoms
!      write(*,*) 'if state', md%fixed_list(j,1), path_integral%num_rp + i
!      if (md%fixed_list(j,1)==path_integral%num_rp+i) then 
!        md%n_constraints = md%n_constraints + &
!          sum(md%fixed_list(j,2:md%space_dim+1))
!        write(*,*)' classical constraint ', md%fixed_list(j,1)
!        x_cnst = x_cnst + md%fixed_list(j,2)
!        y_cnst = y_cnst + md%fixed_list(j,3)
!        z_cnst = z_cnst + md%fixed_list(j,4)
!      endif
!    enddo
!  enddo
!endif
!if (md%centroid) then
!  do i = 1, md%n_centroids
!    md%n_constraints = md%n_constraints + sum(md%centroid_list(i,2:md%space_dim+1)) 
!      x_cnst = x_cnst + md%centroid_list(i,2)*path_integral%rp(i)%num_beads
!      y_cnst = y_cnst + md%centroid_list(i,3)*path_integral%rp(i)%num_beads
!      z_cnst = z_cnst + md%centroid_list(i,4)*path_integral%rp(i)%num_beads
!  enddo
!endif
!if (md%virtual_site) then
!  md%n_constraints = md%n_constraints + md%n_virtual_sites*real(dble(md%space_dim))
!      x_cnst = x_cnst + md%n_virtual_sites
!      y_cnst = y_cnst + md%n_virtual_sites
!      z_cnst = z_cnst + md%n_virtual_sites
!endif

!if (x_cnst + y_cnst+z_cnst /= md%n_constraints) write(*,*)' constraints do not add up '
!write(*,*)' number of constraints in simulation ', md%n_constraints
!!Add 6 constraints for system if num particles > 2
!if (classical%natoms+path_integral%total_beads > 2 ) then
!  md%n_constraints = md%n_constraints + 6
!else
!  md%n_constraints = md%n_constraints + 5
!endif
!write(*,*)' number of constraints in simulation ', md%n_constraints

call MPI_BARRIER(MPI_COMM_WORLD,ierror)
if (para_var%my_rank == 0) then
if ( (md%init_vel .eqv. .true.) ) then

write(*,*)' allocate '
!write(*,*)' classical%natoms ', classical%natoms
!write(*,*)' path_integral%total_beads ', path_integral%total_beads
!write(*,*)' size of path_integral ', SIZE(path_integral%rp)
!allocate(vel_array(md%space_dim*(classical%natoms+path_integral%total_beads)))
!icount = 0
!do irp = 1, path_integral%num_rp
!  do i = 1, path_integral%rp(irp)%num_beads
!    icount = icount + 1
!    vel_array( (icount-1)*md%space_dim+1: icount*md%space_dim) = &
!	    path_integral%rp(irp)%initialized_velocities(i,:)
!  enddo
!enddo
!if (icount /= path_integral%total_beads) write(*,*)' ERROR IN VEL-ARRAY '
!write(*,*)' done pi'

!do i = 1, classical%natoms
!  vel_array( (i-1+path_integral%total_beads)*md%space_dim+1: (i+path_integral%total_beads)*md%space_dim ) = &
!    classical%i_velocity(i,:)
!enddo
!write(*,*)' done classical'


!xmin = -0.004
!xmax = 0.004
!xnum = 50
!write(*,*)' call his'
!allocate(his(0:xnum))
!CALL histogram(path_integral%total_beads*md%space_dim, vel_array(1:path_integral%total_beads*md%space_dim), xmin, xmax, xnum, his) 

!open(199, file='his_pi_before.out')
!delx = (xmax-xmin)/real(xnum-1)
!do i = 0, xnum
!  write(199,*) i, xmin + (i-1)*delx, his(i)
!enddo
!close(199)


  com_vel = 0.0d0
  do i = 1, path_integral%num_rp
        do j = 1, path_integral%rp(i)%num_beads
          com_vel(:) = com_vel(:) + path_integral%rp(i)%initialized_velocities(j,:)
        enddo
  enddo
  do i = 1, classical%natoms
      com_vel(:) = com_vel(:) + classical%i_velocity(i,:)
  enddo
  if (.not. md%constraint) then
    com_vel = com_vel/dble(real(path_integral%total_beads + classical%natoms))
  else
     com_vel(1) = com_vel(1)/(dble(real(path_integral%total_beads + classical%natoms - md%x_cnst)))
     com_vel(2) = com_vel(2)/(dble(real(path_integral%total_beads + classical%natoms - md%y_cnst)))
     com_vel(3) = com_vel(3)/(dble(real(path_integral%total_beads + classical%natoms - md%z_cnst)))
  endif
  
  write(*,*)' com_vel is ', com_vel

! Zero-out the com velocity
  do i = 1, path_integral%num_rp
    if (md%centroid .and. md%fixed) then
      if (.not. ANY(md%centroid_list(:,1)==i) .and. .not. ANY(md%fixed_list(:,1)==i) ) then
        !write(*,*)' rp centroid and fixed', i
        do j = 1, path_integral%rp(i)%num_beads
          path_integral%rp(i)%initialized_velocities(j,:) = path_integral%rp(i)%initialized_velocities(j,:) - com_vel(:)
        enddo
      endif
    else if (md%centroid .and. .not. md%fixed) then 
        if (.not. ANY(md%centroid_list(:,1)==i) ) then
        !write(*,*)' rp centroid', i
          do j = 1, path_integral%rp(i)%num_beads
            path_integral%rp(i)%initialized_velocities(j,:) = path_integral%rp(i)%initialized_velocities(j,:) - com_vel(:)
          enddo
        endif
    else if (.not. md%centroid .and. md%fixed) then 
        if (.not. ANY(md%fixed_list(:,1)==i) ) then
        !write(*,*)' rp fixed', i
          do j = 1, path_integral%rp(i)%num_beads
            path_integral%rp(i)%initialized_velocities(j,:) = path_integral%rp(i)%initialized_velocities(j,:) - com_vel(:)
          enddo
        endif
    else if (.not. md%fixed .and. .not. md%centroid) then
        !write(*,*)' rp not fixed not centroid ', i
        do j = 1, path_integral%rp(i)%num_beads
          path_integral%rp(i)%initialized_velocities(j,:) = path_integral%rp(i)%initialized_velocities(j,:) - com_vel(:)
        enddo
    endif
  enddo
  do i = 1, classical%natoms
    if (md%fixed) then
      if (.not. ANY(md%fixed_list(:,1)==i+path_integral%num_rp)) then
        write(*,*)' classical fixed'
        classical%i_velocity(i,:) = classical%i_velocity(i,:) - com_vel(:)
      endif
    else
        write(*,*)' classical not fixed'
      classical%i_velocity(i,:) = classical%i_velocity(i,:) - com_vel(:)
    endif
  enddo

!Check that it worked.
  com_vel = 0.0d0
  do i = 1, path_integral%num_rp
    do j = 1, path_integral%rp(i)%num_beads
      com_vel(:) = com_vel(:) + path_integral%rp(i)%initialized_velocities(j,:)
    enddo
  enddo
  do i = 1, classical%natoms
    com_vel(:) = com_vel(:) + classical%i_velocity(i,:)
  enddo
  write(*,*)' com_vel is ', com_vel



!Scale the velocities to match temperature

  kin_energy = 0.0d0
  do i = 1, classical%natoms
    !write(*,*)' classical velocities ', classical%i_velocity(i,:), classical%mass(i)
    kin_energy = kin_energy + 0.5d0 * dot_product(classical%i_velocity(i,:),classical%i_velocity(i,:)) * classical%mass(i)
  enddo
  do i = 1, path_integral%num_rp
      !write(*,*)' PI mass ', path_integral%rp(i)%kinetic_mass
    do j = 1, path_integral%rp(i)%num_beads
      !write(*,*)' path_integral velocities ', path_integral%rp(i)%initialized_velocities(j,:)
      kin_energy = kin_energy + 0.5d0*dot_product(path_integral%rp(i)%initialized_velocities(j,:), &
        path_integral%rp(i)%initialized_velocities(j,:))*&
        path_integral%rp(i)%kinetic_mass
    enddo
  enddo

  write(*,*)' md%n_constrains is ', md%n_constraints
  if (classical%natoms+path_integral%total_beads <=2 ) then
    instant_temp = kin_energy/((real(md%space_dim)*(classical%natoms+path_integral%total_beads)-5)*kb/2.0d0)
  else
    instant_temp = kin_energy/((real(md%space_dim)*(classical%natoms+path_integral%total_beads)-(md%n_constraints))*kb/2.0d0)
  endif
  !write(*,*)' kin_energy is ', kin_energy
  lambda = dsqrt(md%temp/instant_temp)

  do i = 1, path_integral%num_rp
    do j = 1, path_integral%rp(i)%num_beads
      path_integral%rp(i)%initialized_velocities(j,:) = path_integral%rp(i)%initialized_velocities(j,:) * lambda
    enddo
  enddo
  do i = 1, classical%natoms
    classical%i_velocity(i,:) = classical%i_velocity(i,:) * lambda
  enddo
endif
endif

root = 0
call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  do i = 1, path_integral%num_rp
    do j = 1, path_integral%rp(i)%num_beads
      !do k = 1, md%space_dim
      call MPI_bcast(path_integral%rp(i)%initialized_velocities(j,1:md%space_dim), md%space_dim, MPI_DOUBLE_PRECISION,root, &
        MPI_COMM_WORLD, ierror)
      !write(*,*)' value ', path_integral%rp(i)%initialized_velocities(j,:)
      !enddo
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    enddo
  enddo
  do i = 1, classical%natoms
    call MPI_bcast(classical%i_velocity(i,:), md%space_dim, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror) 
  enddo

!do i = 0, para_var%size_rank-1
!  if (i == para_var%my_rank) then
!    write(*,*)' on rank ', para_var%my_rank, ' after rescaling '
!    write(*,*)' the velocities are ', path_integral%rp(1)%initialized_velocities(1,:)
!  endif
!  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
!enddo


END SUBROUTINE rescale_velocities
! --------------------------------------------------- !


SUBROUTINE histogram(nsize, val, xmin, xmax, xnum, his)
implicit none

integer    ::nsize, xnum, ind, i
integer, dimension(0:xnum)   :: his
real(kind=8)                 ::xmin, xmax, delx
real(kind=8), dimension(nsize)   ::val

delx = (xmax-xmin)/real(xnum-1)

his = 0
do i = 1, nsize
  ind = int((val(i) - xmin)/delx )
  if (ind < 0 .or. ind > xnum) then
    write(*,*) ind, ' is outside the histogram range '
  else
    his(ind) = his(ind) + 1
  endif
enddo

END SUBROUTINE histogram




! ----------------------------------------------------- !
! Subroutine to initialize the seed for the random ---- !
! number generator. ----------------------------------- !
! ----------------------------------------------------- !
SUBROUTINE init_rnd_seed()
implicit none

integer     ::i,n,clock
integer,dimension(:),allocatable    ::seed

call random_seed(size=n)
allocate(seed(n))

call system_clock(count=clock)

seed = clock + 37 * (/ (i-1,i=1,n) /)
call random_seed(put=seed)

deallocate(seed)

END SUBROUTINE init_rnd_seed
! ----------------------------------------------------- !

END MODULE input_setup
