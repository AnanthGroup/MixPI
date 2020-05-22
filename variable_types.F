MODULE variable_types
implicit none

REAL(KIND=8),PARAMETER, PUBLIC        :: pi=DACOS(-1.D0)
REAL(KIND=8),PARAMETER, PUBLIC        :: kb=3.1668114e-6       !Eh/K
REAL(KIND=8),PARAMETER, PUBLIC        :: KtoH=3.1668114e-6
REAL(KIND=8),PARAMETER, PUBLIC        :: AtoBohr=1.8897259886d0  !Bohr/Angstrom
REAL(KIND=8),PARAMETER, PUBLIC        :: amutome = (1.0e-3)/( 6.022e23 * 9.1093837015e-31)   ! me/(kg/mol)
REAL(KIND=8),PARAMETER, PUBLIC        :: fstota = 1.0d0/(2.418884326e-2)      !t_h/fs

PUBLIC   ::path_integral_env
PUBLIC   ::ring_polymer_type

!-------------------------------------------------!
!A variable type that contains the information for!
!the path-integral atom types. -------------------!
TYPE ring_polymer_type
  INTEGER                               :: num_beads       !number of beads
  REAL(KIND=8)                          :: kinetic_mass, harm_mass !fictious mass used in dynamics/actual mass of complete particle
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:)             :: i_centroid      !initial location of centroid
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)           :: i_position
  REAL(KIND=8)                          :: charge          !charge of complete particle
  REAL(KIND=8),DIMENSION(2)             :: lj_param
  CHARACTER(LEN=100)                    :: particle_type    !atom, special (to indicate atom search or not)
  LOGICAL                               :: coord           !specifies if have a coordinate file for the beads
  REAL(KIND=8)                          :: betaN
  REAL(KINd=8)                          :: beta
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)       ::x_old,x_new
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)     ::init_traj
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)       ::initialized_positions      !the coordinates to begin MD run from
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)       ::initialized_velocities
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)     ::traj
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)       ::x
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)       ::p
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)       ::force
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)       ::norm_p, norm_x
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)       ::rot_matrix
END TYPE ring_polymer_type
  
TYPE path_integral_env
  INTEGER                                            :: num_rp
  LOGICAL                                            :: debug_initialize, initialize
  INTEGER                                            :: init_nsteps, init_print_nsteps
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:)              :: kinetic, potential, interaction_pot, inter_rp_pot
  REAL(KIND=8)                                       :: mc_step
  TYPE(ring_polymer_type),ALLOCATABLE,DIMENSION(:)               :: rp
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:)              :: com_pos, com_vel
  CHARACTER(LEN=20)                                  :: propagation     !none=single bead moves, normal=normal mode, Cayley 
END TYPE path_integral_env

TYPE md_info
  INTEGER                             :: nsteps
  INTEGER                             :: freq_print
  REAL(KIND=8)                        :: dt
  REAL(KIND=8)                        :: temp
  INTEGER                             :: space_dim       !dof in space (1d,2d,3d)
  CHARACTER(LEN=100)                                 :: obtain_force    !cp2k, model - indicates source of the external potential
  REAL(KIND=8)                                       :: model_no        !if model, indicate which model to use
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:)              :: total_energy
END TYPE md_info

TYPE classical_type
  INTEGER                             :: cp2k_env_flag
  INTEGER                             :: natoms 
  REAL(KIND=8)                        :: beta
  CHARACTER(LEN=150)                  :: cp2k_input_file,cp2k_output_file
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)        ::i_position
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:)          ::ff_atom_charges
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)        ::ff_atom_lj
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)          ::force
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)        ::i_velocity
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:)          ::mass
  CHARACTER(LEN=2),ALLOCATABLE,DIMENSION(:)      ::atom_kinds
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)        ::x,p
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)      ::traj
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:)          ::kinetic, potential
END TYPE classical_type

END MODULE variable_types
