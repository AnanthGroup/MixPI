MODULE variable_types
implicit none

REAL(KIND=8),PARAMETER, PUBLIC        :: pi=DACOS(-1.D0)
REAL(KIND=8),PARAMETER, PUBLIC        :: c_light=299792458.0    !speed of light in m/s
REAL(KIND=8),PARAMETER, PUBLIC        :: rhydberg=10973731.568527 !rhydberg constant in 1/m
REAL(KIND=8),PARAMETER, PUBLIC        :: h_planck=6.62606896E-34 !planck constant J*s
REAL(KIND=8),PARAMETER, PUBLIC        :: joule=2.0d0*rhydberg*h_planck*c_light !a.u. to joule
REAL(KIND=8),PARAMETER, PUBLIC        :: kb=1.3806504E-23/joule !3.1668114e-6       !Eh/K
REAL(KIND=8),PARAMETER, PUBLIC        :: KtoH=kb !3.1668114e-6
REAL(KIND=8),PARAMETER, PUBLIC        :: AtoBohr = 1.0d0/(0.52917720859d0) !AtoBohr=1.8897259886d0  !Bohr/Angstrom
REAL(KIND=8),PARAMETER, PUBLIC        :: amutome = 1.660538782e-27/9.10938215e-31!(1.0e-3)/( 6.022e23 * 9.1093837015e-31)   ! me/(kg/mol)
REAL(KIND=8),PARAMETER, PUBLIC        :: fstota = 4.0d0*pi*rhydberg*c_light/(1.0e15) !1.0d0/(2.418884326e-2)      !t_h/fs
REAL(KIND=8),PARAMETER, PUBLIC        :: kjmoltoH = 0.00038088d0

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
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: lang_coeff
  LOGICAL                               :: bond
  CHARACTER(LEN=100)                    :: particle_type    !atom, special (to indicate atom search or not)
  CHARACTER(LEN=10),ALLOCATABLE,DIMENSION(:)                      :: particle_label
  REAL(KIND=8)                          :: betaN
  REAL(KINd=8)                          :: beta
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)       ::x_old,x_new
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)     ::init_traj
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)       ::initialized_positions      !the coordinates to begin MD run from
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)       ::initialized_velocities
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)     ::traj
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)       ::x
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)       ::p
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:)       ::centroid
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)       ::force
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)       ::norm_p, norm_x
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)       ::rot_matrix
  CHARACTER(LEN=160)                            ::psf_file
  INTEGER                                       ::molecule_atom_number
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)              :: com_pos, com_vel
END TYPE ring_polymer_type
  
TYPE path_integral_env
  INTEGER                                            :: num_rp
  INTEGER                                            :: num_molecules
  INTEGER                                            :: total_beads
  REAL(KIND=8)                                       :: kinetic, potential
  TYPE(ring_polymer_type),ALLOCATABLE,DIMENSION(:)               :: rp
END TYPE path_integral_env

TYPE md_info
  INTEGER                             :: nsteps
  INTEGER                             :: freq_print
  REAL(KIND=8)                        :: dt
  REAL(KIND=8)                        :: temp
  INTEGER                             :: space_dim       !dof in space (1d,2d,3d)
  CHARACTER(LEN=100)                                 :: obtain_force    !cp2k, model - indicates source of the external potential
  CHARACTER(LEN=3)                                   :: ensemble
  CHARACTER(LEN=5)                                   :: thermostat
  CHARACTER(LEN=4)                                   :: method
  REAL(KIND=8)                                       :: model_no        !if model, indicate which model to use
  REAL(KIND=8)                                       :: kinetic_primitive, potential_virial, kinetic_virial
  CHARACTER(LEN=20)                                  :: propagation     !none=single bead moves, normal=normal mode, Cayley 
  LOGICAL                                            :: init_vel
  REAL(KIND=8)                                       :: classical_lang, rp_lang
  LOGICAL                                            :: constraint
  INTEGER                                            :: n_constraints
  LOGICAL                                            :: fixed
  INTEGER                                            :: n_fixed_atoms
  LOGICAL                                            :: centroid
  INTEGER                                            :: n_centroids
  LOGICAL                                            :: virtual_site
  CHARACTER(LEN=10)                                  :: vs_type
  INTEGER                                            :: n_virtual_sites
  INTEGER, ALLOCATABLE,DIMENSION(:,:)                :: fixed_list
  INTEGER, ALLOCATABLE,DIMENSION(:,:)                :: centroid_list
  INTEGER, ALLOCATABLE,DIMENSION(:)                  :: virtual_list
  CHARACTER(LEN=6)                                   ::print_level
  LOGICAL                                            ::xx_present = .false.
  REAL(kind=8)                                       ::xx_rcut
  INTEGER                                            ::xx_atom_num
END TYPE md_info

TYPE classical_type
  INTEGER                             :: cp2k_env_flag
  INTEGER                             :: natoms 
  REAL(KIND=8)                        :: beta
  CHARACTER(LEN=90)                   :: cp2k_input_file,cp2k_output_file, cp2k_structure_file, cp2k_conn_file
  LOGICAL                                        ::atom_present
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)        ::i_position
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:)          ::ff_atom_charges
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:)          ::lang_atom_coeff
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)          ::force
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)        ::i_velocity
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:)          ::mass
  CHARACTER(LEN=2),ALLOCATABLE,DIMENSION(:)      ::atom_kinds
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)        ::x,p
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)      ::traj
  REAL(KIND=8)                                   ::kinetic, potential
END TYPE classical_type

TYPE parallel_env
  INTEGER                            ::my_rank
  INTEGER                            ::size_rank
END TYPE

END MODULE variable_types
