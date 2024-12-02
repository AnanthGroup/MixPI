MODULE pimd_mod
use variable_types
implicit none


TYPE(path_integral_env)                ::path_integral
TYPE(classical_type)                   ::classical
TYPE(md_info)                           ::md

contains



! ------------------------------------------------------------- !
! This the main MD loop. It will initialize the bead velocities !
! and then it will run nsteps of the MD intergrator with either !
! the results from the MC initialization or the input file. --- !
! The external potential can either be a model system (current- !
! modes are :                                                   !
! MODEL 1: External HO with omega_ext and mass_ext centered at- !
! 0.0.
! Or calculcated from potential generated via the positions of- !
! classical atoms. This is calculcated using CP2K ------------- !
! ------------------------------------------------------------- !
SUBROUTINE run_md
implicit none

logical         ::model_system = .FALSE.
integer         ::irp, istep, iprint, ibead, iclassical
integer         ::i,j

!!initalize velocities of the beads
!CALL initialize_bead_velocities

if (md%obtain_force == 'model') model_system = .TRUE.
! ------------------------------------------------------------- !
! ALLOCATE THE TRAJECTORY ARRAYS THAT SAVE THE POSITIONS AND -- !
! VELOCITIES AT THE FREQ_PRINT TIME STEPS. -------------------- !
! ALLOCATE THE FORCE, MOMENTA, AND POSITIONS ARRAYS THAT WILL - !
! BE USED AND OVERWRITTEN AT EACH TIME STEP DURING THE MD ----- !
! ROUTINE. ---------------------------------------------------- !
! ------------------------------------------------------------- !
! CALCULATE THE INITAL MOMENTUM FROM THE INITIAL VELOCITIES --- !
! ------------------------------------------------------------- !
! IF USING A NORMAL PROPAGATION, ALLOCATE THE NORMAL MODE ----- !
! POSITION AND MOMENTA ARRAYS. -------------------------------- !
! ------------------------------------------------------------- !
do irp = 1, path_integral%num_rp
  ALLOCATE(path_integral%rp(irp)%traj(0:md%nsteps/md%freq_print,&
    path_integral%rp(irp)%num_beads,2*md%space_dim))
  ALLOCATE(path_integral%com_pos(md%space_dim),path_integral%com_vel(md%space_dim))
  ALLOCATE(path_integral%rp(irp)%force(path_integral%rp(irp)%num_beads,md%space_dim))
  ALLOCATE(path_integral%rp(irp)%p(path_integral%rp(irp)%num_beads,md%space_dim))
  ALLOCATE(path_integral%rp(irp)%x(path_integral%rp(irp)%num_beads,md%space_dim))


  path_integral%rp(irp)%traj = 0.0d0
  path_integral%rp(irp)%force = 0.0d0
  path_integral%rp(irp)%p = 0.0d0
  path_integral%rp(irp)%x = 0.0d0

  path_integral%rp(irp)%p = path_integral%rp(irp)%kinetic_mass &
    *path_integral%rp(irp)%initialized_velocities
  path_integral%rp(irp)%x = path_integral%rp(irp)%initialized_positions

  do ibead = 1, path_integral%rp(irp)%num_beads
    path_integral%rp(irp)%traj(0,ibead,1:md%space_dim)=path_integral%rp(irp)%x(ibead,:)
    path_integral%rp(irp)%traj(0,ibead,md%space_dim+1:2*md%space_dim)&
       =path_integral%rp(irp)%p(ibead,1:md%space_dim)/path_integral%rp(irp)%kinetic_mass
  enddo

  if (trim(path_integral%propagation) == 'normal_fftw' .OR. trim(path_integral%propagation) == 'normal_manual') then
    ALLOCATE(path_integral%rp(irp)%norm_p(path_integral%rp(irp)%num_beads,md%space_dim))
    ALLOCATE(path_integral%rp(irp)%norm_x(path_integral%rp(irp)%num_beads,md%space_dim))
    path_integral%rp(irp)%norm_p = 0.0d0
    path_integral%rp(irp)%norm_x = 0.0d0
  endif

  ALLOCATE(path_integral%potential(0:md%nsteps),path_integral%kinetic(0:md%nsteps))
  ALLOCATE(path_integral%interaction_pot(0:md%nsteps), path_integral%inter_rp_pot(0:md%nsteps))
  path_integral%potential = 0.0d0
  path_integral%kinetic = 0.0d0
  path_integral%interaction_pot = 0.0d0
  path_integral%inter_rp_pot = 0.0d0

enddo

CALL update_force_env
CALL calc_total_energy(0)
iprint = 0
open(301,file='force_classical_print.dat')
open(302,file='force_pimd_print.dat')
open(303,file='momentum_classical_print.dat')
open(304,file='momentum_xcoord_print.dat')
open(305,file='momentum_normcoord_print.dat')
open(306,file='position_classical_print.dat')
open(307,file='position_xcoord_print.dat')
open(308,file='position_normcoord_print.dat')
do istep = 1, md%nsteps
! ------------------------------------------------- !
! EDIT NEEDED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !
! THE MD_RUN FOR A SPECIFIED EXTERNAL POTENTIAL --- !
! USES THE OLD MECHANISMS. NEEDS TO BE MODIFIED TO  !
! MATCH THE NEW CLASSICAL CP2K SUBROUTINES. ------- !
! ------------------------------------------------- !
  if(model_system) THEN
     if(istep == 1) write(*,*)' evaluate the MD steps using model potential ', md%model_no
     do irp = 1, path_integral%num_rp
       path_integral%rp(irp)%force = 0.0d0
     enddo
! HERE WE NEED TO PUT THE ANDERSON THERMOSTAT
     CALL vv_update_pimd
! SECOND PART OF ANDERSON THERMOSTAT
     if(mod(istep,md%freq_print)==0) THEN
       iprint = iprint+1
       do irp = 1, path_integral%num_rp
         do ibead = 1, path_integral%rp(irp)%num_beads
           path_integral%rp(irp)%traj(iprint,ibead,1:md%space_dim)=path_integral%rp(irp)%x(ibead,:)
           path_integral%rp(irp)%traj(iprint,ibead,md%space_dim+1:2*md%space_dim)&
             =path_integral%rp(irp)%p(ibead,1:md%space_dim)
         enddo
       enddo
     endif

! ------------------------------------------------- !
! THE MOLECULAR DYNAMICS LOOP FOR SIMULATIONS WITH  !
! PATH-INTEGRALS AND CLASSICAL ATOMS WHERE THE ---- !
! CLASSICAL FORCES ARE CALCULATED USING THE INTERNAL!
! SUBROUTINES IN THE LINKED CP2K LIBRARIES. ------- !
! ------------------------------------------------- !
  else if (trim(md%obtain_force) == 'cp2k') then
    if (istep == 1) write(*,*)' evalute the MD steps using external potential from CP2K'
    ! ------------------------------------------------ !
    ! BEGIN SUBROUTINE BY ZEROING OUT ALL FORCES FOR - !
    ! BOTH PATH-INTEGRAL AND CLASSICAL DOF ----------- !
    ! ------------------------------------------------ !
    do irp = 1, path_integral%num_rp
      path_integral%rp(irp)%force = 0.0d0
    enddo
    classical%force = 0.d0

    ! ------------------------------------------------ !
    ! EDIT NEEDED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !
    ! THE PILE_MIXEDTIME SUBROUTINE WILL NEED TO BE -- !
    ! CALLED HERE WHEN DOING NVT SIMULATIONS. -------- !
     !CALL pile_mixedtime_1    !call the initial part of the pile thermostat
    ! ------------------------------------------------ !

     ! ----------------------------------------------- !
     ! VELOCITY-VERLET STEP WITH THE CLASSICAL ATOMS - !
     ! FOR A TIMESTEP OF DT/2 ------------------------ !
     ! ----------------------------------------------- !
     write(*,*)' vv_update_classical-1 BAJ error'
     CALL vv_update_classical(md%dt/2.d0,istep)

     if (trim(path_integral%propagation) == 'normal_fftw') then
       write(*,*)' vv_update_pi BAJ error'
       CALL vv_update_pi_fftw(md%dt)
     else if (trim(path_integral%propagation) == 'normal_manual') then
       write(*,*)' vv_update_pi manual BAJ error'
       CALL vv_update_pi_normal_manual(md%dt,istep)
     endif

     write(*,*)' vv_update_pi done BAJ error'
     do irp = 1, path_integral%num_rp            !technically causes a repeat of this calculation
       path_integral%rp(irp)%force = 0.0d0       !either need to write a second vv_update_classical
     enddo                                       !or figure out another way.
     classical%force = 0.0d0

     write(*,*)' vv_update_classical-2 BAJ error'
     CALL vv_update_classical(md%dt/2.d0,istep)

     ! ----------------------------------------------- !
     ! AFTER FINAL CLASSICAL UDPATE FOR EACHS STEP, WE !
     ! NEED TO UPDATE THE FORCE_ENV INSIDE CP2K SO --- !
     ! THAT THE FORCE CALCULATIONS (ALONG WITH ENERGY, !
     ! ETC) WILL ALL BE UPDATED CORRECTLY FOR THE NEXT !
     ! STEP ------------------------------------------ !
     ! EDIT NEEDED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !
     ! ----------------------------------------------- !
     call update_force_env


     ! ----------------------------------------------- !
     ! CALCULCATE THE CONSERVED QUANTITIES (ENERGY FOR !
     ! NVE SIMULATIONS) AND PRINT TO FILE AT EACH STEP !
     ! ----------------------------------------------- !
     CALL calc_total_energy(istep)

     ! ----------------------------------------------- !
     ! BASED UPONNTHE INPUT FREQ_PRINT, WE WILL SAVE- !
     ! THE CLASSICAL AND RP POSITIONS IN A TRAJECTORY- !
     ! VARIABLE. 
     ! WE WILL PROBABLY EDIT THIS EVENTUALLY TO JUST-- !
     ! PRINT THE VALUES (INSTEAD OF SAVING ALL OF ---- !
     ! THESE POSITIONS. ------------------------------ !
     ! ----------------------------------------------- !
     if(mod(istep,md%freq_print)==0) THEN
       iprint = iprint+1
       do irp = 1, path_integral%num_rp
         do ibead = 1, path_integral%rp(irp)%num_beads
           path_integral%rp(irp)%traj(iprint,ibead,1:md%space_dim)=path_integral%rp(irp)%x(ibead,:)
           path_integral%rp(irp)%traj(iprint,ibead,md%space_dim+1:2*md%space_dim)&
             =path_integral%rp(irp)%p(ibead,1:md%space_dim)/path_integral%rp(irp)%kinetic_mass
         enddo
       enddo
       do iclassical = 1, classical%natoms
         classical%traj(iprint,iclassical,1:md%space_dim) = classical%x(iclassical,:)
         classical%traj(iprint,iclassical,md%space_dim+1:2*md%space_dim) = classical%p(iclassical,:)/classical%mass(iclassical)
       enddo
     endif

  endif

enddo
close(301)
close(302)
close(303)

!writes the trajectory file
do irp = 1, path_integral%num_rp
  open(201,file='trajectory_pimd_'//char(irp+48)//'.xyz')
  write(*,*)' print the md traj file '
  do i = 0, md%nsteps/md%freq_print
    write(201,*) path_integral%rp(irp)%num_beads
    write(201,*)' i = ', i
    do j = 1, path_integral%rp(irp)%num_beads
      write(201,'(a2, 6f15.8)') 'H', path_integral%rp(irp)%traj(i,j,1:3)/AtoBohr, &
        path_integral%rp(irp)%traj(i,j,4:6)/(AtoBohr)
    enddo
  enddo
  close(201) 
enddo

open(202,file='trajectory_classical.xyz')
  do i = 1, md%nsteps/md%freq_print
    write(202,*) classical%natoms
    write(202,*)' i = ', i
    do j = 1, classical%natoms
      write(202,'(a2, 6f15.8)') trim(classical%atom_kinds(j)), classical%traj(i,j,1:3)/AtoBohr, &
        classical%traj(i,j,4:6)/(AtoBohr/fstota)
    enddo
  enddo
close(202)

open(203,file='trajectory_full.xyz')
  do i = 0, md%nsteps/md%freq_print
    write(203,*) path_integral%rp(1)%num_beads + classical%natoms
    write(203,*)' i = ', i
    do j = 1, classical%natoms
      write(203,'(a2, 6f15.8)') trim(classical%atom_kinds(j)), classical%traj(i,j,1:3)/AtoBohr, &
        classical%traj(i,j,4:6)/(AtoBohr/fstota)
    enddo
    do j = 1, path_integral%rp(1)%num_beads
      write(203,'(a2, 6f15.8)') 'H', path_integral%rp(1)%traj(i,j,1:3)/AtoBohr, &
        path_integral%rp(1)%traj(i,j,4:6)/(AtoBohr/fstota)
    enddo
  enddo
  close(203)

open(204,file='trajectory_com.xyz')
do i = 0, md%nsteps/md%freq_print
  write(204,*) path_integral%num_rp + classical%natoms
  write(204,*)' i = ', i
  do j = 1, classical%natoms
    write(204,'(a2, 9f15.8)') trim(classical%atom_kinds(j)), classical%traj(i,j,1:3)/(AtoBohr), &
    classical%traj(i,j,4:6), classical%traj(i,j,4:6)/(AtoBohr/fstota)
  enddo
  do irp = 1, path_integral%num_rp
    path_integral%com_pos = 0.0d0
    path_integral%com_vel = 0.0d0
    do j = 1, path_integral%rp(irp)%num_beads
      path_integral%com_pos(1:md%space_dim) = path_integral%com_pos(1:md%space_dim) + &
        path_integral%rp(irp)%traj(i,j,1:md%space_dim)
      path_integral%com_vel(1:md%space_dim) = path_integral%com_vel(1:md%space_dim) + &
        path_integral%rp(irp)%traj(i,j,md%space_dim+1:2*md%space_dim)
    enddo
    path_integral%com_pos(:) = path_integral%com_pos(:)/dble(path_integral%rp(irp)%num_beads)
    path_integral%com_vel(:) = path_integral%com_vel(:)/dble(path_integral%rp(irp)%num_beads)
    write(204,'(a2, 9f15.8)') 'H', path_integral%com_pos(1:md%space_dim)/AtoBohr, &
        path_integral%com_vel(1:md%space_dim)/(AtoBohr), path_integral%com_vel(1:md%space_dim)/(AtoBohr/fstota)
  enddo
enddo

open(203,file='energy.dat')
do i = 0, md%nsteps
  write(203,'(i5,f16.2,7f20.10)')i, i*md%dt,classical%kinetic(i),classical%potential(i), path_integral%kinetic(i), &
    path_integral%potential(i), path_integral%interaction_pot(i), path_integral%inter_rp_pot(i), md%total_energy(i)
enddo
close(203)

END SUBROUTINE run_md
! ----------------------------------------------------!




SUBROUTINE update_force_env
use libcp2k
implicit none

REAL(KIND=8),DIMENSION(md%space_dim*classical%natoms)     ::positions
REAL(KIND=8),DIMENSION(md%space_dim*classical%natoms)     ::velocities
integer                              ::i


do i = 1, classical%natoms
  positions( (i-1)*md%space_dim+1: i*md%space_dim ) = classical%x(i,1:md%space_dim) 
  velocities( (i-1)*md%space_dim+1: i*md%space_dim ) = classical%p(i,1:md%space_dim)/classical%mass(i)
enddo

CALL cp2k_set_positions(classical%cp2k_env_flag,positions,md%space_dim*classical%natoms) 
CALL cp2k_set_velocities(classical%cp2k_env_flag,velocities,md%space_dim*classical%natoms)

END SUBROUTINE update_force_env


SUBROUTINE calc_total_energy(istep)
implicit none

real(kind=8)       ::classical_pot, classical_kin
real(kind=8)       ::pimd_pot, pimd_kin
real(kind=8)       ::interaction_pot, inter_rp_pot
integer            ::istep

CALL calc_classical_energy(classical_pot, classical_kin)

CALL calc_pimd_energy(pimd_pot, pimd_kin)

! ------------------------------------------ !
! EDIT NEEDED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !
! BECAUSE CALC_INTERACTION_POTENTIAL ASSUMES !
! ALL RP HAVE SAME NUMBER OF BEADS. -------- !
! ------------------------------------------ !
CALL calc_interaction_potential(interaction_pot, inter_rp_pot)

classical%potential(istep) = classical_pot
classical%kinetic(istep) = classical_kin
path_integral%potential(istep) = pimd_pot
path_integral%kinetic(istep) = pimd_kin
path_integral%interaction_pot(istep) = interaction_pot
path_integral%inter_rp_pot(istep) = inter_rp_pot
md%total_energy(istep) = classical_pot + classical_kin + pimd_pot + pimd_kin + interaction_pot + inter_rp_pot

END SUBROUTINE calc_total_energy


SUBROUTINE calc_pimd_energy(pot_energy,kin_energy)
implicit none


real(kind=8)          ::pot_energy, kin_energy, omega
integer               ::irp,i,iplus,nbead

pot_energy = 0.0d0
do irp = 1, path_integral%num_rp
  omega = 1.0d0/path_integral%rp(irp)%betaN
  nbead = path_integral%rp(irp)%num_beads
  do i = 1, nbead
    iplus = i+1
    if (i==nbead) iplus = 1
    pot_energy = pot_energy + 0.5d0*path_integral%rp(irp)%harm_mass*omega**2/dble(path_integral%rp(irp)%num_beads) * &
      dot_product(path_integral%rp(irp)%x(i,:)-path_integral%rp(irp)%x(iplus,:),&
      path_integral%rp(irp)%x(i,:)-path_integral%rp(irp)%x(iplus,:))
    write(*,*)' prefactor is ', 0.5d0*path_integral%rp(irp)%harm_mass*omega**2/dble(path_integral%rp(irp)%num_beads) 
    write(*,*)' dot product is ', dot_product(path_integral%rp(irp)%x(i,:)-path_integral%rp(irp)%x(iplus,:),&
      path_integral%rp(irp)%x(i,:)-path_integral%rp(irp)%x(iplus,:))
    kin_energy = kin_energy + 0.5d0*dot_product(path_integral%rp(irp)%p(i,:),path_integral%rp(irp)%p(i,:))/&
      path_integral%rp(irp)%kinetic_mass
  enddo
enddo

END SUBROUTINE calc_pimd_energy

SUBROUTINE calc_classical_energy(pot_energy, kin_energy)
use libcp2k
implicit none

real(kind=8)          ::pot_energy, kin_energy
integer               ::i

CALL cp2k_calc_energy_force(classical%cp2k_env_flag)
CALL cp2k_get_potential_energy(classical%cp2k_env_flag, pot_energy)
! --------------------------------------------- !
! WOULD EVENTUALLY LIKE TO HAVE A CALL -------- !
! TO INTERNAL CP2K FUNCTION CP2K_GET_ENERGY --- !
! BUT UNTIL THEN WE WILL MANUAL CALCULATE THE - !
! KINETIC ENERGY. ----------------------------- !
! --------------------------------------------- !
kin_energy = 0.0d0
do i = 1, classical%natoms
  kin_energy = kin_energy + 0.5d0 * dot_product(classical%p(i,:),classical%p(i,:))/classical%mass(i)
enddo



END SUBROUTINE calc_classical_energy



! --------------------------------------------- !
! The velocity-verlet integrator for the path - !
! integral degrees of freedom. This subroutine  !
! differs from the previous vv_pimd scheme by - !
! using the normal mode propagation. ---------- !
! --------------------------------------------- !
SUBROUTINE vv_update_pi_fftw(dt)
implicit none

real(kind=8)          :: dt
integer               :: irp

write(*,*)' update_pi_p BAJ error'
CALL update_pi_p(dt)

write(*,*)' fftw_cart_to_normal BAJ error'
CALL fftw_cart_to_normal
write(*,*)' update_normal_coord BAJ error'
CALL update_normal_coord(dt)
write(*,*)' fftw_normal_to_cart BAJ error'
CALL fftw_normal_to_cart

do irp = 1, path_integral%num_rp
  path_integral%rp(irp)%force = 0.0d0
enddo
write(*,*)' calc_force_mix BAJ error'
CALL calc_force_mix(.TRUE., .FALSE.)   !update_rp, update_classical
write(*,*)' calc_force_pi BAJ error'
CALL calc_force_pi

CALL update_pi_p(dt)

END SUBROUTINE vv_update_pi_fftw


! ----------------------------------------------- !
! The velocity-verlet algorithm for the path ---- !
! integral degrees of freedome. This subroutine - !
! manual calculates the  matrix for the --------- !
! transformation from the cartesian to normal --- !
! coordinates. ---------------------------------- !
! ----------------------------------------------- !
SUBROUTINE vv_update_pi_normal_manual(dt,icount)
implicit none

real(kind=8)          :: dt
integer               :: irp,icount,i

CALL update_pi_p(dt)
do irp = 1, path_integral%num_rp
  do i = 1, path_integral%rp(irp)%num_beads
    write(304,'(3i4,10f15.8)') icount,irp,i,path_integral%rp(irp)%p(i,:)
  enddo
enddo

write(*,*)' calculcate transform matrix '
CALL manual_normal_rotation(dt)

write(*,*)' transform cart to normal '
CALL manual_cart_to_normal(icount)
write(307,*)' transform cart to normal '
write(305,*)' transform cart to normal '
write(308,*)' transform cart to normal '
do irp = 1, path_integral%num_rp
  do i = 1, path_integral%rp(irp)%num_beads
    write(307,'(3i4,10f15.8)') icount,irp,i,path_integral%rp(irp)%norm_p(i,:)
  enddo
enddo
do irp = 1, path_integral%num_rp
  do i = 1, path_integral%rp(irp)%num_beads
    write(305,'(3i4,10f15.8)') icount,irp,i,path_integral%rp(irp)%x(i,:)
  enddo
enddo
do irp = 1, path_integral%num_rp
  do i = 1, path_integral%rp(irp)%num_beads
    write(308,'(3i4,10f15.8)') icount,irp,i,path_integral%rp(irp)%norm_x(i,:)
  enddo
enddo

write(*,*)' update'
CALL update_normal_coord(dt)
write(307,*)' update coord'
write(308,*)' update coord'
do irp = 1, path_integral%num_rp
  do i = 1, path_integral%rp(irp)%num_beads
    write(307,'(3i4,10f15.8)') icount,irp,i,path_integral%rp(irp)%norm_p(i,:)
  enddo
enddo
do irp = 1, path_integral%num_rp
  do i = 1, path_integral%rp(irp)%num_beads
    write(308,'(3i4,10f15.8)') icount,irp,i,path_integral%rp(irp)%norm_x(i,:)
  enddo
enddo

write(*,*)' transform normal to cart '
CALL manual_normal_to_cart(icount)
write(304,*)' transform normal to cart '
write(305,*)' transform normal to cart '
do irp = 1, path_integral%num_rp
  do i = 1, path_integral%rp(irp)%num_beads
    write(305,'(3i4,10f15.8)') icount,irp,i,path_integral%rp(irp)%x(i,:)
  enddo
enddo
do irp = 1, path_integral%num_rp
  do i = 1, path_integral%rp(irp)%num_beads
    write(304,'(3i4,10f15.8)') icount,irp,i,path_integral%rp(irp)%p(i,:)
  enddo
enddo

do irp = 1, path_integral%num_rp
  path_integral%rp(irp)%force = 0.0d0
enddo
write(*,*)' calc_force_mix BAJ error'
CALL calc_force_mix(.TRUE., .FALSE.)   !update_rp, update_classical
write(*,*)' calc_force_pi BAJ error'
CALL calc_force_pi
do irp = 1, path_integral%num_rp
  do i = 1, path_integral%rp(irp)%num_beads
    write(302,'(3i4,10f15.8)') icount,irp,i,path_integral%rp(irp)%force(i,:)
  enddo
enddo

CALL update_pi_p(dt)
do irp = 1, path_integral%num_rp
  do i = 1, path_integral%rp(irp)%num_beads
    write(304,'(3i4,10f15.8)') icount,irp,i,path_integral%rp(irp)%p(i,:)
  enddo
enddo

END SUBROUTINE vv_update_pi_normal_manual


SUBROUTINE manual_cart_to_normal(istep)
implicit none
integer            ::irp,ispace
integer            ::j,k, i
integer            ::istep

if (istep==1) then
  open(777,file='rot_matrix_check.dat')
  do i = 1, path_integral%rp(1)%num_beads
    write(777,'(100f12.4)') path_integral%rp(1)%rot_matrix(i,:)
  enddo
  close(777)
endif

do irp = 1, path_integral%num_rp
  path_integral%rp(irp)%norm_p = 0.0d0
  path_integral%rp(irp)%norm_x = 0.0d0
  do ispace = 1, md%space_dim
    do j = 1, path_integral%rp(irp)%num_beads
      do k = 1, path_integral%rp(irp)%num_beads
        path_integral%rp(irp)%norm_p(j,ispace) = path_integral%rp(irp)%norm_p(j,ispace) + path_integral%rp(irp)%p(k,ispace) &
          * path_integral%rp(irp)%rot_matrix(k,j)
        path_integral%rp(irp)%norm_x(j,ispace) = path_integral%rp(irp)%norm_x(j,ispace) + path_integral%rp(irp)%x(k,ispace) &
          * path_integral%rp(irp)%rot_matrix(k,j)
      enddo
    enddo
  enddo
enddo

END SUBROUTINE manual_cart_to_normal

SUBROUTINE manual_normal_to_cart(istep)
implicit none
integer              ::irp,ispace
integer              ::j,k, i
integer              ::istep
if (istep==1) then
  open(777,file='rot_matrix_check_2.dat')
  do i = 1, path_integral%rp(1)%num_beads
    write(777,'(100f12.4)') path_integral%rp(1)%rot_matrix(i,:)
  enddo
  close(777)
endif

do irp = 1, path_integral%num_rp
  path_integral%rp(irp)%p = 0.0d0
  path_integral%rp(irp)%x = 0.0d0
  do ispace = 1, md%space_dim
    do j = 1, path_integral%rp(irp)%num_beads
      do k = 1, path_integral%rp(irp)%num_beads
        path_integral%rp(irp)%p(k,ispace) = path_integral%rp(irp)%p(k,ispace) + path_integral%rp(irp)%rot_matrix(k,j) & 
          * path_integral%rp(irp)%norm_p(j,ispace)
        path_integral%rp(irp)%x(k,ispace) = path_integral%rp(irp)%x(k,ispace) + path_integral%rp(irp)%rot_matrix(k,j) & 
          * path_integral%rp(irp)%norm_x(j,ispace)
      enddo
    enddo
  enddo
deallocate(path_integral%rp(irp)%rot_matrix)
enddo

END SUBROUTINE manual_normal_to_cart

SUBROUTINE manual_normal_rotation(dt)
implicit none

real(kind=8)     ::dt
integer          ::irp, ibead, jbead, nbeads
integer          ::i

real(kind=8),allocatable,dimension(:,:)    ::rot_matrix

do irp = 1, path_integral%num_rp
  
  nbeads = path_integral%rp(irp)%num_beads
  allocate(path_integral%rp(irp)%rot_matrix(nbeads,nbeads),rot_matrix(nbeads,nbeads))
  do ibead = 1, nbeads
    if ( (ibead-1) == 0 ) then
      rot_matrix(:,ibead) = 1.0d0/dsqrt(dble(nbeads))
    else if ( (ibead-1) .le. nbeads/2-1) then
      do jbead = 1, nbeads
        rot_matrix(jbead,ibead) = dsqrt(2.0d0/dble(nbeads)) * dcos(2.0*pi*(ibead-1)*(jbead-1)/dble(nbeads))
      enddo
    else if ( (ibead-1) == nbeads/2) then
      do jbead = 1, nbeads
        rot_matrix(jbead,ibead) = dsqrt(1.0d0/dble(nbeads)) * (-1.0d0)**(jbead-1)
      enddo
    else if ( (ibead-1) .le. nbeads-1) then
      do jbead = 1, nbeads
        rot_matrix(jbead,ibead) = dsqrt(2.0d0/dble(nbeads))*dsin(2.0*pi*(ibead-1)*(jbead-1)/dble(nbeads))
      enddo
    endif
  enddo

  open(888,file='rotation_matrix.dat')
  do i = 1, nbeads
    write(888,'(100f10.2)') rot_matrix(i,:)
  enddo
  close(888)

  path_integral%rp(irp)%rot_matrix = rot_matrix
  deallocate(rot_matrix)
  !do i_dim = 1, md%space_dim
  !  path_integral%rp(irp)%p(:,i_dim) = path_integral%rp(irp)%norm_p(:,i_dim)*rot_matrix(
  !MATRIX MUTLIPLICATION TO PROPOGATE THE NORMAL MODE COORDINATES.I
enddo

END SUBROUTINE manual_normal_rotation

! --------------------------------------------- !
! Subroutine that uses the FFT3 to convert the  !
! cartesian representation into the normal mode !
! free ring polymer representation. ----------- !
! --------------------------------------------- !
SUBROUTINE fftw_cart_to_normal
use, intrinsic :: iso_c_binding
implicit none


integer           ::i
integer           ::irp, i_dim, nbeads
type(c_ptr)         ::plan1, plan2
real(kind=8)              ::factor
!real(kind=8),allocatable,dimension(:)     ::copy1a,copy1b
real(c_double),allocatable,dimension(:)     :: copy1a, copy1b

!point to fftw3 library
include '/share/apps/intel/2015u2/mkl/include/fftw/fftw3.f'
!#include "fftw3.f"

do irp = 1, path_integral%num_rp
  nbeads = path_integral%rp(irp)%num_beads
  write(*,*)' norm_p rotate BAJ error'
  do i_dim = 1, md%space_dim
    write(*,*)' allocate BAJ error'
    allocate(copy1a(nbeads), copy1b(nbeads))
    copy1a = 0.d0
    copy1b = 0.d0
    write(*,*)' set 1a BAJ error'
    copy1a = path_integral%rp(irp)%p(:,i_dim)
    write(*,*)' dfftw_plan_dft_r2r_1d BAJ error'
    call dfftw_plan_r2r_1d(plan1, nbeads, copy1a, copy1b, FFTW_FORWARD, FFTW_MEASURE)
   ! plan1 = fftw_plan_r2r_1d(nbeads, copy1a, copy1b, FFTW_FORWARD, FFTW_MEASURE)
    write(*,*)' factor BAJ error'
    write(*,*)' plan is ', plan1, plan2
    factor = dsqrt(1.0/dble(nbeads))

    !write(*,*)' set 1a BAJ error'
    !copy1a = path_integral%rp(irp)%p(:,i_dim)
    write(*,*)' dfftw_execute BAJ error'
    CALL dfftw_execute_r2r(plan1, copy1a, copy1b)   !DIFFERENT THAN JESS'S CODE??
    write(*,*)' set norm_p BAJ error'
    path_integral%rp(irp)%norm_p(:,i_dim) = factor*copy1b
    ! ----------------------------------------- !
    ! EDIT NEEDED !!!!!!!!!!!!!!!!!!!!!!!!!!!!! !
    ! WON'T LET ME DEALLOCATE/DESTORY PLAN? --- !
    write(*,*)' destroy BAJ error'                  
    CALL dfftw_destroy_plan(plan1)
    ! ----------------------------------------- !
    write(*,*)' deallocate BAJ error'
    deallocate(copy1a)
    deallocate(copy1b)
  enddo
  write(*,*)' norm_x rotate BAJ error'
  do i_dim = 1, md%space_dim
    allocate(copy1a(nbeads), copy1b(nbeads))
    call dfftw_plan_r2r_1d(plan1, nbeads, copy1a, copy1b, FFTW_FORWARD, FFTW_MEASURE)
    factor = dsqrt(1.0/dble(nbeads))

    copy1a = path_integral%rp(irp)%x(:,i_dim)
    CALL dfftw_execute_r2r(plan1, copy1a, copy1b)
    path_integral%rp(irp)%norm_x(:,i_dim) = factor*copy1b
    deallocate(copy1a, copy1b)
    CALL dfftw_destroy_plan(plan1)
  enddo
enddo


open(401,file='check_fftw_cart_to_normal.out')
do i=1,path_integral%rp(1)%num_beads
  write(401,'(3f10.2)') path_integral%rp(1)%norm_x(i,1:3)
enddo
close(401)


END SUBROUTINE fftw_cart_to_normal


! --------------------------------------------- !
! Subroutine that uses the FFT3 to convert the  !
! normal representation into the cartesian ---- !
! free ring polymer representation. ----------- !
! --------------------------------------------- !
SUBROUTINE fftw_normal_to_cart
implicit none

integer           ::i
integer           ::irp, i_dim, nbeads
real(kind=8)              ::plan1
real(kind=8)              ::factor
real(kind=8),allocatable,dimension(:)     ::copy1a,copy1b

!point to fftw3 library
include '/share/apps/intel/2015u1/mkl/include/fftw/fftw3.f'

do irp = 1, path_integral%num_rp
  nbeads = path_integral%rp(irp)%num_beads
  do i_dim = 1, md%space_dim
    allocate(copy1a(nbeads), copy1b(nbeads))
    call dfftw_plan_dft_1d(plan1, nbeads, copy1a, copy1b, FFTW_FORWARD, FFTW_MEASURE)
    factor = dsqrt(1.0/dble(nbeads))

    copy1a = path_integral%rp(irp)%norm_p(:,i_dim)
    CALL dfftw_execute_dft(plan1, copy1a, copy1b)   !DIFFERENT THAN JESS'S CODE??
    path_integral%rp(irp)%p(:,i_dim) = factor*copy1b
    deallocate(copy1a, copy1b)
    CALL dfftw_destroy_plan(plan1)
  enddo
  do i_dim = 1, md%space_dim
    allocate(copy1a(nbeads), copy1b(nbeads))
    call dfftw_plan_dft_1d(plan1, nbeads, copy1a, copy1b, FFTW_FORWARD, FFTW_MEASURE)
    factor = dsqrt(1.0/dble(nbeads))

    copy1a = path_integral%rp(irp)%norm_x(:,i_dim)
    CALL dfftw_execute_dft(plan1, copy1a, copy1b)
    path_integral%rp(irp)%x(:,i_dim) = factor*copy1b
    deallocate(copy1a, copy1b)
    CALL dfftw_destroy_plan(plan1)
  enddo
enddo


open(401,file='check_fftw_normal_to_cart.out')
do i=1,path_integral%rp(1)%num_beads
  write(401,'(3f10.2)') path_integral%rp(1)%x(i,1:3)
enddo
close(401)


END SUBROUTINE fftw_normal_to_cart


! --------------------------------------------- !
! Subroutine that updates the momentum of the - !
! path integral beads using the --------------- !
! verlet scheme.                                !
! p_j(t + delt/2)  <- p_j(t) + delt/2 * F_j     !
! --------------------------------------------- !
SUBROUTINE update_pi_p(dt)
implicit none
integer                          ::irp
integer                          ::ibead
real(kind=8)                     ::dt

do irp = 1, path_integral%num_rp
  path_integral%rp(irp)%p(:,:) = path_integral%rp(irp)%p(:,:) + dt * path_integral%rp(irp)%force(:,:)/2.0d0 
enddo

END SUBROUTINE update_pi_p
! -------------------------------------------- !




! ---------------------------------------------- !
! The velocity-verlet interrator for the ------- !
! classical atoms with a time-step of dt ------- !
! The following updates occur:  ---------------- !
! p(t + dt/2) <- p(t) + dt/2 * F(t) ------------ !
! x(t+dt) <- x(t) + dt * p(t+dt/2)/mass -------- !
! p(t+dt) <- p(t + dt/2) + dt/2 * F(t+dt/2) ---- !
! ---------------------------------------------- !
SUBROUTINE vv_update_classical(dt,icount)
implicit none

real(kind=8)            ::dt
integer                 ::irp, icount,i

!classical%force = 0.0d0
!do irp = 1, path_integral%num_rp
!  path_integral%rp(irp)%force = 0.0d0
!enddo

CALL calc_force_classical
write(301,*) 'classical-classical forces '
do i = 1, classical%natoms
  write(301,'(2i4,10f15.8)') icount,i, classical%force(i,:)
enddo
CALL calc_force_mix(.TRUE.,.TRUE.)    !logical whether to update the RP forces
write(302,*)' interaction forces '
do irp = 1, path_integral%num_rp
  do i = 1, path_integral%rp(irp)%num_beads
    write(302,'(3i4,10f15.8)') icount,irp,i,path_integral%rp(irp)%force(i,:)
  enddo
enddo
CALL calc_force_pi

write(301,*) 'classical-classical and interaction forces '
do i = 1, classical%natoms
  write(301,'(2i4,10f15.8)') icount,i, classical%force(i,:)
enddo
write(302,*)' interaction and rp-rp forces '
do irp = 1, path_integral%num_rp
  do i = 1, path_integral%rp(irp)%num_beads
    write(302,'(3i4,10f15.8)') icount,irp,i,path_integral%rp(irp)%force(i,:)
  enddo
enddo

CALL update_classical_p(dt)
do i = 1, classical%natoms
  write(303,'(2i4,10f15.8)') icount,i, classical%p(i,:)
enddo
CALL update_classical_x(dt)
do i = 1, classical%natoms
  write(306,'(2i4,10f15.8)') icount,i, classical%x(i,:)
enddo

classical%force = 0.0d0
CALL update_force_env
CALL calc_force_classical
CALL calc_force_mix(.FALSE.,.TRUE.)

do i = 1, classical%natoms
  write(301,'(2i4,10f15.8)') icount,i, classical%force(i,:)
enddo
do irp = 1, path_integral%num_rp
  do i = 1, path_integral%rp(irp)%num_beads
    write(302,'(3i4,10f15.8)') icount,irp,i,path_integral%rp(irp)%force(i,:)
  enddo
enddo

CALL update_classical_p(dt)
do i = 1, classical%natoms
  write(303,'(2i4,10f15.8)') icount,i, classical%p(i,:)
enddo

END SUBROUTINE vv_update_classical
! ---------------------------------------------- !


! ---------------------------------------------- !
! Subroutine that updates the classical atom --- !
! positions by dt ------------------------------ !
! ---------------------------------------------- !
SUBROUTINE update_classical_x(dt)
implicit none
real(kind=8)       ::dt
integer            ::iatom

do iatom = 1, classical%natoms
  classical%x(iatom,:) = classical%x(iatom,:) + dt * classical%p(iatom,:)/( classical%mass(iatom) )
enddo

END SUBROUTINE update_classical_x

! ---------------------------------------------- !
! Subroutine that updates the classical atom --- !
! momentum by dt/2 ----------------------------- !
! ---------------------------------------------- !
SUBROUTINE update_classical_p(dt)
implicit none

real(kind=8)                ::dt

classical%p(:,:) = classical%p(:,:) + dt * classical%force(:,:)/2.0d0

END SUBROUTINE update_classical_p



! ---------------------------------------------- !
! Subroutine that evaluates the forces on the -- !
! classical atoms using the CP2K subroutines --- !
! ---------------------------------------------- !
SUBROUTINE calc_force_classical
use libcp2k
implicit none

real(kind=8),allocatable,dimension(:)    ::force
integer                                  ::i

ALLOCATE(force(classical%natoms*md%space_dim))
force = 0.0d0

! ------------------------------------------------ !
! CALL CP2K INTERNAL FUNCTION CP2K_GET_FORCES ---- !
! THAT CALCULATES THE FORCES BETWEEN CLASSICAL --- !
! ATOMS ARE RETURNS THEM IN AN ARRAY OF SIZE ----- !
! NATOMS*DIMENSION ------------------------------- !
! ------------------------------------------------ !
CALL cp2k_calc_energy_force(classical%cp2k_env_flag)
CALL cp2k_get_forces(classical%cp2k_env_flag, force, classical%natoms*md%space_dim)

! ------------------------------------------------ !
! CONVERTS THE ARRAY FORCES(1:NATOMS*DIMENSION) -- !
! INTO ARRAY CLASSICAL%FORCE(1:NATOMS,1:DIMENSION) !
! WHCIH IS JUST MORE INTUITIVE. ------------------ !
! IMPORTANT NOTE: IF CLASSICAL$FORCE ARRAY HAS NOT !
! BEEN ZERO'D BEFORE CALLING SUBROUTINE, IT WILL - !
! BE ADDED TO THE EXISTING FORCES THERE. --------- !
! ------------------------------------------------ !
do i = 1, classical%natoms
  classical%force(i,1:md%space_dim) = classical%force(i,1:md%space_dim) + force( (i-1)*md%space_dim+1:i*md%space_dim )
enddo

DEALLOCATE(force)

END SUBROUTINE calc_force_classical
! --------------------------------------------- !


! --------------------------------------------- !
! Subroutine that will calculate the forces --- !
! between the beads and the classical atoms. -- !
! --------------------------------------------- !
SUBROUTINE calc_force_mix(update_pi,update_classical)
implicit none

integer               ::irp
integer               ::iclassical, ibead, i
real(kind=8),allocatable,dimension(:,:)    ::force_c,force_pi
logical               ::update_pi,update_classical

!anything else? some sort of cut off?

ALLOCATE(force_c(classical%natoms,md%space_dim))
force_c = 0.0d0

write(*,*)' classical epsilon is ', classical%ff_atom_lj(:,1)
write(*,*)' classical sigma is ', classical%ff_atom_lj(:,2)
write(*,*)' RP epsilon is ', path_integral%rp(1)%lj_param(1)/(dble(path_integral%rp(1)%num_beads)**2)
write(*,*)' RP sigma is ', path_integral%rp(1)%lj_param(2)

do irp = 1, path_integral%num_rp
  ALLOCATE(force_pi(path_integral%rp(irp)%num_beads,md%space_dim))
  force_pi = 0.0d0

  do iclassical = 1, classical%natoms
    do ibead = 1, path_integral%rp(irp)%num_beads
      ! -------------------------------------------- !
      ! Coulombic forces between classical atom and  !
      ! RP beads. ---------------------------------- !
      ! -------------------------------------------- !

      ! ------------------------------------------- !
      ! EDIT NEEDED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! If atom and RP beads are oppositely charged !
      ! we need to include a "cutoff" value to ---- !
      ! prevent the particles from crashing into -- !
      ! each other. ------------------------------- !
      ! ------------------------------------------- !
      if(update_classical) then
        force_c(iclassical,1:md%space_dim) = force_c(iclassical,1:md%space_dim) + & 
          coulomb_force( classical%ff_atom_charges(iclassical),path_integral%rp(irp)%charge/dble(path_integral%rp(irp)%num_beads), &
          classical%x(iclassical,1:md%space_dim), path_integral%rp(irp)%x(ibead,1:md%space_dim) )
      endif
      if(update_pi) then
        force_pi(ibead,1:md%space_dim) = force_pi(ibead,1:md%space_dim) - &
          coulomb_force( classical%ff_atom_charges(iclassical),path_integral%rp(irp)%charge/dble(path_integral%rp(irp)%num_beads), &
          classical%x(iclassical,1:md%space_dim), path_integral%rp(irp)%x(ibead,1:md%space_dim) )
      endif



      ! ------------------------------------------- !
      ! Dispersion forces between classical atoms - !
      ! and RP beads. Currently, only have a LJ --- !
      ! force option. ----------------------------- !
      ! ------------------------------------------- !

      if(update_classical) then
        force_c(iclassical,1:md%space_dim) = force_c(iclassical,1:md%space_dim) + &
          lj_force( classical%ff_atom_lj(iclassical,1), classical%ff_atom_lj(iclassical,2), &
          path_integral%rp(irp)%lj_param(1)/(dble(path_integral%rp(irp)%num_beads)**2), path_integral%rp(irp)%lj_param(2), &
          classical%x(iclassical,1:md%space_dim), path_integral%rp(irp)%x(ibead,1:md%space_dim) )
      endif

      if (update_pi) then
        force_pi(ibead,1:md%space_dim) = force_pi(ibead,1:md%space_dim) - &
        lj_force( classical%ff_atom_lj(iclassical,1), classical%ff_atom_lj(iclassical,2), &
        path_integral%rp(irp)%lj_param(1)/(dble(path_integral%rp(irp)%num_beads)**2), path_integral%rp(irp)%lj_param(2), &
        classical%x(iclassical,1:md%space_dim), path_integral%rp(irp)%x(ibead,1:md%space_dim) )
      endif

    enddo
  enddo

  do i = 1, path_integral%rp(irp)%num_beads
    path_integral%rp(irp)%force(i,1:md%space_dim) = path_integral%rp(irp)%force(i,1:md%space_dim)&
      + force_pi(i,1:md%space_dim)
  enddo

  DEALLOCATE(force_pi)
enddo

do i = 1, classical%natoms
  classical%force(i,1:md%space_dim) = classical%force(i,1:md%space_dim) + force_c(i,1:md%space_dim)
enddo
DEALLOCATE(force_c)

END SUBROUTINE calc_force_mix


FUNCTION lj_force(eps1, sig1, eps2, sig2, r1, r2) RESULT(force)
implicit none

real(kind=8)     ::eps1,sig1,eps2,sig2
real(kind=8),dimension(md%space_dim)    ::r1,r2
real(kind=8),dimension(md%space_dim)    ::force
real(kind=8)     ::distance, eps, sig

distance = dsqrt( dot_product(r1-r2,r1-r2) )
write(*,*)' distance is ', distance
write(*,*)' eps is ', eps1, eps2
write(*,*)' sig is ', sig1, sig2
! ------------------------------------------------ !
! Use the Lorentz-Berthelot rules for combinations !
! ------------------------------------------------ !
sig = (sig1 + sig2)/2.0d0
eps = dsqrt(eps1*eps2)

force = 24.0d0*eps/distance * ( 2.0d0*(sig/distance)**12 - (sig/distance)**6 ) * (r1-r2)/distance

write(*,*)' force is ', force

END FUNCTION lj_force


FUNCTION lj_pot(eps1, sig1, eps2, sig2, r1, r2) RESULT(pot)
implicit none

real(kind=8)     ::eps1,sig1,eps2,sig2
real(kind=8),dimension(md%space_dim)    ::r1,r2
real(kind=8)     ::pot
real(kind=8)     ::distance, eps, sig

distance = dsqrt( dot_product(r1-r2,r1-r2) )
! ------------------------------------------------ !
! Use the Lorentz-Berthelot rules for combinations !
! ------------------------------------------------ !
sig = (sig1 + sig2)/2.0d0
eps = dsqrt(eps1*eps2)

pot = 4.0d0*eps * ( (sig/distance)**12 - (sig/distance)**6 )

END FUNCTION lj_pot



FUNCTION coulomb_force(q1, q2, r1, r2) RESULT(force)
implicit none

real(kind=8)       ::q1,q2
real(kind=8),dimension(md%space_dim)    ::r1,r2
real(kind=8),dimension(md%space_dim)    ::force
real(kind=8)       ::distance

distance = dsqrt( dot_product(r1-r2, r1-r2) )

force = ( q1*q2/distance) * (r1 - r2)/distance

END FUNCTION coulomb_force


FUNCTION coulomb_pot(q1, q2, r1, r2) RESULT(pot)
implicit none

real(kind=8)       ::q1,q2
real(kind=8),dimension(md%space_dim)    ::r1,r2
real(kind=8)       ::pot
real(kind=8)       ::distance

distance = dsqrt( dot_product(r1-r2, r1-r2) )

pot = ( q1*q2/distance)

END FUNCTION coulomb_pot



SUBROUTINE calc_force_pi
implicit none

if (trim(path_integral%propagation) == 'NONE') then
  CALL calc_force_pi_internal
endif

if (path_integral%num_rp > 1) then
  CALL calc_force_pi_pi
endif

END SUBROUTINE calc_force_pi



SUBROUTINE calc_force_pi_pi
implicit none
integer                 ::irp1,irp2
real(kind=8),dimension(md%space_dim)      ::force
integer                 ::ibead

do irp1 = 1, path_integral%num_rp
  do irp2 = irp1+1, path_integral%num_rp

    ! ----------------------------------------- !
    ! EDIT NEEDED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !
    ! WE CURRENTLY ASSUME SAME NUMBER OF BEADS  !
    ! FOR DIFFERENT RP. ----------------------- !
    ! ----------------------------------------- !
    do ibead = 1, path_integral%rp(irp1)%num_beads
      force = coulomb_force(path_integral%rp(irp1)%charge/dble(path_integral%rp(irp1)%num_beads), &
        path_integral%rp(irp2)%charge/dble(path_integral%rp(irp2)%num_beads), &
        path_integral%rp(irp1)%x(ibead,1:md%space_dim), path_integral%rp(irp2)%x(ibead,1:md%space_dim) )
      force = force + lj_force(path_integral%rp(irp1)%lj_param(1)/dble(path_integral%rp(irp1)%num_beads), &
        path_integral%rp(irp1)%lj_param(2), path_integral%rp(irp2)%lj_param(2)/dble(path_integral%rp(irp2)%num_beads), &
        path_integral%rp(irp2)%lj_param(2), path_integral%rp(irp1)%x(ibead,1:md%space_dim), &
        path_integral%rp(irp2)%x(ibead,1:md%space_dim) )
      path_integral%rp(irp1)%force(ibead,1:md%space_dim) = path_integral%rp(irp1)%force(ibead,1:md%space_dim) + force
      path_integral%rp(irp2)%force(ibead,1:md%space_dim) = path_integral%rp(irp2)%force(ibead,1:md%space_dim) - force
    enddo
  enddo
enddo


END SUBROUTINE calc_force_pi_pi



! ----------------------------------------------- !
! The velocity-verlet intergrator for the path -- !
! integrator variables. This is the numeric ----- !
! propagator that stays in the cartesian -------- !
! representation. This means a smaller time step- !
! will need to be used when compared to other --- !
! options. -------------------------------------- !
! ----------------------------------------------- !
SUBROUTINE vv_update_pimd
implicit none
integer                 ::irp

!write(*,*)' here we will update the RP positions and velocities '

do irp = 1, path_integral%num_rp
  path_integral%rp(irp)%force = 0.0d0
enddo
write(301,*)' internal force update '
CALL calc_force_pi_internal    !calculcates force on beads from bead spring terms
write(301,*)' external potential force update '
CALL calc_force_pi_model    !calculcates force on beads from external model potential
if (path_integral%num_rp > 1) write(*,*)' do not currently support multiple RP '
if (trim(path_integral%propagation) == 'none') then
  write(302,*)' momentum update '
  CALL update_p_pimd_numeric
  write(303,*)' position update '
  CALL update_x_pimd_numeric
else if (trim(path_integral%propagation) == 'normal') then
  CALL convert_bead_to_normal
  !CALL update_p_pimd_normal
else
  write(*,*)' do not support this type of bead propagation'
  STOP
endif

do irp = 1, path_integral%num_rp
  path_integral%rp(irp)%force = 0.0d0
enddo
CALL calc_force_pi_internal
CALL calc_force_pi_model
CALL update_p_pimd_numeric

END SUBROUTINE vv_update_pimd
! --------------------------------------------- !

! --------------------------------------------- !
! Subroutine that converts the momenta and ---- !
! positions into normal mode coordinates. ----- !
! --------------------------------------------- !
SUBROUTINE convert_bead_to_normal
implicit none
integer                       ::irp
integer                       ::j,k
real(kind=8),allocatable,dimension(:,:)   ::rotation_matrix

do irp = 1, path_integral%num_rp
  allocate(rotation_matrix(path_integral%rp(irp)%num_beads,path_integral%rp(irp)%num_beads))
  do j = 0, path_integral%rp(irp)%num_beads-1
    do k = 0, path_integral%rp(irp)%num_beads-1
      if (k==0) then
        rotation_matrix(j+1,k+1) = dsqrt(1.0d0/path_integral%rp(irp)%num_beads)
      else if (k .lt. path_integral%rp(irp)%num_beads/2 - 1) then
        rotation_matrix(j+1,k+1) = dsqrt(2.0d0/path_integral%rp(irp)%num_beads)*&
          dcos(2.0d0*pi*j*k/dble(path_integral%rp(irp)%num_beads))
      else if (k==path_integral%rp(irp)%num_beads/2) then
        rotation_matrix(j+1,k+1) = dsqrt(1.0d0/path_integral%rp(irp)%num_beads) * (-1.0d0)**j
      else
        rotation_matrix(j+1,k+1) = dsqrt(2.0d0/path_integral%rp(irp)%num_beads)*&
          dsin(2.0d0*pi*j*k/dble(path_integral%rp(irp)%num_beads))
      endif
    enddo
  enddo
  path_integral%rp(irp)%norm_p(:,:) = 0.0d0
  path_integral%rp(irp)%norm_x(:,:) = 0.0d0
  do k = 0, path_integral%rp(irp)%num_beads-1   !over the normal modes   
    do j = 1, path_integral%rp(irp)%num_beads
      path_integral%rp(irp)%norm_p(k,:) = path_integral%rp(irp)%norm_p(k,:)&
        +path_integral%rp(irp)%p(j,:) * rotation_matrix(j,k+1)
      path_integral%rp(irp)%norm_x(k,:) = path_integral%rp(irp)%norm_x(k,:)&
        +path_integral%rp(irp)%x(j,:) * rotation_matrix(j,k+1)
    enddo
  enddo
enddo

END SUBROUTINE convert_bead_to_normal

    


! --------------------------------------------- !
! Subroutine that updates the momentum of the - !
! path integral beads using the numeric velocity!
! verlet scheme.                                !
! p_j(t + delt/2)  <- p_j(t) + delt/2 * F_j     !
! --------------------------------------------- !
SUBROUTINE update_p_pimd_numeric
implicit none
integer                          ::irp
integer                          ::ibead

do irp = 1, path_integral%num_rp
  write(302,*)' initial momentum'
  do ibead = 1, path_integral%rp(irp)%num_beads
      write(302,'(i4,3f15.8)') ibead,path_integral%rp(irp)%p(ibead,1:md%space_dim)
  enddo
  path_integral%rp(irp)%p(:,:) = path_integral%rp(irp)%p(:,:) + md%dt * path_integral%rp(irp)%force(:,:)/2.0d0 
  do ibead = 1, path_integral%rp(irp)%num_beads
      write(302,'(i4,3f15.8)') ibead,path_integral%rp(irp)%p(ibead,1:md%space_dim)
  enddo
enddo

END SUBROUTINE update_p_pimd_numeric
! -------------------------------------------- !


! -------------------------------------------- !
! Subroutine that updates the position of the- !
! path integral beads using the numeric ------ !
! velocity verlet scheme.                      !
! x_j(t+delt) <- x_(t) + delt * p_(t+delt/2)/  !
!      mass_kinetic/N_beads                    !
! -------------------------------------------- !
SUBROUTINE update_x_pimd_numeric
implicit none
integer                          ::irp
integer                          ::ibead

do irp = 1, path_integral%num_rp
  path_integral%rp(irp)%x(:,:) = path_integral%rp(irp)%x(:,:) + md%dt * path_integral%rp(irp)%p(:,:)&
    /path_integral%rp(irp)%kinetic_mass
  do ibead = 1, path_integral%rp(irp)%num_beads
      write(303,'(i4,3f15.8)') ibead,path_integral%rp(irp)%x(ibead,1:md%space_dim)
  enddo
enddo


END SUBROUTINE update_x_pimd_numeric
! -------------------------------------------- !



! -------------------------------------------- !
! Subroutine to calculcate the external force- !
! on the ring polymer. This currently uses the !
! specified model systems to calculate the---- !
! forces. ------------------------------------ !
! -------------------------------------------- !
SUBROUTINE calc_force_pi_model
implicit none     
integer          ::irp,ibead
real(kind=8)     ::omega_ext, mass_ext

if (md%model_no == 1) THEN
  omega_ext = 24.0d0
  mass_ext = 1.0d0
  ! MODEL 1: External HO centered at 0.0
  ! V_ext({q_i}) = 0.5*mass_ext*omega_ext^2 sum_i (q_i)^2
  do irp = 1, path_integral%num_rp
    do ibead = 1, path_integral%rp(irp)%num_beads
      path_integral%rp(irp)%force(ibead,:) = path_integral%rp(irp)%force(ibead,:) &
        - mass_ext/dble(path_integral%rp(irp)%num_beads)& 
        *omega_ext**2 * path_integral%rp(irp)%x(ibead,:)
      write(301,'(i4,3f15.8)') ibead,path_integral%rp(irp)%force(ibead,1:md%space_dim)
    enddo
    write(301,*)
  enddo
else
  write(*,*)' this potential number is not a defined model potential'
endif

END SUBROUTINE calc_force_pi_model
! ------------------------------------------- !


! ---------------------------------------------- !
!Subroutine that calculcates the internal forces !
!on the ring-polymer beads due to the bead-bead- !
!spring terms. --------------------------------- !
! ---------------------------------------------- !
SUBROUTINE calc_force_pi_internal     !THE BEAD-BEAD PREFACTOR IS NOT RIGHT!!!
implicit none

integer           ::irp,ibead
integer           ::ibead_minus,ibead_plus

do irp = 1, path_integral%num_rp
  do ibead = 1, path_integral%rp(irp)%num_beads
    ibead_minus = ibead-1
    if (ibead == 1) ibead_minus = path_integral%rp(irp)%num_beads
    ibead_plus = ibead+1
    if(ibead == path_integral%rp(irp)%num_beads) ibead_plus = 1
    path_integral%rp(irp)%force(ibead,:) = path_integral%rp(irp)%force(ibead,:) &
      - (path_integral%rp(irp)%harm_mass/dble(path_integral%rp(irp)%num_beads))&
      /(path_integral%rp(irp)%betaN**2*path_integral%rp(irp)%num_beads)&     !CHECK!!!
      *(2.0d0*path_integral%rp(irp)%x(ibead,:) - path_integral%rp(irp)%x(ibead_minus,:) &
      - path_integral%rp(irp)%x(ibead_plus,:) )
    write(301,'(i4,3f15.8)') ibead,path_integral%rp(irp)%force(ibead,1:md%space_dim)
  enddo
  write(301,*)
enddo

END SUBROUTINE calc_force_pi_internal
! -------------------------------------------- !


! -------------------------------------------- !
! Reads in the variables from the md_run.inp - !
! input file. Generally concerned with the --- !
! details of the MD run.                       !
! nsteps = number of MD steps                  !
! freq_print = how often to save the positions !
!       and momenta to the traj array--------- !
! dt = time step in fs ----------------------- !
! space_dim = dimension of space (1D, 2D, 3D)- !
! cp2k_input_file = absolute path to the cp2k- !
!       input file if using CP2K to calculcate !
!       forces.                                !
! cp2k_output_file = absolute path to the cp2k !
!       output file.                           !
! -------------------------------------------- !
SUBROUTINE read_in_md
implicit none


open(101,file='md_run.inp')
read(101,*) md%nsteps
read(101,*) md%freq_print
read(101,*) md%dt
read(101,*) md%space_dim
read(101,*) md%obtain_force
read(101,*) md%model_no
read(101,*) md%temp
close(101)

! -------------------------------------------- !
! WE ASSUME THAT TIME WAS READ IN IN FS. HERE  !
! WE CONVERT THIS TIME INTO ATOMIC UNITS. ---- !
! -------------------------------------------- !
!convert time into atomic units

!md%dt = md%dt/(2.418884326e-2)
md%dt = md%dt * fstota
write(*,*)'atomic time is ', md%dt

ALLOCATE(md%total_energy(0:md%nsteps))
md%total_energy = 0.0d0

END SUBROUTINE read_in_md
! --------------------------------------------- !



! ---------------------------------------------- !
! Subroutine to initalize the velocoties of the- !
! beads for each RP. They are random numbers---- !
! generaed by a normal distribution to mimic the !
! Maxwell-Boltzmann with mean = 0 and ---------- !
! sigma = sqrt(1/kinetic_mass*betaN) in each --- !
! direction.------------------------------------ !
SUBROUTINE initialize_bead_velocities
implicit none

integer                  ::irp, ibead, i
real(kind=8)             ::rnd_vel

do irp = 1, path_integral%num_rp

! ---------------------------------------------- !
! ALLOCATES THE INITALIZE_VELOCITY ARRAYS ------ !
! ---------------------------------------------- !
  ALLOCATE(path_integral%rp(irp)%initialized_velocities(path_integral%rp(irp)%num_beads, &
md%space_dim))
  path_integral%rp(irp)%initialized_velocities = 0.0d0

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


! --------------------------------------------- !
! FILE THAT WILL PRINT THE INITALIZE BEAD ----- !
! POSITIONS AND VELOCITIES TO AN XYZ FILE. ---- !
! CANNOT USE THIS IN THE NEW LOCATION. -------- !
! THE POSITIONS ARE NOT YET INITALIZED. ------- !
! EDIT NEEDED. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !
! --------------------------------------------- !
!write(*,*) 'print the velocities-BAJ error'
!open(203,file='check_initialize.dat')
!do irp = 1, path_integral%num_rp
!  do ibead = 1, path_integral%rp(irp)%num_beads
!    write(*,*) 'wirte the velocities-BAJ error'
!    write(203,'(A4, 6f15.8)') 'H', path_integral%rp(irp)%initialized_positions(ibead,1:md%space_dim), &
!path_integral%rp(irp)%initialized_velocities(ibead,1:md%space_dim)
!    write(*,*) 'done wirte the velocities-BAJ error'
!  enddo
!  write(203,*)
!enddo

END SUBROUTINE initialize_bead_velocities
! --------------------------------------------- !

SUBROUTINE initialize_classical_velocities
implicit none

real(kind=8)         ::rnd_vel
integer              ::iatom,i

ALLOCATE(classical%i_velocity(classical%natoms,md%space_dim))
classical%i_velocity = 0.0d0

do iatom = 1, classical%natoms
  do i = 1, md%space_dim
    CALL rnd_gauss(0.0d0,dsqrt(1.0d0/( classical%beta * classical%mass(iatom) )), rnd_vel)
    classical%i_velocity(iatom,i) = rnd_vel
  enddo
enddo

END SUBROUTINE initialize_classical_velocities
  



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
SUBROUTINE pimd_initialize
implicit none

integer                 ::i,j
integer                 ::ibead
character(len=6)        ::char_h


!For now, will not try any "fancy" reads. Modify this later.
! ------------------------------------------------------- !
! READ IN THE GENERAL INFORMATION ABOUT EACH OF THE RING- !
! POLYMERS IN THE SYSTEM. SEE NOTES IN PIMD.INP FOR MORE- !
! INFO ON EACH VARIABLE. -------------------------------- !
! ------------------------------------------------------- ! 
OPEN(101,FILE='pimd.inp')
read(101,*) path_integral%num_rp
ALLOCATE(path_integral%rp(path_integral%num_rp))
do i = 1, path_integral%num_rp
  read(101,*) path_integral%rp(i)%num_beads
  read(101,*) path_integral%rp(i)%particle_type
  read(101,*) path_integral%rp(i)%kinetic_mass
  read(101,*) path_integral%rp(i)%harm_mass
  read(101,*) path_integral%rp(i)%charge
  read(101,*) path_integral%rp(i)%lj_param(1:2)
  read(101,*) path_integral%rp(i)%coord
  path_integral%rp(i)%kinetic_mass = path_integral%rp(i)%harm_mass/dble(path_integral%rp(i)%num_beads)
enddo
CLOSE(101)

! -------------------------------------------------------- !
! READ IN THE PARAMETERS FOR THE PIMD SPECIFIC VARIABLES.- !
! IN PARTICULAR, THE INFORMATION ON THE INITALIZATION AND- !
! INITIAL VALUES FOR THE RING POLYMERS. ------------------ !
! -------------------------------------------------------- !
OPEN(101,FILE='system_parameters.inp') 
read(101,*) path_integral%initialize
read(101,*) path_integral%debug_initialize
read(101,*) path_integral%init_nsteps
read(101,*) path_integral%init_print_nsteps
read(101,*) path_integral%mc_step
read(101,*) path_integral%propagation
CLOSE(101)

! -------------------------------------------------------- !
! ALLOCATES THE INITIAL POSITION ARRAYS FOR ALL BEAD ----- !
! AND CENTROID POSITIONS --------------------------------- !
! READS IN EITHER THE INITIAL COORDINATES OR INITIAL ----- !
! CENTROID POSITIONS FROM FILES COORD_RP-NBEAD.XYZ OR ---- !
! CENTROID_RP-NBEAD.XYZ. --------------------------------- !
! -------------------------------------------------------- !
do i = 1, path_integral%num_rp
  ALLOCATE(path_integral%rp(i)%i_position(path_integral%rp(i)%num_beads,md%space_dim))
  path_integral%rp(i)%i_position = 0.0d0
  ALLOCATE(path_integral%rp(i)%i_centroid(md%space_dim))
  path_integral%rp(i)%i_centroid = 0.0d0
enddo

do i = 1, path_integral%num_rp
  if(path_integral%rp(i)%coord .EQV. .TRUE.) THEN
    OPEN(101,FILE='coord_rp-'//char(i+48)//'.xyz')
    do ibead = 1, path_integral%rp(i)%num_beads
      read(101,*) char_h, path_integral%rp(i)%i_position(ibead,:)
    enddo
    CLOSE(101)
    do j = 1, md%space_dim
      path_integral%rp(i)%i_centroid(j) = sum(path_integral%rp(i)%i_position(:,j))
      write(*,*)' the initial calculcated i_centroid is ', path_integral%rp(i)%i_centroid(j)
    enddo
  else
    OPEN(101,FILE='centroid_rp-'//char(i+48)//'.xyz')
      read(101,*) path_integral%rp(i)%i_centroid(:)
    CLOSE(101)
  endif
enddo

! ------------------------------------------------------ !
! CONVERTS INPUT PARAMETERS INTO ATOMIC UNITS. HERE I -- !
! ASSUME THAT THE INPUTS ARE IN THE FOLLOWING UNITS: --- !
! MASSES: AMU (UNLESS PARTICLE=='SPECIAL') ------------- !
! DT: FS
! LENNARD-JONES PARAMETERS: K_e (EPS), A (SIG)
! POSITIONS: A
! ALSO ASSUME THAT CHARGE RIGHT NOW IS IN ATOMIC UNITS!!!!
! ------------------------------------------------------ !
do i = 1, path_integral%num_rp
  path_integral%rp(i)%i_position = path_integral%rp(i)%i_position * AtoBohr
  path_integral%rp(i)%i_centroid = path_integral%rp(i)%i_centroid * AtoBohr
  if (trim(path_integral%rp(i)%particle_type) == 'ATOM') then
    path_integral%rp(i)%harm_mass = path_integral%rp(i)%harm_mass * amutome
    path_integral%rp(i)%kinetic_mass = path_integral%rp(i)%kinetic_mass * amutome
  endif
  path_integral%rp(i)%lj_param(1) = path_integral%rp(i)%lj_param(1) * KtoH
  path_integral%rp(i)%lj_param(2) = path_integral%rp(i)%lj_param(2) * AtoBohr
enddo


! ------------------------------------------------------- !
! CALCULATES BETA_N VALUES FOR EACH OF THE RING-POLYMERS  !
! ------------------------------------------------------- !
do i = 1, path_integral%num_rp
  path_integral%rp(i)%betaN = (1.0d0/(kb * md%temp)) / dble(path_integral%rp(i)%num_beads)
  path_integral%rp(i)%beta = (1.0d0/(kb * md%temp))
enddo
write(*,*)' num_beads is ', path_integral%rp(1)%num_beads
write(*,*)' md temp is ', md%temp
write(*,*)' betaN is ', path_integral%rp(1)%betaN

! ------------------------------------------------------- !
! INITALIZE THE RANDOM NUMBER GENERATOR. ---------------- !
! ------------------------------------------------------- !
write(*,*)' init seed - BAJ error'
call init_rnd_seed()

! ------------------------------------------------------- !
! INITALIZE THE VELOCITIES OF THE BEADS ----------------- !
! ------------------------------------------------------- !
!initalize velocities of the bea
write(*,*)' init vel-BAJ error'
CALL initialize_bead_velocities

open(222,file='initial_bead_velocities.dat')
do i = 1, path_integral%num_rp
  write(222,*)' particle number ', i
  do j = 1, path_integral%rp(i)%num_beads
    write(222,'(3f10.6)') path_integral%rp(i)%initialized_velocities(j,:)
  enddo
enddo
close(222)
!open(222,file='initial_bead_velocities.dat')
!do i = 1, path_integral%num_rp
!  read(222,*)
!  do j = 1, path_integral%rp(i)%num_beads
!    read(222,*) path_integral%rp(i)%initialized_velocities(j,:)
!  enddo
!enddo
!close(222)

END SUBROUTINE pimd_initialize
! --------------------------------------------------- !


! ------------------------------------------------------ !
! Subroutine to initalize the ring-polymer configuration !
! using a series of MC steps in the "frozen" external -- !
! potential that is either the model potential or the -- !
! CP2K geneated potential. The centroid of the ring ---- !
! polymer is confined to the original location of the -- !
! centroid specified in the input files (either -------- !
! centroid_rp.xyz or coord_rp.xyz ---------------------- !
! ------------------------------------------------------ !
SUBROUTINE initialize_mc
implicit none

integer                         ::i, j
integer                         ::irp
real(kind=8)                    ::rnd
integer                         ::kbead,kdim
integer                         ::accept, iprint, istep
real(kind=8)                    ::rnd1
real(kind=8)                    ::pot_old, pot_new
real(kind=8),dimension(:)       ::centroid_n(md%space_dim), centroid_check(md%space_dim)

if (path_integral%debug_initialize) then
  write(*,*)' will print off ', int(path_integral%init_nsteps/path_integral%init_print_nsteps)+1
  do irp = 1, path_integral%num_rp
    ALLOCATE(path_integral%rp(irp)%init_traj(path_integral%rp(irp)%num_beads, &
      0:path_integral%init_nsteps/path_integral%init_print_nsteps, md%space_dim))
    path_integral%rp(irp)%init_traj = 0.0d0
  enddo
endif




do irp = 1, path_integral%num_rp
  ALLOCATE(path_integral%rp(irp)%x_old(path_integral%rp(irp)%num_beads, md%space_dim))
  ALLOCATE(path_integral%rp(irp)%x_new(path_integral%rp(irp)%num_beads, md%space_dim))
  path_integral%rp(irp)%x_old = 0.0d0
  path_integral%rp(irp)%x_new = 0.0d0
  write(*,*)' the centroid is ', path_integral%rp(irp)%i_centroid(:)
  if (.NOT. path_integral%rp(irp)%coord) THEN
    do i = 1, path_integral%rp(irp)%num_beads
      do j = 1, md%space_dim
        path_integral%rp(irp)%x_old(i,j) = path_integral%rp(irp)%i_centroid(j) + path_integral%mc_step*&
          dcos(2.0d0*pi*i/dble(path_integral%rp(irp)%num_beads))
      enddo
    enddo
  else
    do i = 1, path_integral%rp(irp)%num_beads
      path_integral%rp(irp)%x_old(i,:) = path_integral%rp(irp)%i_position(i,:)
    enddo
  endif
  
  if (path_integral%debug_initialize) THEN
    do i = 1, path_integral%rp(irp)%num_beads
      path_integral%rp(irp)%init_traj(i,0,:) = path_integral%rp(irp)%x_old(i,:)
    enddo
  endif
  path_integral%rp(irp)%x_new(i,:) = path_integral%rp(irp)%x_old(i,:)
enddo

pot_old = 0.0d0
call calc_potential_mc(pot_old)

iprint = 0
accept = 0.0d0
open(201,file='potential_energy_mc.dat')
do istep = 1, path_integral%init_nsteps
  call random_number(rnd1)
  irp = int(rnd1*path_integral%num_rp)+1
  write(*,*)' irp is ', irp
  path_integral%rp(irp)%x_new = path_integral%rp(irp)%x_old
  do i = 1, path_integral%rp(irp)%num_beads
    call random_number(rnd1)
    kbead = int(rnd1*path_integral%rp(irp)%num_beads)+1
    write(*,*)' rnd1 is ', rnd1
    write(*,*)' kbead ', kbead
    call random_number(rnd1)
    kdim = int(rnd1*md%space_dim)+1
    write(*,*)' rnd2 is ', rnd1
    write(*,*)' kdim ', kdim
   
    
    call random_number(rnd1)
    path_integral%rp(irp)%x_new(kbead,kdim) = path_integral%rp(irp)%x_new(kbead,kdim) + path_integral%mc_step*dble(rnd1-0.5d0)
    write(*,*)' rnd3 is ', rnd1
    write(*,*)' x_new is ',  path_integral%rp(irp)%x_new(kbead,kdim)
  enddo

  !write(*,*)' the initial centroid is ', path_integral%rp(irp)%i_centroid(i)
  do i = 1, md%space_dim
    centroid_n(i) = sum(path_integral%rp(irp)%x_new(:,i))/dble(path_integral%rp(irp)%num_beads)
    !do j = 1, path_integral%rp(irp)%num_beads
    !  write(*,*)' so the original RP is ', x_new(j,i)
    !enddo
    !write(*,*)' on step ', istep
    !write(*,*)' new centroid is ', centroid_n(i)/dble(path_integral%rp(irp)%num_beads)


    path_integral%rp(irp)%x_new(:,i) = path_integral%rp(irp)%x_new(:,i) - centroid_n(i) + path_integral%rp(irp)%i_centroid(i)


    do j = 1, path_integral%rp(irp)%num_beads
      write(*,*)' so the shifted RP is ', path_integral%rp(irp)%x_new(j,i)
    enddo
    centroid_check(i) = sum(path_integral%rp(irp)%x_new(:,i))/dble(path_integral%rp(irp)%num_beads)
    write(*,*)' to check shift, centroid is ', centroid_check(i)
  enddo

  call calc_potential_mc(pot_new)

  call random_number(rnd1)
  write(*,*)' rnd1 is ', rnd1
  if (rnd1 .LT. exp((-1.0d0)*path_integral%rp(irp)%betaN*(pot_new-pot_old)) ) THEN
    write(*,*)' move accepted '
    accept = accept+1
    path_integral%rp(irp)%x_old = path_integral%rp(irp)%x_new
    pot_old = pot_new
  endif

  if (path_integral%debug_initialize) THEN
    if(mod(istep,path_integral%init_print_nsteps)==0) THEN
      iprint = iprint+1
      !write(*,*)' iprint is ', iprint, pot_old
      write(201,'(2i12,f20.6)') istep, iprint, pot_old
      do i = 1, path_integral%rp(irp)%num_beads
        path_integral%rp(irp)%init_traj(i,iprint,:) = path_integral%rp(irp)%x_old(i,:)
      enddo
    endif
  endif

enddo
close(201)

write(*,*)' accept is ', accept
write(*,*)' ratios is ', dble(accept)/dble(path_integral%init_nsteps)
  

!writes the trajectory file
do irp = 1, path_integral%num_rp
  open(201,file='initialize_traj_mc'//char(irp+48)//'.xyz')
  write(*,*)' print the mc traj file '
  do i = 0, path_integral%init_nsteps/path_integral%init_print_nsteps
    write(201,*) path_integral%rp(irp)%num_beads
    write(201,*)' i = ', i
    do j = 1, path_integral%rp(irp)%num_beads
      write(201,'(a2, 3f12.4)') 'H', path_integral%rp(irp)%init_traj(j,i,:)     !find way to change...
    enddo
  enddo
  close(201) 
enddo

do irp = 1, path_integral%num_rp
  ALLOCATE(path_integral%rp(irp)%initialized_positions(path_integral%rp(irp)%num_beads,md%space_dim))
  path_integral%rp(irp)%initialized_positions = 0.0d0
  do i = 1, path_integral%rp(irp)%num_beads
    path_integral%rp(irp)%initialized_positions(i,:) = path_integral%rp(irp)%init_traj(i,&
      path_integral%init_nsteps/path_integral%init_print_nsteps,:)
  enddo

  DEALLOCATE(path_integral%rp(irp)%x_old, path_integral%rp(irp)%x_new)
  DEALLOCATE(path_integral%rp(irp)%init_traj)
enddo

END SUBROUTINE initialize_mc
! --------------------------------------------------- !




! --------------------------------------------------- !
! Subroutine that evaluates the potential of the ring !
! polymer at a given ring-polymer configuration. 
! REWRITE THIS SUBROUTINE TO BE MORE GENERAL. INIDIVUALLY !
! CALL THE INTERNAL VERSUS EXTERNAL POTENTIAL AND ALSO    !
! CHANGE SO THAT CP2K POTENTIALS CAN BE CALCULCATED AS    !
! WELL. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! --------------------------------------------------- !
SUBROUTINE calc_potential_mc(pot)
implicit none

integer                       ::irp
real(kind=8)                  ::pot

integer              ::nbead
integer              ::i, iplus

real(kind=8)         ::omega
real(kind=8)         ::omega_ext, mass_ext

real(kind=8)         ::inter_pot, inter_rp_pot


pot = 0.0d0

!calc spring terms

do irp = 1, path_integral%num_rp
  nbead = path_integral%rp(irp)%num_beads
  omega = 1.0d0/path_integral%rp(irp)%betaN
  do i = 1, nbead
    iplus = i+1
    if (i==nbead) iplus = 1
    pot = pot + 0.5d0*path_integral%rp(irp)%harm_mass*omega**2/dble(path_integral%rp(irp)%num_beads)& 
      * dot_product(path_integral%rp(irp)%x_new(i,:) &
      -path_integral%rp(irp)%x_new(iplus,:),path_integral%rp(irp)%x_new(i,:)-path_integral%rp(irp)%x_new(iplus,:))
  enddo
enddo

!write(*,*)' the harmonic mass is '
!write(*,*) path_integral%rp(irp)%harm_mass
!write(*,*)' bead-bead prefactor is '
!write(*,*) 0.5d0*(path_integral%rp(irp)%harm_mass/dble(path_integral%rp(irp)%num_beads))/&
!  (path_integral%rp(irp)%betaN**2*&
!  dble(path_integral%rp(irp)%num_beads))

!write(*,*)' the bead-potential is ', pot

if(trim(md%obtain_force) == 'model') THEN
  if (md%model_no == 1) THEN
  !the model potential V(x) =  0.5 omega_ext**2 x^2 in atomic units
    omega_ext = 24.0d0
    mass_ext = 1.0d0
    do irp = 1, path_integral%num_rp
      do i = 1, nbead
        pot = pot + 0.5d0*mass_ext/dble(path_integral%rp(irp)%num_beads)&                
            *omega_ext**2*dot_product(path_integral%rp(irp)%x(i,:),path_integral%rp(irp)%x(i,:))
      enddo
    enddo
   else
     write(*,*)' this potential number is not a defined model potential'
   endif
else if (trim(md%obtain_force) == 'cp2k') THEN
  CALL calc_interaction_potential_mc(inter_pot,inter_rp_pot)
  pot = pot + inter_pot + inter_rp_pot
  write(*,*)' the potential_mc is '
  write(*,*) pot, inter_pot, inter_rp_pot
else
  write(*,*)' this method of potential evalulation is not currently implemented '
endif


END SUBROUTINE calc_potential_mc
! ----------------------------------------------------- !



! ----------------------------------------------------- !
! SUBROUTINE TO CALCULATE THE INTERACTION POTENTIAL --- !
! BETWEEN THE BEADS IN THE RING POLYMER AND THE ------- !
! CLASSICAL ATOMS. ------------------------------------ !
! ----------------------------------------------------- !
SUBROUTINE calc_interaction_potential(inter_pot, inter_rp_pot)
implicit none

real(kind=8)              ::inter_pot, inter_rp_pot
integer                       ::irp, jrp, iclassical, ibead

!1. calculates the coulombic potential using
!       i. x(irp) (bead positions of "specified" ring polymer
!       ii. classical%i_positions (classical atom positions)
!       iii. classical%ff_atom_charges
!       iv. path_integral%rp(irp)%charge/nbead

!EDIT NEEDED!!!!!!!!!!!!!!!!!!!!!!!!!!! !
! GOING TO CURRENTLY ASSUME THAT THERE IS !
! THE SAME NUMBER OF BEADS IN EACH RP     !

inter_pot = 0.0d0
inter_rp_pot = 0.0d0
do irp = 1, path_integral%num_rp
  do ibead = 1, path_integral%rp(irp)%num_beads
    do iclassical = 1, classical%natoms
      inter_pot = inter_pot + coulomb_pot( classical%ff_atom_charges(iclassical), &
        path_integral%rp(irp)%charge/dble(path_integral%rp(irp)%num_beads), classical%x(iclassical,:), &
        path_integral%rp(irp)%x(ibead,:))
      inter_pot = inter_pot + lj_pot( classical%ff_atom_lj(iclassical,1), classical%ff_atom_lj(iclassical,2), &
        path_integral%rp(irp)%lj_param(1)/(dble(path_integral%rp(irp)%num_beads)**2), &
        path_integral%rp(irp)%lj_param(2), classical%x(iclassical,:),path_integral%rp(irp)%x(ibead,:))
        write(*,*)' inter_pot is ', ibead, iclassical, inter_pot
      do jrp = irp+1, path_integral%num_rp
        inter_rp_pot = inter_rp_pot + coulomb_pot( &
          path_integral%rp(irp)%charge/dble(path_integral%rp(irp)%num_beads), &
          path_integral%rp(jrp)%charge/dble(path_integral%rp(jrp)%num_beads), &
          path_integral%rp(irp)%x(ibead,:), path_integral%rp(jrp)%x(ibead,:))
        inter_rp_pot = inter_rp_pot + lj_pot( &
          path_integral%rp(irp)%lj_param(1)/dble(path_integral%rp(irp)%num_beads), &
          path_integral%rp(irp)%lj_param(2)/dble(path_integral%rp(irp)%num_beads), &
          path_integral%rp(jrp)%lj_param(1)/dble(path_integral%rp(jrp)%num_beads), &
          path_integral%rp(jrp)%lj_param(2)/dble(path_integral%rp(jrp)%num_beads), &
          path_integral%rp(irp)%x(ibead,:), path_integral%rp(jrp)%x(ibead,:) )
      enddo
    enddo
  enddo
enddo

END SUBROUTINE calc_interaction_potential

! ----------------------------------------------------- !
! SUBROUTINE TO CALCULATE THE INTERACTION POTENTIAL --- !
! BETWEEN THE BEADS IN THE RING POLYMER AND THE ------- !
! CLASSICAL ATOMS. ------------------------------------ !
! ----------------------------------------------------- !
SUBROUTINE calc_interaction_potential_mc(inter_pot, inter_rp_pot)
implicit none

real(kind=8)              ::inter_pot, inter_rp_pot
integer                       ::irp, jrp, iclassical, ibead

!1. calculates the coulombic potential using
!       i. x(irp) (bead positions of "specified" ring polymer
!       ii. classical%i_positions (classical atom positions)
!       iii. classical%ff_atom_charges
!       iv. path_integral%rp(irp)%charge/nbead

!EDIT NEEDED!!!!!!!!!!!!!!!!!!!!!!!!!!! !
! GOING TO CURRENTLY ASSUME THAT THERE IS !
! THE SAME NUMBER OF BEADS IN EACH RP     !

inter_pot = 0.0d0
inter_rp_pot = 0.0d0
do irp = 1, path_integral%num_rp
  do ibead = 1, path_integral%rp(irp)%num_beads
    do iclassical = 1, classical%natoms
      inter_pot = inter_pot + coulomb_pot( classical%ff_atom_charges(iclassical), &
        path_integral%rp(irp)%charge/dble(path_integral%rp(irp)%num_beads), classical%x(iclassical,:), &
        path_integral%rp(irp)%x_new(ibead,:))
      inter_pot = inter_pot + lj_pot( classical%ff_atom_lj(iclassical,1), classical%ff_atom_lj(iclassical,2), &
        path_integral%rp(irp)%lj_param(1)/dble(path_integral%rp(irp)%num_beads), &
        path_integral%rp(irp)%lj_param(2), classical%x(iclassical,:),path_integral%rp(irp)%x_new(ibead,:))
      do jrp = irp+1, path_integral%num_rp
        inter_rp_pot = inter_rp_pot + coulomb_pot( &
          path_integral%rp(irp)%charge/dble(path_integral%rp(irp)%num_beads), &
          path_integral%rp(jrp)%charge/dble(path_integral%rp(jrp)%num_beads), &
          path_integral%rp(irp)%x_new(ibead,:), path_integral%rp(jrp)%x_new(ibead,:))
        inter_rp_pot = inter_rp_pot + lj_pot( &
          path_integral%rp(irp)%lj_param(1)/dble(path_integral%rp(irp)%num_beads), &
          path_integral%rp(irp)%lj_param(2)/dble(path_integral%rp(irp)%num_beads), &
          path_integral%rp(jrp)%lj_param(1)/dble(path_integral%rp(jrp)%num_beads), &
          path_integral%rp(jrp)%lj_param(2)/dble(path_integral%rp(jrp)%num_beads), &
          path_integral%rp(irp)%x_new(ibead,:), path_integral%rp(jrp)%x_new(ibead,:) )
      enddo
    enddo
  enddo
enddo

END SUBROUTINE calc_interaction_potential_mc





! ------------------------------------------------------ !
! Subroutine that will deallocate/remove all of the pimd !
! structures to prevent any memory leakage, etc. ------- !
! ------------------------------------------------------ !
SUBROUTINE pimd_destory
implicit none

integer                     ::i


do i = 1, path_integral%num_rp
  DEALLOCATE(path_integral%rp(i)%i_position)
  DEALLOCATE(path_integral%rp(i)%i_centroid)
  DEALLOCATE(path_integral%rp(i)%initialized_positions)
  DEALLOCATE(path_integral%rp(i)%initialized_velocities)
  DEALLOCATE(path_integral%rp(i)%traj)
  DEALLOCATE(path_integral%rp(i)%force)
  DEALLOCATE(path_integral%rp(i)%p)
  DEALLOCATE(path_integral%rp(i)%x)
enddo

DEALLOCATE(path_integral%rp)

END SUBROUTINE pimd_destory
! ----------------------------------------------------- !


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

SUBROUTINE update_normal_coord(dt)
implicit none

real(kind=8)     ::omega_n, omega_k
real(kind=8)     ::dt
integer          ::irp, ibead, nbeads, i_dim

real(kind=8),dimension(2,2)    ::free_rp_prop
real(kind=8),dimension(2)      ::old_norm, new_norm


do irp = 1, path_integral%num_rp
  !omega_n = (path_integral%rp(irp)%harm_mass/dble(path_integral%rp(irp)%num_beads))/(path_integral%rp(irp)%betaN**2*path_integral%rp(irp)%num_beads)   !CHECK!!
  !omega_n = dsqrt(path_integral%rp(irp)%num_beads/(path_integral%rp(irp)%betaN))   !CHECK!!
  omega_n = dsqrt(path_integral%rp(irp)%harm_mass/(path_integral%rp(irp)%kinetic_mass*&
    dble(path_integral%rp(irp)%num_beads)))/path_integral%rp(irp)%betaN
  
  nbeads = path_integral%rp(irp)%num_beads
  do ibead = 1, nbeads
    omega_k =2.0d0*omega_n*dsin( (ibead-1)*pi/dble(nbeads) )
    if (ibead==1) then
      free_rp_prop(1,1) = 1.0d0
      free_rp_prop(2,2) = 1.0d0
      free_rp_prop(1,2) = 0.0d0
      free_rp_prop(2,1) = dt/path_integral%rp(irp)%kinetic_mass   !CHECK-the NBEADS
    else
      free_rp_prop(1,1) = dcos(omega_k*dt)
      free_rp_prop(2,2) = free_rp_prop(1,1)
      free_rp_prop(1,2) = -path_integral%rp(irp)%kinetic_mass*omega_k*dsin(omega_k*dt)              !CHECK-the NBEADS
      free_rp_prop(2,1) = (1.0d0/(path_integral%rp(irp)%kinetic_mass * omega_k))*dsin(omega_k*dt)
    endif

    do i_dim = 1, md%space_dim
      old_norm(1) = path_integral%rp(irp)%norm_p(ibead,i_dim)
      old_norm(2) = path_integral%rp(irp)%norm_x(ibead,i_dim)

      new_norm = matmul( free_rp_prop, old_norm )

      path_integral%rp(irp)%norm_p(ibead,i_dim) = new_norm(1)
      path_integral%rp(irp)%norm_x(ibead,i_dim) = new_norm(2)
    enddo
  enddo
enddo


END SUBROUTINE update_normal_coord




END MODULE pimd_mod
