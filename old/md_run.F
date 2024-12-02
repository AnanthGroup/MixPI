MODULE md_run
use variable_types
implicit none

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
SUBROUTINE run_md(classical, path_integral, md, para_var)
implicit none

TYPE(path_integral_env)                ::path_integral
TYPE(classical_type)                   ::classical
TYPE(md_info)                           ::md
TYPE(parallel_env)                     ::para_var
logical         ::model_system = .FALSE.
integer         ::irp, istep, iprint, ibead, iclassical, icount
integer         ::i,j

integer         ::irate, itime_istart, itime_iend
integer         ::trate, itime_tstart, itime_tend
integer         ::crate, itime_cstart, itime_cend
real(kind=8)    ::run_temp, run_tot, run_pot, run_kin
real(kind=8),dimension(1:3)     ::com, com_check


! ------------------------------------------------------ !
! OPEN CONTINUALLY RUNNING FILES ----------------------- !
! ------------------------------------------------------ !
!trajectory files
write(*,*)' start md_run '
open(676, file='trajectory_run_full.dat')
open(677, file='trajectory_run_com.dat')
open(678, file='trajectory_run_centroid.dat')
!energy files
open(203,file='energy_acutal.dat')
open(204,file='energy_virials.dat')
open(205,file='energy_average.dat')
run_kin = 0.0d0
run_pot = 0.0d0
run_temp = 0.0d0
icount = 0
if (trim(adjustl(md%print_level)) == 'high') open(801,file='forces.dat')


! ------------------------------------------- !
! ALLOCATE AND SET-UP PATH-INTEGRAL VARIABLES !
! AND ARRAYS. ------------------------------- !
! ------------------------------------------- !
do irp = 1, path_integral%num_rp
  ALLOCATE(path_integral%rp(irp)%force(path_integral%rp(irp)%num_beads,md%space_dim))
  ALLOCATE(path_integral%rp(irp)%p(path_integral%rp(irp)%num_beads,md%space_dim))
  ALLOCATE(path_integral%rp(irp)%x(path_integral%rp(irp)%num_beads,md%space_dim))
  ALLOCATE(path_integral%rp(irp)%centroid(md%space_dim))

  path_integral%rp(irp)%force = 0.0d0
  path_integral%rp(irp)%p = 0.0d0
  path_integral%rp(irp)%x = 0.0d0
  path_integral%rp(irp)%centroid = 0.0d0

  path_integral%rp(irp)%p = path_integral%rp(irp)%kinetic_mass &
    *path_integral%rp(irp)%initialized_velocities
  path_integral%rp(irp)%x = path_integral%rp(irp)%initialized_positions
  do i = 1, md%space_dim
    path_integral%rp(irp)%centroid(i) = sum(path_integral%rp(irp)%x(:,i))&
    /dble(real(path_integral%rp(irp)%num_beads))
  enddo
  if (trim(adjustl(md%print_level)) == 'high') then
    ALLOCATE(path_integral%rp(irp)%traj(0:md%nsteps/md%freq_print,&
      path_integral%rp(irp)%num_beads,2*md%space_dim))
    path_integral%rp(irp)%traj = 0.0d0
    do ibead = 1, path_integral%rp(irp)%num_beads
      path_integral%rp(irp)%traj(0,ibead,1:md%space_dim)=path_integral%rp(irp)%x(ibead,:)
      path_integral%rp(irp)%traj(0,ibead,md%space_dim+1:2*md%space_dim)&
         =path_integral%rp(irp)%p(ibead,1:md%space_dim)/path_integral%rp(irp)%kinetic_mass
    enddo
  endif

  if (trim(md%propagation) == 'normal_fftw' .OR. trim(md%propagation) == 'normal_manual') then
    ALLOCATE(path_integral%rp(irp)%norm_p(path_integral%rp(irp)%num_beads,md%space_dim))
    ALLOCATE(path_integral%rp(irp)%norm_x(path_integral%rp(irp)%num_beads,md%space_dim))
    path_integral%rp(irp)%norm_p = 0.0d0
    path_integral%rp(irp)%norm_x = 0.0d0
  endif
  ALLOCATE(path_integral%rp(irp)%com_pos(0:md%nsteps/md%freq_print, md%space_dim))
  ALLOCATE(path_integral%rp(irp)%com_vel(0:md%nsteps/md%freq_print, md%space_dim))
enddo
path_integral%potential = 0.0d0
path_integral%kinetic = 0.0d0


! ------------------------------------------------------!
! INITIALIZES THE FORCE ARRAYS AND GENERATES THE -------!
! ENERGY AND FORCES FOR THE INITIAL POSITIONS ----------!
! ------------------------------------------------------!

CALL update_force_env(classical, path_integral, md)
do irp = 1, path_integral%num_rp
  path_integral%rp(irp)%force = 0.0d0
enddo
classical%force = 0.d0
CALL calc_force_cp2k(classical, path_integral, md)
CALL calc_force_pi(path_integral, md, 0)

CALL calc_total_energy(classical, path_integral, md, 0)
iprint = 0


! -------------------------------------------------------!
! PRINTS OUT THE ENERGIES AND TRAJECTORY FILES FOR THE --!
! INITIAL STRUCTURES ------------------------------------!
! -------------------------------------------------------!
! energy files
write(203,'(i5,f16.2,7f20.10)')0, 0*md%dt/fstota,classical%kinetic,classical%potential, path_integral%kinetic, &
    path_integral%potential, md%kinetic_virial, md%potential_virial, md%kinetic_primitive
if (classical%natoms+path_integral%total_beads <= 2) then
  write(204,'(i5,2f16.2,4f20.10)')0, 0*md%dt/fstota,(classical%kinetic+path_integral%kinetic)/&
    ((real(md%space_dim)*(classical%natoms+path_integral%total_beads)-5)*kb/2.0d0), &
    md%kinetic_virial, md%potential_virial, md%kinetic_primitive, md%kinetic_virial + md%potential_virial
else 
  write(204,'(i5,2f16.2,4f20.10)')0, 0*md%dt/fstota,(classical%kinetic+path_integral%kinetic)/&
    ((real(md%space_dim)*(classical%natoms+path_integral%total_beads)-(3+md%n_constraints))*kb/2.0d0), &
    md%kinetic_virial, md%potential_virial, md%kinetic_primitive, md%kinetic_virial + md%potential_virial
endif

! trajectory files
write(676,*) path_integral%total_beads + classical%natoms
write(676,*)' i = ', iprint
do irp = 1, path_integral%num_rp
  do j = 1, path_integral%rp(irp)%num_beads
    write(676,'(a2, 3f15.8)') trim(path_integral%rp(irp)%particle_label(j)&
        (len(trim(path_integral%rp(irp)%particle_label(j)))-1:)), &
        path_integral%rp(irp)%x(j,:)/AtoBohr!, path_integral%rp(irp)%traj(i,j,4:6)/(AtoBohr/fstota)
  enddo
enddo
do j = 1, classical%natoms
  write(676,'(a2, 3f15.8)') trim(classical%atom_kinds(j)), classical%x(j,:)/AtoBohr!, classical%traj(i,j,4:6)/(AtoBohr/fstota)
enddo

com = 0.0d0
do irp = 1, path_integral%num_rp
  do j = 1, path_integral%rp(irp)%num_beads
    com(1:3) = com(1:3) + path_integral%rp(irp)%x(j,1:3)/AtoBohr
  enddo
enddo
do j = 1, classical%natoms
  com(1:3) = com(1:3) + classical%x(j,1:3)/AtoBohr
enddo
com = com/dble(real(path_integral%total_beads + classical%natoms))

write(677,*) path_integral%total_beads + classical%natoms
write(677,*)' i = ', iprint
do irp = 1, path_integral%num_rp
  do j = 1, path_integral%rp(irp)%num_beads
    write(677,'(a2, 3f15.8)') trim(path_integral%rp(irp)%particle_label(j)&
        (len(trim(path_integral%rp(irp)%particle_label(j)))-1:)), &
        path_integral%rp(irp)%x(j,:)/AtoBohr -com !, path_integral%rp(irp)%traj(i,j,4:6)/(AtoBohr/fstota)
  enddo
enddo
do j = 1, classical%natoms
  write(677,'(a2, 3f15.8)') trim(classical%atom_kinds(j)), classical%x(j,:)/AtoBohr -com !, classical%traj(i,j,4:6)/(AtoBohr/fstota)
enddo


! ------------------------------------------------- !
! THE MOLECULAR DYNAMICS LOOP FOR SIMULATIONS WITH  !
! PATH-INTEGRALS AND CLASSICAL ATOMS WHERE THE ---- !
! CLASSICAL FORCES ARE CALCULATED USING THE INTERNAL!
! SUBROUTINES IN THE LINKED CP2K LIBRARIES. ------- !
! ------------------------------------------------- !
do istep = 1, md%nsteps
  !write(*,*)' step ', istep
  if (istep == 1 .and. para_var%my_rank == 0) then
    call system_clock(count_rate=irate)
    call system_clock(itime_istart)
  else if (istep == 2 .and. para_var%my_rank == 0) then
    call system_clock(count_rate=trate)
    call system_clock(itime_tstart)
  endif
  ! ------------------------------------------------ !
  ! BEGIN SUBROUTINE BY ZEROING OUT ALL FORCES FOR - !
  ! BOTH PATH-INTEGRAL AND CLASSICAL DOF ----------- !
  ! ------------------------------------------------ !
  do irp = 1, path_integral%num_rp
    path_integral%rp(irp)%force = 0.0d0
  enddo
  classical%force = 0.d0

    ! ------------------------------------------------ !
    ! IF DOING NVT SIMULATIONS, THE FIRST CALL TO THE- !
    ! THERMOSTAT. CURRENT ACCEPTED THERMOSTATS ARE --- !
    ! PILE: PATH-INTEGRAL LANGEVIN EQUATION ---------- !
    ! PILE-G: PATH-INTEGRAL LANGEVIN EQUATION FOR JUST !
    ! THE CENTROID VELOCITIES. ----------------------- !
    ! ------------------------------------------------ !
  if (md%ensemble == 'nvt') then 
    if (istep == 1 .and. para_var%my_rank == 0) then
      call system_clock(count_rate=crate)
      call system_clock(itime_cstart)
    endif
    if (trim(adjustl(md%thermostat)) == 'pile') CALL thermostat_pile(classical, path_integral, md, md%dt, para_var)    !call the initial part of the pile thermostat
    if (trim(adjustl(md%thermostat)) == 'pileg') CALL thermostat_pile_g(classical, path_integral, md, md%dt, para_var)    !call the initial part of the pile thermostat
    if (istep == 1 .and. para_var%my_rank == 0) then
      call system_clock(itime_cend)
      open(111,file='runtime_info.dat',position='append')
      write(111,*)' time for thermostat step ', itime_cend-itime_cstart, real(itime_cend-itime_cstart)/real(crate)
      close(111)
    endif
  endif

     ! ----------------------------------------------- !
     ! VELOCITY-VERLET STEP WITH THE CLASSICAL ATOMS - !
     ! FOR A TIMESTEP OF DT/2 ------------------------ !
     ! ----------------------------------------------- !
  if (classical%natoms > 0) then
    CALL vv_update_classical(classical, path_integral, md, md%dt/2.0d0, istep, para_var)
  endif

  if (trim(md%propagation) == 'normal_fftw' .and. path_integral%num_rp > 0) then
     !  write(*,*)' vv_update_pi BAJ error'
     !  CALL vv_update_pi_fftw(md%dt)     !THIS (CALL TO FFTW) NEEDS TO BE FIXED
    write(*,*)' option not currently supported! '
  else if (trim(md%propagation) == 'normal_manual') then
    ! -----------------------------------------------------------!
    !uses a normal mode propoagtor for the PIs. Calculcates the -!
    !propagation/rotation matrix manually (not using FFTW) ------!
    ! -----------------------------------------------------------!
    classical%force = 0.0d0
    do irp = 1, path_integral%num_rp
      path_integral%rp(irp)%force = 0.0d0
    enddo
    CALL vv_update_pi_normal_manual(classical, path_integral, md, md%dt,istep, para_var)
  else if (trim(md%propagation) == 'none' .and. path_integral%num_rp > 0) then
    ! -----------------------------------------------------------!
    ! uses a standard velocity-verlet integrator to step each ---!
    ! individual bead. ------------------------------------------!
    ! -----------------------------------------------------------!
    classical%force = 0.0d0
    do irp = 1, path_integral%num_rp
      path_integral%rp(irp)%force = 0.0d0
    enddo
    CALL vv_update_pi(classical,path_integral,md, md%dt,istep)
  else if (path_integral%num_rp > 0) then
    write(*,*)' this propagation scheme does not match an option give '
    STOP
  endif

  ! ---------------------------------------------------------!
  ! ZEROS OUT FORCES AND RECALCULCATES FORCES AT NEW UPDATED-!
  ! POSITIONS -----------------------------------------------!
  do irp = 1, path_integral%num_rp            
    path_integral%rp(irp)%force = 0.0d0       
  enddo                                 
  classical%force = 0.0d0

  if (classical%natoms > 0) then
    CALL vv_update_classical(classical, path_integral, md, md%dt/2.d0, istep, para_var)
  endif

    ! --------------------------------------------- !
    ! Second call to the thermostat --------------- !
    ! --------------------------------------------- !
  if (md%ensemble == 'nvt') then 
    if (trim(adjustl(md%thermostat)) == 'pile') CALL thermostat_pile(classical, path_integral, md, md%dt, para_var)    !call the initial part of the pile thermostat
    if (trim(adjustl(md%thermostat)) == 'pileg') CALL thermostat_pile_g(classical, path_integral, md, md%dt, para_var)    !call the initial part of the pile thermostat
  endif
     ! ----------------------------------------------- !
     ! AFTER FINAL CLASSICAL UDPATE FOR EACHS STEP, WE !
     ! NEED TO UPDATE THE FORCE_ENV INSIDE CP2K SO --- !
     ! THAT THE FORCE CALCULATIONS (ALONG WITH ENERGY, !
     ! ETC) WILL ALL BE UPDATED CORRECTLY FOR THE NEXT !
     ! STEP ------------------------------------------ !
     ! EDIT NEEDED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !
     ! ----------------------------------------------- !
     !call update_force_env(classical, path_integral, md)


     ! ----------------------------------------------- !
     ! CALCULCATE THE CONSERVED QUANTITIES (ENERGY FOR !
     ! NVE SIMULATIONS) AND PRINT TO FILE AT EACH STEP !
     ! ----------------------------------------------- !
  CALL update_force_env(classical, path_integral, md)
  do irp = 1, path_integral%num_rp
    path_integral%rp(irp)%force = 0.0d0
  enddo
  classical%force = 0.d0
  CALL calc_force_cp2k(classical, path_integral, md)
  CALL calc_force_pi(path_integral, md, 0)

  CALL calc_total_energy(classical, path_integral, md, istep)
  

  ! --------------------------------------------!
  ! PRINT THE TRAJECTORY AND ENERGIES IF THE ---!
  ! PRINT FREQUENCY IS MET.---------------------!
  if(mod(istep,md%freq_print)==0) THEN
    write(203,'(i5,f16.2,7f20.10)')istep, istep*md%dt/fstota,classical%kinetic,classical%potential, path_integral%kinetic, &
        path_integral%potential, md%kinetic_virial, md%potential_virial, md%kinetic_primitive
    if (classical%natoms+path_integral%total_beads <= 2) then
      write(204,'(i5,2f16.2,4f20.10)')istep, istep*md%dt/fstota,(classical%kinetic+path_integral%kinetic)/&
        ((real(md%space_dim)*(classical%natoms+path_integral%total_beads)-5)*kb/2.0d0), &
        md%kinetic_virial, md%potential_virial, md%kinetic_primitive, md%kinetic_virial+md%potential_virial
    else 
      write(204,'(i5,2f16.2,4f20.10)')istep, istep*md%dt/fstota,(classical%kinetic+path_integral%kinetic)/&
        ((real(md%space_dim)*(classical%natoms+path_integral%total_beads)-(3.0+md%n_constraints))*kb/2.0d0), &
        md%kinetic_virial, md%potential_virial, md%kinetic_primitive, md%kinetic_virial+md%potential_virial
    endif
  endif
  if (istep >= int(md%nsteps*0.4)) then
    if (icount == 0) then
      run_kin = md%kinetic_virial
      run_pot = md%potential_virial
      run_tot = md%kinetic_virial+md%potential_virial
      run_temp = (classical%kinetic+path_integral%kinetic)/((real(md%space_dim)*&
          (classical%natoms+path_integral%total_beads-1))*kb/2.0d0)
      icount = icount + 1
    else
      run_kin = (run_kin*real(icount-1) + md%kinetic_virial )/real(icount)
      run_pot = (run_pot*real(icount-1) + md%potential_virial)/real(icount)
      run_tot = (run_tot*real(icount-1) + (md%kinetic_virial+md%potential_virial))/real(icount)
      run_temp = (run_temp*real(icount-1) + (classical%kinetic+path_integral%kinetic)/((real(md%space_dim)*&
          (classical%natoms+path_integral%total_beads)-(3.0+md%n_constraints))*kb/2.0d0))/real(icount)
     if(mod(istep,md%freq_print)==0) write(205,'(i5,2f16.2,3f20.10)')istep, istep*md%dt/fstota,run_temp, run_kin, run_pot, run_tot
     icount = icount+1
    endif
  endif

     ! ----------------------------------------------- !
     ! BASED UPONNTHE INPUT FREQ_PRINT, WE WILL SAVE- !
     ! THE CLASSICAL AND RP POSITIONS IN A TRAJECTORY- !
     ! VARIABLE IF PRINT_LEVEL IS HIGH --------------- ! 
     ! ----------------------------------------------- !
  if(mod(istep,md%freq_print)==0) THEN
     if (trim(adjustl(md%print_level)) == 'high') then
       iprint = iprint+1
       do irp = 1, path_integral%num_rp
         do ibead = 1, path_integral%rp(irp)%num_beads
           path_integral%rp(irp)%traj(iprint,ibead,1:md%space_dim)=path_integral%rp(irp)%x(ibead,:)
           path_integral%rp(irp)%traj(iprint,ibead,md%space_dim+1:2*md%space_dim)&
               =path_integral%rp(irp)%p(ibead,1:md%space_dim)/path_integral%rp(irp)%kinetic_mass
         enddo
         path_integral%rp(irp)%com_pos(iprint,:) = path_integral%rp(irp)%centroid
       enddo
       do iclassical = 1, classical%natoms
         classical%traj(iprint,iclassical,1:md%space_dim) = classical%x(iclassical,:)
         classical%traj(iprint,iclassical,md%space_dim+1:2*md%space_dim) = classical%p(iclassical,:)/classical%mass(iclassical)
       enddo
     endif

     write(676,*) path_integral%total_beads + classical%natoms
     write(676,*)' i = ', iprint
     do irp = 1, path_integral%num_rp
       do j = 1, path_integral%rp(irp)%num_beads
         write(676,'(a2, 3f15.8)') trim(path_integral%rp(irp)%particle_label(j)&
             (len(trim(path_integral%rp(irp)%particle_label(j)))-1:)), &
             path_integral%rp(irp)%x(j,:)/AtoBohr!, path_integral%rp(irp)%traj(i,j,4:6)/(AtoBohr/fstota)
       enddo
     enddo
     do j = 1, classical%natoms
       write(676,'(a2, 3f15.8)') trim(classical%atom_kinds(j)), classical%x(j,:)/AtoBohr!, classical%traj(i,j,4:6)/(AtoBohr/fstota)
     enddo

     com = 0.0d0
     do irp = 1, path_integral%num_rp
       do j = 1, path_integral%rp(irp)%num_beads
         com(1:3) = com(1:3) + path_integral%rp(irp)%x(j,1:3)/AtoBohr
       enddo
     enddo
     do j = 1, classical%natoms
       com(1:3) = com(1:3) + classical%x(j,1:3)/AtoBohr
     enddo
     com = com/dble(real(path_integral%total_beads + classical%natoms))

     write(677,*) path_integral%total_beads + classical%natoms
     write(677,*)' i = ', iprint
     do irp = 1, path_integral%num_rp
       do j = 1, path_integral%rp(irp)%num_beads
         write(677,'(a2, 3f15.8)') trim(path_integral%rp(irp)%particle_label(j)&
             (len(trim(path_integral%rp(irp)%particle_label(j)))-1:)), &
             path_integral%rp(irp)%x(j,:)/AtoBohr -com !, path_integral%rp(irp)%traj(i,j,4:6)/(AtoBohr/fstota)
       enddo
     enddo
     do j = 1, classical%natoms
       write(677,'(a2, 3f15.8)') trim(classical%atom_kinds(j)), classical%x(j,:)/AtoBohr -com !, classical%traj(i,j,4:6)/(AtoBohr/fstota)
     enddo

     write(678,*) path_integral%num_rp + classical%natoms
     write(678,*)' i = ', iprint
     do irp = 1, path_integral%num_rp
         !com(1) = sum(path_integral%rp(irp)%x(:,1))/path_integral%rp(irp)%num_beads
         !com(2) = sum(path_integral%rp(irp)%x(:,2))/path_integral%rp(irp)%num_beads
         !com(3) = sum(path_integral%rp(irp)%x(:,3))/path_integral%rp(irp)%num_beads
         !write(678,'(a2, 3f15.8)') trim(path_integral%rp(irp)%particle_label(1)&
         !    (len(trim(path_integral%rp(irp)%particle_label(1)))-1:)), &
         !    com(1:3)/AtoBohr !, path_integral%rp(irp)%traj(i,j,4:6)/(AtoBohr/fstota)
         write(678,'(a2, 3f15.8)') trim(path_integral%rp(irp)%particle_label(1)&
             (len(trim(path_integral%rp(irp)%particle_label(1)))-1:)), &
             path_integral%rp(irp)%centroid(1:3)/AtoBohr !, path_integral%rp(irp)%traj(i,j,4:6)/(AtoBohr/fstota)
     enddo
     do j = 1, classical%natoms
       write(678,'(a2, 3f15.8)') trim(classical%atom_kinds(j)), classical%x(j,:)/AtoBohr !, classical%traj(i,j,4:6)/(AtoBohr/fstota)
     enddo
  endif

  if (istep == 1 .and. para_var%my_rank == 0) then
    call system_clock(itime_iend)
    open(111,file='runtime_info.dat',position='append')
      write(111,*)' time for one step ', itime_iend-itime_istart, real(itime_iend-itime_istart)/real(irate)
    close(111)
  else if (istep == 2 .and. para_var%my_rank == 0) then
    call system_clock(itime_tend)
    open(111,file='runtime_info.dat',status="old", position="append")
      write(111,*)' time for one and two step ', itime_tend-itime_tstart, real(itime_tend-itime_tstart)/real(trate)
    close(111)
  endif
enddo
if (trim(adjustl(md%print_level)) == 'high') close(801)
close(675)
close(676)
close(677)
close(678)
close(203)
close(204)
close(205)

open(678, file='final_structure.xyz')
write(678,*) path_integral%total_beads + classical%natoms
write(678,*)' i = ', iprint
do irp = 1, path_integral%num_rp
  do j = 1, path_integral%rp(irp)%num_beads
    write(678,'(a2, 3f15.8)') trim(path_integral%rp(irp)%particle_label(j)&
             (len(trim(path_integral%rp(irp)%particle_label(j)))-1:)), &
             path_integral%rp(irp)%x(j,:)/AtoBohr !, path_integral%rp(irp)%traj(i,j,4:6)/(AtoBohr/fstota)
  enddo
enddo
do j = 1, classical%natoms
  write(678,'(a2, 3f15.8)') trim(classical%atom_kinds(j)), classical%x(j,:)/AtoBohr !, classical%traj(i,j,4:6)/(AtoBohr/fstota)
enddo
close(678)
open(678, file='final_rp-vel.xyz')
do irp = 1, path_integral%num_rp
  do j = 1, path_integral%rp(irp)%num_beads
    write(678,'(3f15.8)') path_integral%rp(irp)%p(j,:)/path_integral%rp(irp)%kinetic_mass !, path_integral%rp(irp)%traj(i,j,4:6)/(AtoBohr/fstota)
  enddo
enddo
close(678)
open(678, file='final_atom-vel.xyz')
do j = 1, classical%natoms
  write(678,'(3f15.8)') classical%p(j,:)/classical%mass(j) !, classical%traj(i,j,4:6)/(AtoBohr/fstota)
enddo
close(678)

if (trim(adjustl(md%print_level)) == 'high') then
!writes the trajectory file
  do irp = 1, path_integral%num_rp
    open(201,file='trajectory_pimd_'//char(irp+48)//'.xyz')
    write(*,*)' print the md traj file '
    do i = 0, md%nsteps/md%freq_print
      write(201,*) path_integral%rp(irp)%num_beads
      write(201,*)' i = ', i
      do j = 1, path_integral%rp(irp)%num_beads
        write(201,'(a2, 6f15.8)') trim(path_integral%rp(irp)%particle_label(j)), path_integral%rp(irp)%traj(i,j,1:3)/AtoBohr, &
            path_integral%rp(irp)%traj(i,j,4:6)/(AtoBohr)
      enddo
    enddo
    close(201) 
  enddo
endif


END SUBROUTINE run_md
! ----------------------------------------------------!


SUBROUTINE update_force_env(classical, path_integral, md)
use libcp2k
implicit none

TYPE(path_integral_env)                ::path_integral
TYPE(classical_type)                   ::classical
TYPE(md_info)                          ::md
REAL(KIND=8),DIMENSION(md%space_dim*(classical%natoms+path_integral%total_beads))     ::positions
REAL(KIND=8),DIMENSION(md%space_dim*(classical%natoms+path_integral%total_beads))     ::velocities
integer                              ::i, irp, icount, ierror

include 'mpif.h'

icount = 0
do irp = 1, path_integral%num_rp
  do i = 1, path_integral%rp(irp)%num_beads
    icount = icount+1
    positions( (icount-1)*md%space_dim+1: icount*md%space_dim) = path_integral%rp(irp)%x(i,:)
    velocities( (icount-1)*md%space_dim+1: icount*md%space_dim) = path_integral%rp(irp)%p(i,:)/path_integral%rp(irp)%kinetic_mass
  enddo
enddo
if (icount /= path_integral%total_beads) WRITE(*,*)' ERROR IN FORCE UPDATE ENV PI INDEX '
do i = 1, classical%natoms
  positions( (i-1+path_integral%total_beads)*md%space_dim+1: (i+path_integral%total_beads)*md%space_dim ) = &
    classical%x(i,1:md%space_dim) 
  velocities( (i-1+path_integral%total_beads)*md%space_dim+1: (i+path_integral%total_beads)*md%space_dim ) = &
    classical%p(i,1:md%space_dim)/classical%mass(i)
enddo

call MPI_BARRIER(MPI_COMM_WORLD,ierror)
CALL cp2k_set_positions(classical%cp2k_env_flag,positions,md%space_dim*(classical%natoms+path_integral%total_beads)) 
call MPI_BARRIER(MPI_COMM_WORLD,ierror)
!CALL cp2k_set_velocities(classical%cp2k_env_flag,velocities,md%space_dim*(classical%natoms+path_integral%total_beads))

END SUBROUTINE update_force_env


! ----------------------------------------- !
SUBROUTINE calc_total_energy(classical, path_integral, md, istep)
implicit none


TYPE(path_integral_env)                ::path_integral
TYPE(classical_type)                   ::classical
TYPE(md_info)                          ::md
real(kind=8)       ::cp2k_pot, cp2k_kin
real(kind=8)       ::pimd_pot, pimd_kin
real(kind=8)       ::tot_energy, kinetic_virial
integer            ::istep, ierror


CALL calc_cp2k_energy(classical, cp2k_pot, cp2k_kin)

CALL calc_pimd_energy(path_integral, pimd_pot, pimd_kin)


if (md%xx_present) then
  call calc_xx_energy(classical, path_integral, md, cp2k_pot)
endif

classical%potential = cp2k_pot
classical%kinetic = cp2k_kin
path_integral%potential = pimd_pot
path_integral%kinetic = pimd_kin

CALL primitive_total_energy_estimator(path_integral, classical, cp2k_pot, tot_energy, md)
CALL kinetic_virial_estimator(path_integral, classical, md, kinetic_virial)
md%kinetic_primitive = tot_energy
md%kinetic_virial = kinetic_virial
md%potential_virial = cp2k_pot

END SUBROUTINE calc_total_energy


! ------------------------------------------------- !
SUBROUTINE calc_cp2k_energy(classical, pot_energy, kin_energy)
use libcp2k
implicit none

TYPE(classical_type)  ::classical
real(kind=8)          ::pot_energy, kin_energy
integer               ::i, ierror

include 'mpif.h'

call MPI_BARRIER(MPI_COMM_WORLD,ierror)
CALL cp2k_calc_energy_force(classical%cp2k_env_flag)
call MPI_BARRIER(MPI_COMM_WORLD,ierror)

CALL cp2k_get_potential_energy(classical%cp2k_env_flag, pot_energy)
!write(*,*)'classical potential energy is ', pot_energy
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



END SUBROUTINE calc_cp2k_energy

! ----------------------------------------- !

SUBROUTINE kinetic_virial_estimator(path_integral, classical, md, kinetic_virial)
implicit none

TYPE(path_integral_env)                ::path_integral
TYPE(classical_type)                   ::classical
TYPE(md_info)                          ::md
real(kind=8)                           ::kinetic_virial, extra_term, bead_term, bead_term_local
integer                                ::nbead, irp, i
real(kind=8),dimension(md%space_dim)   ::centroid

extra_term = 0.0d0
bead_term = 0.0d0
do irp = 1, path_integral%num_rp
  bead_term_local = 0.0d0
  nbead = path_integral%rp(irp)%num_beads
  !write(*,*)' nbead is ', nbead
  do i = 1, md%space_dim
    centroid(i) = sum(path_integral%rp(irp)%x(:,i))&
    /dble(real(path_integral%rp(irp)%num_beads))
  enddo
  !write(*,*)' centroid is ', centroid(:)
  do i = 1, nbead
	  bead_term_local = bead_term_local + dot_product(path_integral%rp(irp)%x(i,:)-centroid(:),path_integral%rp(irp)%force(i,:))
	  !write(*,*)' dot product is ', dot_product(path_integral%rp(irp)%x(i,:)-centroid(:),path_integral%rp(irp)%force(i,:))
  enddo
  !write(*,*)' bead_term_local is ', bead_term_local
  bead_term = bead_term + bead_term_local/(2.0d0*dble(real(nbead)))
enddo

!write(*,*)' bead_term is ', bead_term
  
if (path_integral%num_rp > 0) then 
  extra_term = md%space_dim*dble(real(classical%natoms - path_integral%num_rp))/(2.0d0*path_integral%rp(1)%beta)
else
  extra_term = md%space_dim*dble(real(classical%natoms))/(2.0d0*classical%beta)
endif

!write(*,*)' extra_term is ', extra_term

kinetic_virial = bead_term + extra_term

END SUBROUTINE kinetic_virial_estimator

! ----------------------------------------- !

SUBROUTINE primitive_total_energy_estimator(path_integral, classical, classical_pot, total_energy, md)
implicit none

TYPE(path_integral_env)  ::path_integral
TYPE(classical_type)                   ::classical
TYPE(md_info)           ::md
real(kind=8)          ::pot_energy, omega, extra_term
real(kind=8)          ::classical_pot, total_energy
integer               ::irp,i,iplus,nbead

pot_energy = 0.0d0
extra_term = 0.0d0
do irp = 1, path_integral%num_rp
  omega = ((dble(real(path_integral%rp(irp)%num_beads)))/path_integral%rp(irp)%beta)&
        *dsqrt(path_integral%rp(irp)%harm_mass/path_integral%rp(irp)%kinetic_mass)
  nbead = path_integral%rp(irp)%num_beads
  do i = 1, nbead
    iplus = i+1
    if (i==nbead) iplus = 1
    pot_energy = pot_energy - 0.5d0*path_integral%rp(irp)%kinetic_mass*omega**2 * &      !BAJ: Feb8 made change
      dot_product(path_integral%rp(irp)%x(i,:)-path_integral%rp(irp)%x(iplus,:),&
      path_integral%rp(irp)%x(i,:)-path_integral%rp(irp)%x(iplus,:))
  enddo
enddo

if (path_integral%num_rp > 0) then 
  extra_term = dble(3.0d0*real(path_integral%total_beads + classical%natoms-md%n_constraints))/(2.0d0*path_integral%rp(1)%beta)
else
  extra_term = dble(3.0d0*real(classical%natoms-md%n_constraints))/(2.0d0*classical%beta)
endif

total_energy = pot_energy + extra_term     !BAJ: March27 made change above to include
                                                           !constraint and
                                                           !space_dim
                                                           !information. 
!total_energy = classical_pot - pot_energy + extra_term


END SUBROUTINE primitive_total_energy_estimator


! ------------------------------------------------- !

SUBROUTINE calc_pimd_energy(path_integral, pot_energy,kin_energy)
implicit none

TYPE(path_integral_env)  ::path_integral
real(kind=8)          ::pot_energy, kin_energy, omega
integer               ::irp,i,iplus,nbead

pot_energy = 0.0d0
kin_energy = 0.0d0
do irp = 1, path_integral%num_rp
  !omega = 1.0d0/path_integral%rp(irp)%beta
  !omega = dsqrt(dble(real(path_integral%rp(irp)%num_beads)))/path_integral%rp(irp)%beta
  omega = ((dble(real(path_integral%rp(irp)%num_beads)))/path_integral%rp(irp)%beta)&
        *dsqrt(path_integral%rp(irp)%harm_mass/path_integral%rp(irp)%kinetic_mass)
  nbead = path_integral%rp(irp)%num_beads
  do i = 1, nbead
    iplus = i+1
    if (i==nbead) iplus = 1
    pot_energy = pot_energy + 0.5d0*path_integral%rp(irp)%kinetic_mass*omega**2 * &      !BAJ: Feb8 made change
      dot_product(path_integral%rp(irp)%x(i,:)-path_integral%rp(irp)%x(iplus,:),&
      path_integral%rp(irp)%x(i,:)-path_integral%rp(irp)%x(iplus,:))
    !write(*,*)' prefactor is ', 0.5d0*path_integral%rp(irp)%harm_mass*omega**2/dble(path_integral%rp(irp)%num_beads) 
    !write(*,*)' dot product is ', dot_product(path_integral%rp(irp)%x(i,:)-path_integral%rp(irp)%x(iplus,:),&
    !  path_integral%rp(irp)%x(i,:)-path_integral%rp(irp)%x(iplus,:))
    kin_energy = kin_energy + 0.5d0*dot_product(path_integral%rp(irp)%p(i,:),path_integral%rp(irp)%p(i,:))/&
      (path_integral%rp(irp)%kinetic_mass)
  enddo
enddo

END SUBROUTINE calc_pimd_energy


! ---------------------------------------------- !
! The velocity-verlet interrator for the ------- !
! classical atoms with a time-step of dt ------- !
! The following updates occur:  ---------------- !
! p(t + dt/2) <- p(t) + dt/2 * F(t) ------------ !
! x(t+dt) <- x(t) + dt * p(t+dt/2)/mass -------- !
! p(t+dt) <- p(t + dt/2) + dt/2 * F(t+dt/2) ---- !
! ---------------------------------------------- !
SUBROUTINE vv_update_pi(classical, path_integral, md, dt, step)
implicit none

TYPE(classical_type)    ::classical
TYPE(path_integral_env) ::path_integral
TYPE(md_info)           ::md
real(kind=8)            ::dt
integer                 ::irp,i, j, step

!classical%force = 0.0d0
!do irp = 1, path_integral%num_rp
!  path_integral%rp(irp)%force = 0.0d0
!enddo

CALL update_force_env(classical, path_integral, md)
CALL calc_force_cp2k(classical, path_integral, md)
CALL calc_force_pi(path_integral, md, 0)

CALL update_pi_p(path_integral, dt)
CALL update_pi_x(path_integral, md, dt)


classical%force = 0.0d0
do irp = 1, path_integral%num_rp
  path_integral%rp(irp)%force = 0.0d0
enddo
CALL update_force_env(classical, path_integral, md)
CALL calc_force_cp2k(classical, path_integral, md)
CALL calc_force_pi(path_integral, md, 0)


CALL update_pi_p(path_integral, dt)

END SUBROUTINE vv_update_pi

! ---------------------------------------------- !
! The velocity-verlet interrator for the ------- !
! classical atoms with a time-step of dt ------- !
! The following updates occur:  ---------------- !
! p(t + dt/2) <- p(t) + dt/2 * F(t) ------------ !
! x(t+dt) <- x(t) + dt * p(t+dt/2)/mass -------- !
! p(t+dt) <- p(t + dt/2) + dt/2 * F(t+dt/2) ---- !
! ---------------------------------------------- !
SUBROUTINE vv_update_classical(classical, path_integral, md, dt, step, para_var)
implicit none

TYPE(classical_type)    ::classical
TYPE(path_integral_env) ::path_integral
TYPE(md_info)           ::md
TYPE(parallel_env)                     ::para_var
real(kind=8)            ::dt
integer                 ::irp,i, j, step
integer         ::crate, itime_cstart, itime_cend

!classical%force = 0.0d0
!do irp = 1, path_integral%num_rp
!  path_integral%rp(irp)%force = 0.0d0
!enddo

    if (step == 1 .and. para_var%my_rank == 0) then
      call system_clock(count_rate=crate)
      call system_clock(itime_cstart)
    endif
CALL update_force_env(classical, path_integral, md)
    if (step == 1 .and. para_var%my_rank == 0) then
      call system_clock(itime_cend)
      open(111,file='runtime_info.dat',position='append')
      write(111,*)' time for update environment step ', itime_cend-itime_cstart, real(itime_cend-itime_cstart)/real(crate)
      close(111)
    endif
    if (step == 1 .and. para_var%my_rank == 0) then
      call system_clock(count_rate=crate)
      call system_clock(itime_cstart)
    endif
CALL calc_force_cp2k(classical, path_integral, md)
    if (step == 1 .and. para_var%my_rank == 0) then
      call system_clock(itime_cend)
      open(111,file='runtime_info.dat',position='append')
      write(111,*)' time for calc force cp2k step ', itime_cend-itime_cstart, real(itime_cend-itime_cstart)/real(crate)
      close(111)
    endif
CALL calc_force_pi(path_integral, md, 0)

if (trim(adjustl(md%print_level)) == 'high') then
  write(801,*) sum(path_integral%rp(:)%num_beads) + classical%natoms
  write(801,*)' dt = ', dt*(step-1)
  do irp = 1, path_integral%num_rp
    do j = 1, path_integral%rp(irp)%num_beads
      write(801,'(a2, 6f15.8)') trim(path_integral%rp(irp)%particle_label(j)), path_integral%rp(irp)%force(j,:)
    enddo
  enddo
  do j = 1, classical%natoms
    write(801,'(a2, 6f15.8)') trim(classical%atom_kinds(j)), classical%force(j,:)
  enddo
endif

CALL update_classical_p(classical, dt)
CALL update_classical_x(classical, dt)


classical%force = 0.0d0
do irp = 1, path_integral%num_rp
  path_integral%rp(irp)%force = 0.0d0
enddo
CALL update_force_env(classical, path_integral, md)
CALL calc_force_cp2k(classical, path_integral, md)
CALL calc_force_pi(path_integral, md, 0)


CALL update_classical_p(classical, dt)

END SUBROUTINE vv_update_classical
! ---------------------------------------------- !


! ---------------------------------------------- !
! Subroutine that evaluates the forces on the -- !
! classical atoms using the CP2K subroutines --- !
! ---------------------------------------------- !
SUBROUTINE calc_force_cp2k(classical, path_integral, md)
use libcp2k
implicit none

TYPE(classical_type)                     ::classical
TYPE(path_integral_env)                  ::path_integral
TYPE(md_info)                            ::md
real(kind=8),allocatable,dimension(:)    ::force
integer                                  ::i, irp, icount, ierror

include 'mpif.h'

ALLOCATE(force((classical%natoms+path_integral%total_beads)*md%space_dim))
force = 0.0d0

! ------------------------------------------------ !
! CALL CP2K INTERNAL FUNCTION CP2K_GET_FORCES ---- !
! THAT CALCULATES THE FORCES BETWEEN CLASSICAL --- !
! ATOMS ARE RETURNS THEM IN AN ARRAY OF SIZE ----- !
! NATOMS*DIMENSION ------------------------------- !
! ------------------------------------------------ !
call MPI_BARRIER(MPI_COMM_WORLD,ierror)
CALL cp2k_calc_energy_force(classical%cp2k_env_flag)
call MPI_BARRIER(MPI_COMM_WORLD,ierror)
CALL cp2k_get_forces(classical%cp2k_env_flag, force, (classical%natoms+path_integral%total_beads)*md%space_dim)

! ------------------------------------------------ !
! CONVERTS THE ARRAY FORCES(1:NATOMS*DIMENSION) -- !
! INTO ARRAY CLASSICAL%FORCE(1:NATOMS,1:DIMENSION) !
! WHCIH IS JUST MORE INTUITIVE. ------------------ !
! IMPORTANT NOTE: IF CLASSICAL$FORCE ARRAY HAS NOT !
! BEEN ZERO'D BEFORE CALLING SUBROUTINE, IT WILL - !
! BE ADDED TO THE EXISTING FORCES THERE. --------- !
! ------------------------------------------------ !
icount = 0
do irp = 1, path_integral%num_rp
  do i = 1, path_integral%rp(irp)%num_beads
    icount = icount+1
    path_integral%rp(irp)%force(i,:) = path_integral%rp(irp)%force(i,:) + force( (icount-1)*md%space_dim+1: icount*md%space_dim)
  enddo
enddo
if (icount /= path_integral%total_beads) WRITE(*,*)' ERROR IN FORCE CALCULATION ENV PI INDEX '
do i = 1, classical%natoms
  classical%force(i,1:md%space_dim) = classical%force(i,1:md%space_dim) + force( (i-1+path_integral%total_beads)&
     *md%space_dim+1:(i+path_integral%total_beads)*md%space_dim )
enddo

DEALLOCATE(force)

if (md%xx_present) then
  call calc_xx_exclusions(classical, path_integral, md)
endif

END SUBROUTINE calc_force_cp2k
! --------------------------------------------- !


SUBROUTINE calc_xx_exclusions(classical, path_integral, md)
implicit none

TYPE(classical_type)                     ::classical
TYPE(path_integral_env)                  ::path_integral
TYPE(md_info)                            ::md
integer                                  ::irp, ixx, ibead, iclassical
real(kind=8)                             ::dist
real(kind=8),dimension(md%space_dim)     ::dist_vec, force


!open(999,file='xx_force_update.dat',position='append')
do irp = 1, path_integral%num_rp
  if (trim(path_integral%rp(irp)%particle_label(1) &
    (len(trim(path_integral%rp(irp)%particle_label(1)))-1:)) /= 'CO') then
  if (path_integral%rp(irp)%charge > 0.0d0) then
    do ixx = 1, path_integral%rp(md%xx_atom_num)%num_beads
      dist_vec(:) = path_integral%rp(md%xx_atom_num)%x(ixx,:) - path_integral%rp(irp)%centroid(:)
      dist = dsqrt(dot_product(dist_vec, dist_vec))
      if (dist < md%xx_rcut) then
        !write(999,*)' update forces for PI', irp
        do ibead = 1, path_integral%rp(irp)%num_beads
          dist_vec(:) = path_integral%rp(md%xx_atom_num)%x(ixx,:) - path_integral%rp(irp)%x(ibead,:)
          dist = dsqrt(dot_product(dist_vec, dist_vec))
          force(:) = ((path_integral%rp(md%xx_atom_num)%charge * path_integral%rp(irp)%charge) &
            /dist**2) * (dist_vec(:)/dist)
          !write(999,*)' update force is ', force(:)
          !write(999,*)' old forces are '
          !write(999,*) path_integral%rp(md%xx_atom_num)%force(ixx,:)
          !write(999,*) path_integral%rp(irp)%force(ibead,:)
          path_integral%rp(md%xx_atom_num)%force(ixx,:) = path_integral%rp(md%xx_atom_num)%force(ixx,:) + force
          path_integral%rp(irp)%force(ibead,:) = path_integral%rp(irp)%force(ibead,:) - force
          !write(999,*)' new forces are '
          !write(999,*) path_integral%rp(md%xx_atom_num)%force(ixx,:)
          !write(999,*) path_integral%rp(irp)%force(ibead,:)
        enddo
      endif
    enddo
  endif
  endif
enddo
do iclassical = 1, classical%natoms
  if (trim(classical%atom_kinds(iclassical)) /= 'CO') then
  if (classical%ff_atom_charges(iclassical) > 0.0d0) then
    do ixx = 1, path_integral%rp(md%xx_atom_num)%num_beads
      dist_vec(:) = path_integral%rp(md%xx_atom_num)%x(ixx,:) - classical%x(iclassical,:)
      dist = dsqrt(dot_product(dist_vec, dist_vec))
      if (dist < md%xx_rcut) then
        !write(999,*)' update forces for classical ', iclassical
          force(:) = ((path_integral%rp(md%xx_atom_num)%charge * classical%ff_atom_charges(iclassical)) &
            /dist**2) * (dist_vec(:)/dist)
          !write(999,*)' update force is ', force(:)
          !write(999,*)' old forces are '
          !write(999,*) path_integral%rp(md%xx_atom_num)%force(ixx,:)
          !write(999,*) classical%force(iclassical,:)
          path_integral%rp(md%xx_atom_num)%force(ixx,:) = path_integral%rp(md%xx_atom_num)%force(ixx,:) + force
          classical%force(iclassical,:) = classical%force(iclassical,:) - force
          !write(999,*)' new forces are '
          !write(999,*) path_integral%rp(md%xx_atom_num)%force(ixx,:)
          !write(999,*) classical%force(iclassical,:)
      endif
    enddo
  endif
  endif
enddo
      
close(999) 




END SUBROUTINE calc_xx_exclusions

SUBROUTINE calc_xx_energy(classical, path_integral, md, tot_pot)
implicit none

TYPE(classical_type)                     ::classical
TYPE(path_integral_env)                  ::path_integral
TYPE(md_info)                            ::md
integer                                  ::irp, ixx, ibead, iclassical
real(kind=8)                             ::dist
real(kind=8)                             ::tot_pot, pot
real(kind=8),dimension(md%space_dim)     ::dist_vec


do irp = 1, path_integral%num_rp
  if (trim(path_integral%rp(irp)%particle_label(1) &
    (len(trim(path_integral%rp(irp)%particle_label(1)))-1:)) /= 'CO') then
  if (path_integral%rp(irp)%charge > 0.0d0) then
    do ixx = 1, path_integral%rp(md%xx_atom_num)%num_beads
      dist_vec(:) = path_integral%rp(md%xx_atom_num)%x(ixx,:) - path_integral%rp(irp)%centroid(:)
      dist = dsqrt(dot_product(dist_vec, dist_vec))
      if (dist < md%xx_rcut) then
        do ibead = 1, path_integral%rp(irp)%num_beads
          dist_vec(:) = path_integral%rp(md%xx_atom_num)%x(ixx,:) - path_integral%rp(irp)%x(ibead,:)
          dist = dsqrt(dot_product(dist_vec, dist_vec))
          pot = - ((path_integral%rp(md%xx_atom_num)%charge * path_integral%rp(irp)%charge) &
            /dist**2) + ((path_integral%rp(md%xx_atom_num)%charge * path_integral%rp(irp)%charge) &
            /md%xx_rcut**2) 
          tot_pot = tot_pot + pot 
        enddo
      endif
    enddo
  endif
  endif
enddo
do iclassical = 1, classical%natoms
  if (trim(classical%atom_kinds(iclassical)) /= 'CO') then
  if (classical%ff_atom_charges(iclassical) > 0.0d0) then
    do ixx = 1, path_integral%rp(md%xx_atom_num)%num_beads
      dist_vec(:) = path_integral%rp(md%xx_atom_num)%x(ixx,:) - classical%x(iclassical,:)
      dist = dsqrt(dot_product(dist_vec, dist_vec))
      if (dist < md%xx_rcut) then
          pot = - ((path_integral%rp(md%xx_atom_num)%charge * classical%ff_atom_charges(iclassical)) &
            /dist**2) + ((path_integral%rp(md%xx_atom_num)%charge * classical%ff_atom_charges(iclassical)) &
            /md%xx_rcut**2) 
          tot_pot = tot_pot + pot 
      endif
    enddo
  endif
  endif
enddo
      



END SUBROUTINE calc_xx_energy


SUBROUTINE calc_force_pi(path_integral, md, flag)
implicit none

TYPE(path_integral_env)    ::path_integral
TYPE(md_info)              ::md
integer                    ::flag
real(kind=8)               ::com
integer                    ::icount

if (trim(md%propagation) == 'none' .or. flag /= 0) then
  CALL calc_force_pi_internal(path_integral, md)
endif
if (md%centroid) then
  icount = 1
  do while (md%centroid_list(icount,1) <= path_integral%num_rp)
    if (md%centroid_list(icount,2) == 1) then
      com = sum(path_integral%rp(md%centroid_list(icount,1))%force(:,1))/&
        dble(real(path_integral%rp(md%centroid_list(icount,1))%num_beads))
      path_integral%rp(md%centroid_list(icount,1))%force(:,1) = path_integral%rp(md%centroid_list(icount,1))%force(:,1) - com
    endif
    if (md%centroid_list(icount,3) == 1) then
      com = sum(path_integral%rp(md%centroid_list(icount,1))%force(:,2))/&
        dble(real(path_integral%rp(md%centroid_list(icount,1))%num_beads))
      path_integral%rp(md%centroid_list(icount,1))%force(:,2) = path_integral%rp(md%centroid_list(icount,1))%force(:,2) - com
    endif
    if (md%centroid_list(icount,4) == 1) then
      com = sum(path_integral%rp(md%centroid_list(icount,1))%force(:,3))/&
        dble(real(path_integral%rp(md%centroid_list(icount,1))%num_beads))
      path_integral%rp(md%centroid_list(icount,1))%force(:,3) = path_integral%rp(md%centroid_list(icount,1))%force(:,3) - com
    endif
    icount = icount+1
    if (icount > md%n_centroids) EXIT
  enddo
endif
  

END SUBROUTINE calc_force_pi

! ---------------------------------------------- !
!Subroutine that calculcates the internal forces !
!on the ring-polymer beads due to the bead-bead- !
!spring terms. --------------------------------- !
! ---------------------------------------------- !
SUBROUTINE calc_force_pi_internal(path_integral,md)     !THE BEAD-BEAD PREFACTOR IS NOT RIGHT!!!
implicit none

TYPE(path_integral_env)    ::path_integral
TYPE(md_info)              ::md
integer           ::irp,ibead
integer           ::ibead_minus,ibead_plus
real(kind=8)      ::omega

do irp = 1, path_integral%num_rp
  omega = ((dble(real(path_integral%rp(irp)%num_beads)))/path_integral%rp(irp)%beta)&
        *dsqrt(path_integral%rp(irp)%harm_mass/path_integral%rp(irp)%kinetic_mass)
  do ibead = 1, path_integral%rp(irp)%num_beads
    ibead_minus = ibead-1
    if (ibead == 1) ibead_minus = path_integral%rp(irp)%num_beads
    ibead_plus = ibead+1
    if(ibead == path_integral%rp(irp)%num_beads) ibead_plus = 1
    path_integral%rp(irp)%force(ibead,:) = path_integral%rp(irp)%force(ibead,:) &    !BAJ: Feb8 Made these changes
      - ((path_integral%rp(irp)%harm_mass * path_integral%rp(irp)%num_beads**2) &    !BAJ: Made changes
      /path_integral%rp(irp)%beta**2) &     !CHECK!!!
      *(2.0d0*path_integral%rp(irp)%x(ibead,:) - path_integral%rp(irp)%x(ibead_minus,:) &
      - path_integral%rp(irp)%x(ibead_plus,:) )
    !path_integral%rp(irp)%force(ibead,:) = path_integral%rp(irp)%force(ibead,:) &    !BAJ: Feb8 Made these changes
    !  - path_integral%rp(irp)%kinetic_mass * omega**2&     !CHECK!!!
    !  *(2.0d0*path_integral%rp(irp)%x(ibead,:) - path_integral%rp(irp)%x(ibead_minus,:) &
    !  - path_integral%rp(irp)%x(ibead_plus,:) )
  enddo
enddo

END SUBROUTINE calc_force_pi_internal
! -------------------------------------------- !

! ---------------------------------------------- !
! Subroutine that updates the classical atom --- !
! positions by dt ------------------------------ !
! ---------------------------------------------- !
SUBROUTINE update_pi_x(path_integral, md, dt)
implicit none
TYPE(path_integral_env)    :: path_integral
TYPE(md_info)             ::md
real(kind=8)       ::dt
integer            ::irp, i

do irp = 1, path_integral%num_rp
  path_integral%rp(irp)%x(:,:) = path_integral%rp(irp)%x(:,:) + dt * &
      path_integral%rp(irp)%p(:,:)/path_integral%rp(irp)%kinetic_mass
enddo
do irp = 1, path_integral%num_rp
  do i = 1, md%space_dim
    path_integral%rp(irp)%centroid(i) = sum(path_integral%rp(irp)%x(:,i))&
    /dble(real(path_integral%rp(irp)%num_beads))
  enddo
enddo

END SUBROUTINE update_pi_x


! ---------------------------------------------- !
! Subroutine that updates the classical atom --- !
! momentum by dt/2 ----------------------------- !
! ---------------------------------------------- !
SUBROUTINE update_classical_p(classical, dt)
implicit none

TYPE(classical_type)        ::classical
real(kind=8)                ::dt

classical%p(:,:) = classical%p(:,:) + dt * classical%force(:,:)/2.0d0

END SUBROUTINE update_classical_p



! ---------------------------------------------- !
! Subroutine that updates the classical atom --- !
! positions by dt ------------------------------ !
! ---------------------------------------------- !
SUBROUTINE update_classical_x(classical, dt)
implicit none
TYPE(classical_type)   ::classical
real(kind=8)       ::dt
integer            ::iatom

do iatom = 1, classical%natoms
  classical%x(iatom,:) = classical%x(iatom,:) + dt * classical%p(iatom,:)/( classical%mass(iatom) )
enddo

END SUBROUTINE update_classical_x

! ---------------------------------------------- !
! Subroutine that updates the classical atom --- !
! positions by dt ------------------------------ !
! ---------------------------------------------- !
SUBROUTINE update_pimd_x(path_integral, dt)
implicit none
TYPE(path_integral_env)   ::path_integral
real(kind=8)       ::dt
integer            ::irp, ibead

do irp = 1, path_integral%num_rp
  do ibead = 1, path_integral%rp(irp)%num_beads
    path_integral%rp(irp)%x(ibead,:) = path_integral%rp(irp)%x(ibead,:) + dt * path_integral%rp(irp)%p(ibead,:)&
      /(path_integral%rp(irp)%kinetic_mass)
  enddo
enddo

END SUBROUTINE update_pimd_x

! ----------------------------------------------- !
! The velocity-verlet algorithm for the path ---- !
! integral degrees of freedom. This subroutine - !
! manually calculates the  matrix for the --------- !
! transformation from the cartesian to normal --- !
! coordinates. ---------------------------------- !
! ----------------------------------------------- !
SUBROUTINE vv_update_pi_normal_manual(classical, path_integral, md, dt,icount, para_var)
implicit none

TYPE(classical_type)      ::classical
TYPE(path_integral_env)   ::path_integral
TYPE(md_info)             ::md
TYPE(parallel_env)                     ::para_var
real(kind=8)          :: dt
integer               :: irp,icount,i, j
integer         ::crate, itime_cstart, itime_cend

    if (icount == 1 .and. para_var%my_rank == 0) then
      call system_clock(count_rate=crate)
      call system_clock(itime_cstart)
    endif
CALL update_force_env(classical, path_integral, md)
    if (icount == 1 .and. para_var%my_rank == 0) then
      call system_clock(itime_cend)
      open(111,file='runtime_info.dat',position='append')
      write(111,*)' time for update force env step ', itime_cend-itime_cstart, real(itime_cend-itime_cstart)/real(crate)
      close(111)
    endif
    if (icount == 1 .and. para_var%my_rank == 0) then
      call system_clock(count_rate=crate)
      call system_clock(itime_cstart)
    endif
CALL calc_force_cp2k(classical, path_integral, md)
    if (icount == 1 .and. para_var%my_rank == 0) then
      call system_clock(itime_cend)
      open(111,file='runtime_info.dat',position='append')
      write(111,*)' time for calc_force_cp2k step ', itime_cend-itime_cstart, real(itime_cend-itime_cstart)/real(crate)
      close(111)
    endif
    if (icount == 1 .and. para_var%my_rank == 0) then
      call system_clock(count_rate=crate)
      call system_clock(itime_cstart)
    endif
CALL calc_force_pi(path_integral, md, 0)
    if (icount == 1 .and. para_var%my_rank == 0) then
      call system_clock(itime_cend)
      open(111,file='runtime_info.dat',position='append')
      write(111,*)' time for calc_force_pi step ', itime_cend-itime_cstart, real(itime_cend-itime_cstart)/real(crate)
      close(111)
    endif

!write(888,*) sum(path_integral%rp(:)%num_beads) + classical%natoms
!write(888,*)' dt = ', dt*(icount-1)
!do irp = 1, path_integral%num_rp
!  do j = 1, path_integral%rp(irp)%num_beads
!    write(888,'(a2, 6f15.8)') trim(path_integral%rp(irp)%particle_label(j)), path_integral%rp(irp)%force(j,:)
!  enddo
!enddo
!do j = 1, classical%natoms
!  write(888,'(a2, 6f15.8)') trim(classical%atom_kinds(j)), classical%force(j,:)
!enddo

    if (icount == 1 .and. para_var%my_rank == 0) then
      call system_clock(count_rate=crate)
      call system_clock(itime_cstart)
    endif
CALL update_pi_p(path_integral, dt)

CALL manual_normal_rotation(path_integral)

!write(*,*)' transform cart to normal '
CALL manual_cart_to_normal(path_integral, md, icount)

!write(*,*)' update'
CALL update_normal_coord(path_integral, md, dt)
!CALL update_pimd_x(path_integral, dt)

!write(*,*)' transform normal to cart '
CALL manual_normal_to_cart(path_integral, md, icount)
    if (icount == 1 .and. para_var%my_rank == 0) then
      call system_clock(itime_cend)
      open(111,file='runtime_info.dat',position='append')
      write(111,*)' time for propagate normal mode step ', itime_cend-itime_cstart, real(itime_cend-itime_cstart)/real(crate)
      close(111)
    endif

do irp = 1, path_integral%num_rp
  do i = 1, md%space_dim
    path_integral%rp(irp)%centroid(i) = sum(path_integral%rp(irp)%x(:,i))&
    /dble(real(path_integral%rp(irp)%num_beads))
  enddo
enddo

classical%force = 0.0d0
do irp = 1, path_integral%num_rp
  path_integral%rp(irp)%force = 0.0d0
enddo
!write(*,*)' calc_force_mix BAJ error'
CALL update_force_env(classical, path_integral, md)
CALL calc_force_cp2k(classical, path_integral, md)
CALL calc_force_pi(path_integral, md, 0)

CALL update_pi_p(path_integral, dt)

!write(888,*) sum(path_integral%rp(:)%num_beads) + classical%natoms
!write(888,*)' dt = ', dt*(icount-1)
!do irp = 1, path_integral%num_rp
!  do j = 1, path_integral%rp(irp)%num_beads
!    write(888,'(a2, 6f15.8)') trim(path_integral%rp(irp)%particle_label(j)), path_integral%rp(irp)%force(j,:)
!  enddo
!enddo
!do j = 1, classical%natoms
!  write(888,'(a2, 6f15.8)') trim(classical%atom_kinds(j)), classical%force(j,:)
!enddo
!write(889,*) sum(path_integral%rp(:)%num_beads) + classical%natoms
!write(889,*)' dt = ', dt*(icount-1)
!do irp = 1, path_integral%num_rp
!  do j = 1, path_integral%rp(irp)%num_beads
!    write(889,'(a2, 6f15.8)') trim(path_integral%rp(irp)%particle_label(j)), path_integral%rp(irp)%x(j,:), &
!      path_integral%rp(irp)%p(j,:)
!  enddo
!enddo
!do j = 1, classical%natoms
!  write(889,'(a2, 6f15.8)') trim(classical%atom_kinds(j)), classical%x(j,:),classical%p(j,:)
!enddo

END SUBROUTINE vv_update_pi_normal_manual

! --------------------------------------------- !
! Subroutine that updates the momentum of the - !
! path integral beads using the --------------- !
! verlet scheme.                                !
! p_j(t + delt/2)  <- p_j(t) + delt/2 * F_j     !
! --------------------------------------------- !
SUBROUTINE update_pi_p(path_integral, dt)
implicit none

TYPE(path_integral_env)           ::path_integral
integer                          ::irp
integer                          ::ibead
real(kind=8)                     ::dt

do irp = 1, path_integral%num_rp
  path_integral%rp(irp)%p(:,:) = path_integral%rp(irp)%p(:,:) + dt * path_integral%rp(irp)%force(:,:)/2.0d0 
enddo

END SUBROUTINE update_pi_p
! -------------------------------------------- !

SUBROUTINE manual_normal_rotation(path_integral)
implicit none

TYPE(path_integral_env)   ::path_integral
integer          ::irp, ibead, jbead, nbeads
integer          ::i

real(kind=8),allocatable,dimension(:,:)    ::rot_matrix

do irp = 1, path_integral%num_rp
  
  nbeads = path_integral%rp(irp)%num_beads
  allocate(path_integral%rp(irp)%rot_matrix(nbeads,nbeads),rot_matrix(nbeads,nbeads))
  do ibead = 1, nbeads
    if ( (ibead-1) == 0 ) then
      rot_matrix(:,ibead) = 1.0d0/dsqrt(dble(real(nbeads)))
    else if ( (ibead-1) .le. nbeads/2-1) then
      do jbead = 1, nbeads
        rot_matrix(jbead,ibead) = dsqrt(2.0d0/dble(real(nbeads))) * dcos(2.0*pi*(ibead-1)*(jbead-1)/dble(real(nbeads)))
      enddo
    else if ( (ibead-1) == nbeads/2) then
      do jbead = 1, nbeads
        rot_matrix(jbead,ibead) = dsqrt(1.0d0/dble(real(nbeads))) * (-1.0d0)**(jbead-1)
      enddo
    else if ( (ibead-1) .le. nbeads-1) then
      do jbead = 1, nbeads
        rot_matrix(jbead,ibead) = dsqrt(2.0d0/dble(real(nbeads)))*dsin(2.0*pi*(ibead-1)*(jbead-1)/dble(real(nbeads)))
      enddo
    endif
  enddo

  !open(888,file='rotation_matrix.dat')
  !do i = 1, nbeads
  !  write(888,'(100f10.2)') rot_matrix(i,:)
  !enddo
  !close(888)

  path_integral%rp(irp)%rot_matrix = rot_matrix
  deallocate(rot_matrix)
  !do i_dim = 1, md%space_dim
  !  path_integral%rp(irp)%p(:,i_dim) = path_integral%rp(irp)%norm_p(:,i_dim)*rot_matrix(
  !MATRIX MUTLIPLICATION TO PROPOGATE THE NORMAL MODE COORDINATES.I
enddo

END SUBROUTINE manual_normal_rotation

SUBROUTINE manual_cart_to_normal(path_integral, md, istep)
implicit none

TYPE(path_integral_env)   ::path_integral
TYPE(md_info)             ::md
integer            ::irp,ispace
integer            ::j,k, i
integer            ::istep

!if (istep==1) then
  !open(777,file='rot_matrix_check.dat')
  !do i = 1, path_integral%rp(1)%num_beads
  !  write(777,'(100f12.4)') path_integral%rp(1)%rot_matrix(i,:)
  !enddo
  !close(777)
!endif

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

SUBROUTINE update_normal_coord(path_integral, md, dt)
implicit none

TYPE(path_integral_env)    ::path_integral
TYPE(md_info)              ::md
real(kind=8)     ::omega_n, omega_k
real(kind=8)     ::dt
integer          ::irp, ibead, nbeads, i_dim

real(kind=8),dimension(2,2)    ::free_rp_prop
real(kind=8),dimension(2)      ::old_norm, new_norm


do irp = 1, path_integral%num_rp
  !omega_n = (path_integral%rp(irp)%harm_mass/dble(path_integral%rp(irp)%num_beads))/(path_integral%rp(irp)%betaN**2*path_integral%rp(irp)%num_beads)   !CHECK!!
  !omega_n = dsqrt(path_integral%rp(irp)%num_beads/(path_integral%rp(irp)%betaN))   !CHECK!!
  !BAJ EDIT on FEB 27. CHECK!
  !omega_n = dsqrt(path_integral%rp(irp)%harm_mass/(path_integral%rp(irp)%kinetic_mass*&
  !  dble(path_integral%rp(irp)%num_beads)))/path_integral%rp(irp)%betaN
  omega_n = (dble(real(path_integral%rp(irp)%num_beads))/path_integral%rp(irp)%beta)&
        *dsqrt(path_integral%rp(irp)%harm_mass/path_integral%rp(irp)%kinetic_mass)
  nbeads = path_integral%rp(irp)%num_beads
  do ibead = 1, nbeads
    omega_k =2.0d0*omega_n*dsin( (ibead-1)*pi/dble(real(nbeads)) )
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


SUBROUTINE manual_normal_to_cart(path_integral, md, istep)
implicit none
TYPE(path_integral_env)   ::path_integral
TYPE(md_info)             ::md
integer              ::irp,ispace
integer              ::j,k, i
integer              ::istep
!if (istep==1) then
!  open(777,file='rot_matrix_check_2.dat')
!  do i = 1, path_integral%rp(1)%num_beads
!    write(777,'(100f12.4)') path_integral%rp(1)%rot_matrix(i,:)
!  enddo
!  close(777)
!endif

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




SUBROUTINE britta_velocities(classical, path_integral, md)
implicit none

TYPE(path_integral_env)         ::path_integral
TYPE(classical_type)            ::classical
TYPE(md_info)                   ::md
real(kind=8), allocatable, dimension(:)    ::vel_array
integer                         ::i, irp
integer                         ::icount
integer, allocatable,dimension(:)  ::his
real(kind=8)                    ::xmin, xmax, delx
integer                         ::xnum

allocate(vel_array(md%space_dim*(classical%natoms+path_integral%total_beads)))
vel_array = 0.0d0

open(199, file='pi_data.out')
do irp = 1, path_integral%num_rp
  do i = 1, path_integral%rp(irp)%num_beads
    write(199,*) path_integral%rp(irp)%p(i,:), path_integral%rp(irp)%kinetic_mass
  enddo
enddo
close(199)

open(198,file='classical_data.out')
do i = 1, classical%natoms
  write(198, *) classical%p(i,1:md%space_dim), classical%mass(i)
enddo
close(198)
icount = 0
do irp = 1, path_integral%num_rp
  do i = 1, path_integral%rp(irp)%num_beads
    icount = icount + 1
    vel_array( (icount-1)*md%space_dim+1: icount*md%space_dim) = path_integral%rp(irp)%p(i,:)/path_integral%rp(irp)%kinetic_mass
  enddo
enddo
if (icount /= path_integral%total_beads) write(*,*)' ERROR IN VEL-ARRAY '

do i = 1, classical%natoms
  vel_array( (i-1+path_integral%total_beads)*md%space_dim+1: (i+path_integral%total_beads)*md%space_dim ) = &
    classical%p(i,1:md%space_dim)/classical%mass(i)
enddo


xmin = -0.004
xmax = 0.004
xnum = 50
allocate(his(0:xnum))
CALL histogram(path_integral%total_beads*md%space_dim, vel_array(1:path_integral%total_beads*md%space_dim), xmin, xmax, xnum, his) 

open(199, file='his_pi.out')
delx = (xmax-xmin)/real(xnum-1)
do i = 0, xnum
  write(199,*) i, xmin + (i-1)*delx, his(i)
enddo
close(199)


END SUBROUTINE britta_velocities


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



SUBROUTINE thermostat_pile_g(classical, path_integral, md, dt, para_var)
implicit none

TYPE(classical_type)   ::classical
TYPE(path_integral_env)   ::path_integral
TYPE(md_info)             ::md
TYPE(parallel_env)     ::para_var
real(kind=8)        :: dt
integer             :: i, j, k, icount
real(kind=8)        :: rnd1, c1, c2
real(kind=8)        :: com

real(kind=8)        :: kconst, sum_rnd, rnd_0
real(kind=8)        :: alpha2, alpha
integer             :: root, ierror
include 'mpif.h'

if (para_var%my_rank == 0) then
do i = 1, classical%natoms
  do j = 1, md%space_dim
    CALL rnd_gauss(0.0d0, 1.0d0, rnd1)
    c1 = dexp(-1.0d0*dt*0.5d0*classical%lang_atom_coeff(i))
    c2 = dsqrt( 1 - c1**2)
    classical%p(i,j) = classical%p(i,j)*c1 + dsqrt(classical%mass(i)/classical%beta)*c2*rnd1
  enddo
enddo
if (md%n_fixed_atoms > 0) then
  icount = md%n_fixed_atoms
  do while (md%fixed_list(icount,1) > path_integral%num_rp)
    if (md%fixed_list(icount,2) == 1) then
      classical%p(md%fixed_list(icount,1)-path_integral%num_rp,1) = 0.0d0
    endif
    if (md%fixed_list(icount,2) == 1) then
      classical%p(md%fixed_list(icount,1)-path_integral%num_rp,2) = 0.0d0
    endif
    if (md%fixed_list(icount,2) == 1) then
      classical%p(md%fixed_list(icount,1)-path_integral%num_rp,3) = 0.0d0
    endif
    icount = icount - 1
  enddo
endif
!calculcate transform matrix
CALL manual_normal_rotation(path_integral)
!transform cart to normal
CALL manual_cart_to_normal(path_integral, md, 1)
kconst = 0.0d0
do i = 1, path_integral%num_rp
    kconst = kconst+ dot_product(path_integral%rp(i)%norm_p(1,:), path_integral%rp(i)%norm_p(1,:)) &
      /(2.0d0 * path_integral%rp(i)%kinetic_mass)
enddo
c1 = dexp(-1.0d0*dt*0.5d0*path_integral%rp(1)%lang_coeff(1))
sum_rnd = 0.0d0
do i = 1, path_integral%num_rp
  do k = 1, md%space_dim
    CALL rnd_gauss(0.0d0, 1.0d0, rnd1)
    sum_rnd = sum_rnd + rnd1**2
    if (i==1) rnd_0 = rnd1
  enddo
enddo
alpha2 = c1 +( ( (1-c1)*sum_rnd)/(2.0d0*path_integral%rp(1)%beta * kconst)) + &
   2.0d0*rnd_0*dsqrt( (c1*(1-c1))/(2.0*path_integral%rp(1)%beta * kconst) )
alpha = sign(dsqrt(alpha2), rnd_0 + dsqrt( (2.0*path_integral%rp(1)%beta * kconst)/(1-c1)))
do i = 1, path_integral%num_rp
  do j = 1, path_integral%rp(i)%num_beads
    do k = 1, md%space_dim
      if (j==1) then
        path_integral%rp(i)%norm_p(j,k) = alpha * path_integral%rp(i)%norm_p(j,k)
      else
        CALL rnd_gauss(0.0d0, 1.0d0, rnd1)
        c1 = dexp(-1.0d0*dt*0.5d0*path_integral%rp(i)%lang_coeff(j))
        c2 = dsqrt( 1 - c1**2)
        path_integral%rp(i)%norm_p(j,k) = path_integral%rp(i)%norm_p(j,k)*c1 + &
          dsqrt(path_integral%rp(i)%kinetic_mass/path_integral%rp(i)%beta)*c2*rnd1     !BAJ CHANGE FROM BEFORE
      endif
    enddo
  enddo
enddo
CALL manual_normal_to_cart(path_integral, md, 1)
if (md%centroid) then
  icount = 1
  do while (md%centroid_list(icount,1) <= path_integral%num_rp)
    if (md%centroid_list(icount,2) == 1) then
      com = sum(path_integral%rp(md%centroid_list(icount,1))%p(:,1))/&
        dble(real(path_integral%rp(md%centroid_list(icount,1))%num_beads))
      path_integral%rp(md%centroid_list(icount,1))%p(:,1) = path_integral%rp(md%centroid_list(icount,1))%p(:,1) - com
    endif
    if (md%centroid_list(icount,3) == 1) then
      com = sum(path_integral%rp(md%centroid_list(icount,1))%p(:,2))/&
        dble(real(path_integral%rp(md%centroid_list(icount,1))%num_beads))
      path_integral%rp(md%centroid_list(icount,1))%p(:,2) = path_integral%rp(md%centroid_list(icount,1))%p(:,2) - com
    endif
    if (md%centroid_list(icount,4) == 1) then
      com = sum(path_integral%rp(md%centroid_list(icount,1))%p(:,3))/&
        dble(real(path_integral%rp(md%centroid_list(icount,1))%num_beads))
      path_integral%rp(md%centroid_list(icount,1))%p(:,3) = path_integral%rp(md%centroid_list(icount,1))%p(:,3) - com
    endif
    icount = icount+1
    if (icount > md%n_centroids) EXIT
  enddo
endif
if (md%n_fixed_atoms > 0) then
  icount = 1
  do while (md%fixed_list(icount,1) <= path_integral%num_rp)
  write(*,*)' fixed RP atoms'
    if (md%fixed_list(icount,2) == 1) then
      path_integral%rp(md%fixed_list(icount,1))%p(:,1) = 0.0d0
    endif
    if (md%fixed_list(icount,3) == 1) then
      path_integral%rp(md%fixed_list(icount,1))%p(:,2) = 0.0d0
    endif
    if (md%fixed_list(icount,4) == 1) then
      path_integral%rp(md%fixed_list(icount,1))%p(:,3) = 0.0d0
    endif
    icount = icount+1
    if (icount > md%n_fixed_atoms) EXIT
  enddo
endif
endif

root = 0
call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  do i = 1, path_integral%num_rp
    do j = 1, path_integral%rp(i)%num_beads
      !do k = 1, md%space_dim
      !write(*,*)' i, j, ', i, j, para_var%my_rank
      !write(*,*)' value ', path_integral%rp(i)%initialized_velocities(j,:)
      call MPI_bcast(path_integral%rp(i)%p(j,1:md%space_dim), md%space_dim, MPI_DOUBLE_PRECISION,root, &
        MPI_COMM_WORLD, ierror)
      !write(*,*)' value ', path_integral%rp(i)%initialized_velocities(j,:)
      !enddo
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    enddo
  enddo
  do i = 1, classical%natoms
    call MPI_bcast(classical%p(i,:), md%space_dim, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror) 
  enddo

!do i = 0, para_var%size_rank-1
!  if (i == para_var%my_rank) then
!    write(*,*)' on rank ', para_var%my_rank, ' after rescaling '
!    write(*,*)' the velocities are ', path_integral%rp(1)%p(1,:)
!  endif
!  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
!enddo



END SUBROUTINE thermostat_pile_g
    

SUBROUTINE thermostat_pile(classical, path_integral, md, dt, para_var)
implicit none

TYPE(classical_type)   ::classical
TYPE(path_integral_env)   ::path_integral
TYPE(md_info)             ::md
TYPE(parallel_env)     ::para_var
real(kind=8)        :: dt, com
integer             :: i, j, k, icount
real(kind=8)        :: rnd1, c1, c2

integer             :: root, ierror
include 'mpif.h'

if (para_var%my_rank == 0) then
do i = 1, classical%natoms
  do j = 1, md%space_dim
    CALL rnd_gauss(0.0d0, 1.0d0, rnd1)
    c1 = dexp(-1.0d0*dt*0.5d0*classical%lang_atom_coeff(i))
    c2 = dsqrt( 1 - c1**2)
    classical%p(i,j) = classical%p(i,j)*c1 + dsqrt(classical%mass(i)/classical%beta)*c2*rnd1
  enddo
enddo
if (md%n_fixed_atoms > 0) then
  icount = md%n_fixed_atoms
  do while (md%fixed_list(icount,1) > path_integral%num_rp)
    if (md%fixed_list(icount,2) == 1) then
      classical%p(md%fixed_list(icount,1)-path_integral%num_rp,1) = 0.0d0
    endif
    if (md%fixed_list(icount,2) == 1) then
      classical%p(md%fixed_list(icount,1)-path_integral%num_rp,2) = 0.0d0
    endif
    if (md%fixed_list(icount,2) == 1) then
      classical%p(md%fixed_list(icount,1)-path_integral%num_rp,3) = 0.0d0
    endif
    icount = icount - 1
  enddo
endif
if (path_integral%num_rp > 0) then
!calculcate transform matrix
CALL manual_normal_rotation(path_integral)
!transform cart to normal
CALL manual_cart_to_normal(path_integral, md, 1)
do i = 1, path_integral%num_rp
  do j = 1, path_integral%rp(i)%num_beads
    do k = 1, md%space_dim
      CALL rnd_gauss(0.0d0, 1.0d0, rnd1)
      c1 = dexp(-1.0d0*dt*0.5d0*path_integral%rp(i)%lang_coeff(j))
      c2 = dsqrt( 1 - c1**2)
      path_integral%rp(i)%norm_p(j,k) = path_integral%rp(i)%norm_p(j,k)*c1 + &
        dsqrt(path_integral%rp(i)%kinetic_mass/path_integral%rp(i)%beta)*c2*rnd1     !BAJ CHANGE FROM BEFORE
    enddo
  enddo
enddo
CALL manual_normal_to_cart(path_integral, md, 1)
if (md%centroid) then
  icount = 1
  do while (md%centroid_list(icount,1) <= path_integral%num_rp)
    if (md%centroid_list(icount,2) == 1) then
      com = sum(path_integral%rp(md%centroid_list(icount,1))%p(:,1))/&
        dble(real(path_integral%rp(md%centroid_list(icount,1))%num_beads))
      path_integral%rp(md%centroid_list(icount,1))%p(:,1) = path_integral%rp(md%centroid_list(icount,1))%p(:,1) - com
    endif
    if (md%centroid_list(icount,3) == 1) then
      com = sum(path_integral%rp(md%centroid_list(icount,1))%p(:,2))/&
        dble(real(path_integral%rp(md%centroid_list(icount,1))%num_beads))
      path_integral%rp(md%centroid_list(icount,1))%p(:,2) = path_integral%rp(md%centroid_list(icount,1))%p(:,2) - com
    endif
    if (md%centroid_list(icount,4) == 1) then
      com = sum(path_integral%rp(md%centroid_list(icount,1))%p(:,3))/&
        dble(real(path_integral%rp(md%centroid_list(icount,1))%num_beads))
      path_integral%rp(md%centroid_list(icount,1))%p(:,3) = path_integral%rp(md%centroid_list(icount,1))%p(:,3) - com
    endif
    icount = icount+1
    if (icount > md%n_centroids) EXIT
  enddo
endif
if (md%n_fixed_atoms > 0) then
  icount = 1
  do while (md%fixed_list(icount,1) <= path_integral%num_rp)
  write(*,*)' fixed RP atoms'
    if (md%fixed_list(icount,2) == 1) then
      path_integral%rp(md%fixed_list(icount,1))%p(:,1) = 0.0d0
    endif
    if (md%fixed_list(icount,3) == 1) then
      path_integral%rp(md%fixed_list(icount,1))%p(:,2) = 0.0d0
    endif
    if (md%fixed_list(icount,4) == 1) then
      path_integral%rp(md%fixed_list(icount,1))%p(:,3) = 0.0d0
    endif
    icount = icount+1
    if (icount > md%n_fixed_atoms) EXIT
  enddo
endif
endif
endif

root = 0
call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  do i = 1, path_integral%num_rp
    do j = 1, path_integral%rp(i)%num_beads
      !do k = 1, md%space_dim
      !write(*,*)' i, j, ', i, j, para_var%my_rank
      !write(*,*)' value ', path_integral%rp(i)%initialized_velocities(j,:)
      call MPI_bcast(path_integral%rp(i)%p(j,1:md%space_dim), md%space_dim, MPI_DOUBLE_PRECISION,root, &
        MPI_COMM_WORLD, ierror)
      !write(*,*)' value ', path_integral%rp(i)%initialized_velocities(j,:)
      !enddo
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    enddo
  enddo
  do i = 1, classical%natoms
    call MPI_bcast(classical%p(i,:), md%space_dim, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror) 
  enddo

!do i = 0, para_var%size_rank-1
!  if (i == para_var%my_rank) then
!    write(*,*)' on rank ', para_var%my_rank, ' after rescaling '
!    write(*,*)' the velocities are ', path_integral%rp(1)%p(1,:)
!  endif
!  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
!enddo



END SUBROUTINE thermostat_pile
    
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


END MODULE md_run
