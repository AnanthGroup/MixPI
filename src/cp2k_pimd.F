! ---------------------------------- !
! -----CP2K_PIMD.F------------------ !
! The main program that mixes the -- !
! path integral representation with  !
! classical MD integration by CP2K - !
! ---------------------------------- !
! B. Johnson   12/15/19 ------------ !
! ---------------------------------- !
! Major changes made by B. Johnson - !
! when CP2K version 9 update. Also - !
! moved all inter molecule (both --- !
! classical and RP force calculations!
! to the CP2K subroutines. --------- !
! B. Johnson 01/20/23 -------------- !
! ---------------------------------- !

PROGRAM cp2k_pimd
use input_setup
use input_destroy
use variable_types
use md_run
implicit none

integer     ::i,j
integer                               ::irp
TYPE(path_integral_env)                ::path_integral
TYPE(classical_type)                   ::classical
TYPE(md_info)                           ::md
TYPE(parallel_env)                     ::para_var

integer     :: rate, itime_start, itime_end

include 'mpif.h'

! ------------------------------------------- !
!read in MD information from md_run.inp
!right now you need to confirm that any
!repeat info is repeated in both the CP2K.inp
!and md_run.inp
!write(*,*)' read in the MD input '
!CALL read_in_md(md)
! ------------------------------------------- !

! ------------------------------------------- !
!initalize and read in the classical atom 
!variables.
!if (trim(md%obtain_force) == 'cp2k') then
  !write(*,*)' GO '
  ! ----------------------------------- !
  ! set up CP2K
  CALL set_up_cp2k(para_var)
  if (para_var%my_rank == 0) then
    CALL system_clock(count_rate=rate)
    call system_clock(itime_start)
  endif
! ------------------------------------------- !
!read in MD information from md_run.inp
!right now you need to confirm that any
!repeat info is repeated in both the CP2K.inp
!and md_run.inp
  !write(*,*)' read in the MD input '
  CALL read_in_md(md)
! ------------------------------------------- !
  ! ----------------------------------- !
  CALL cp2k_pimd_initialize(path_integral, classical, md, para_var)
!else
!  write(*,*)' this obtain force option is not currently supported '
!endif
! ------------------------------------------- !


! ------------------------------------------- !
!call the pimd md routine.
!this will conduct a MD run for the ring polymer
!and classical atoms (if present).
write(*,*) 'run_md '
CALL run_md(classical, path_integral, md, para_var)
! ------------------------------------------- !




! ------------------------------------------- !
! if used classical atoms, deallocate all required
! arrays.
if (trim(md%obtain_force) == 'cp2k') then
  CALL classical_atom_destroy(classical)
  ! ------------------------------------------- !
  !deallocates all cp2k arrays and removes the 
  !force_env information.
    CALL destroy_cp2k
  ! ------------------------------------------- !
endif
! ------------------------------------------- !


! ------------------------------------------- !
!deallocates all RP arrays
call pimd_destroy(path_integral,md)
! ------------------------------------------- !

if (para_var%my_rank == 0) then
  call system_clock(itime_end)

  open(111,file='runtime_info.dat',status="old", position="append")
    write(111,*)' complete run time is ', itime_end - itime_start, real(itime_end-itime_start)/real(rate)
  close(111)
endif


END PROGRAM cp2k_pimd 

 
