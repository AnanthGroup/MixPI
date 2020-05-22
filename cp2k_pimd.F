! ---------------------------------- !
! -----CP2K_PIMD.F------------------ !
! The main program that mixes the -- !
! path integral representation with  !
! classical MD integration by CP2K - !
! ---------------------------------- !
! B. Johnson   12/15/19 ------------ !
! ---------------------------------- !

PROGRAM cp2k_pimd
use cp2k_run_forces
use pimd_mod
implicit none

integer     ::i,j

integer                               ::irp

! ----------------------------------- !
! set up CP2K
write(*,*)' linking to CP2K and creating force environment '
CALL set_up_cp2k
! ----------------------------------- !


! ------------------------------------------- !
!read in MD information from md_run.inp
!right now you need to confirm that any
!repeat info is repeated in both the CP2K.inp
!and md_run.inp
write(*,*)' read in the MD input '
CALL read_in_md
! ------------------------------------------- !


! ------------------------------------------- !
!initialize and read in the PI variables
write(*,*)' initialize the path integral variables '
CALL pimd_initialize
! ------------------------------------------- !



! ------------------------------------------- !
!initalize and read in the classical atom 
!variables.
write(*,*)' now initalize classical atoms '
if (trim(md%obtain_force) == 'cp2k') then
  write(*,*)' GO '
  CALL classical_atom_initialize
endif
! ------------------------------------------- !


! ------------------------------------------- !
!if initialize, complete a series of mc steps 
!to initialize the positions of the ring polymer
!
!if not, set the initial positions to either
!the values read in the position files or all
!beads are initialize at the centroid
!EDIT NEEDED!!!!!!!!!!!!!!!!!!!!!!!!!
!Either need to initialize all rp at same time
!or figure out how only doing one at a time would
!work....seems like you need to move all of them!!
if(path_integral%initialize) THEN
  do irp = 1, path_integral%num_rp
    CALL initialize_mc
  enddo
else
  do i = 1, path_integral%num_rp
    ALLOCATE(path_integral%rp(i)%initialized_positions(path_integral%rp(i)%num_beads,md%space_dim))
    if(path_integral%rp(i)%coord) then
      do j = 1, path_integral%rp(i)%num_beads
        path_integral%rp(i)%initialized_positions(j,:) = path_integral%rp(i)%i_position(j,:)
      enddo
    else
      do j = 1, path_integral%rp(i)%num_beads
        path_integral%rp(i)%initialized_positions(j,:) = path_integral%rp(i)%i_centroid(:)
      enddo
    endif
  enddo
endif
! ------------------------------------------- !




! ------------------------------------------- !
!call the pimd md routine.
!this will conduct a MD run for the ring polymer
!and classical atoms (if present).
CALL run_md
! ------------------------------------------- !




! ------------------------------------------- !
! if used classical atoms, deallocate all required
! arrays.
if (trim(md%obtain_force) == 'cp2k') then
  CALL classical_atom_destroy
endif
! ------------------------------------------- !


! ------------------------------------------- !
!deallocates all RP arrays
call pimd_destory
! ------------------------------------------- !



! ------------------------------------------- !
!deallocates all cp2k arrays and removes the 
!force_env information.
CALL destroy_cp2k
! ------------------------------------------- !


END PROGRAM cp2k_pimd 

 
