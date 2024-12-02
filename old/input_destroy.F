MODULE input_destroy
use libcp2k
use variable_types 
implicit none


contains


SUBROUTINE destroy_cp2k
! ------------------------------------------------ !
! CALL cp2 k_finalize() from the CP2K library ----- !
! This subroutine deallocates/nullifies the ------ !
! structures to prevent leakage. Pairs with ------ !
! cp2k_init() ------------------------------------ !
CALL cp2k_finalize()
! ------------------------------------------------ !

END SUBROUTINE destroy_cp2k





SUBROUTINE classical_atom_destroy(classical)
implicit none
TYPE(classical_type)       ::classical
integer               :: ierror
include 'mpif.h'

CALL cp2k_destroy_force_env(classical%cp2k_env_flag)
DEALLOCATE(classical%i_position)

!CALL MPI_FINALIZE(ierror)

END SUBROUTINE


! ------------------------------------------------------ !
! Subroutine that will deallocate/remove all of the pimd !
! structures to prevent any memory leakage, etc. ------- !
! ------------------------------------------------------ !
SUBROUTINE pimd_destroy(path_integral,md)
implicit none

TYPE(path_integral_env)     ::path_integral
TYPE(md_info)                          ::md
integer                     ::i


do i = 1, path_integral%num_rp
  DEALLOCATE(path_integral%rp(i)%i_position)
  DEALLOCATE(path_integral%rp(i)%i_centroid)
  DEALLOCATE(path_integral%rp(i)%initialized_positions)
  DEALLOCATE(path_integral%rp(i)%initialized_velocities)
  DEALLOCATE(path_integral%rp(i)%force)
  DEALLOCATE(path_integral%rp(i)%p)
  DEALLOCATE(path_integral%rp(i)%x)
  if (trim(adjustl(md%print_level)) == 'high') DEALLOCATE(path_integral%rp(i)%traj)
enddo

DEALLOCATE(path_integral%rp)

END SUBROUTINE pimd_destroy
! ----------------------------------------------------- !






END MODULE input_destroy
