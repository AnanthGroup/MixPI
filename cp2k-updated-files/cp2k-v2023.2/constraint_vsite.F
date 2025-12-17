!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to handle the virtual site constraint/restraint
!> \par History
!>      Teodoro Laino [tlaino] 12.2008 - Preparing for VIRTUAL SITE constraints
!>                                       (patch by Marcel Baer)
! **************************************************************************************************
MODULE constraint_vsite
   USE atomic_kind_list_types,          ONLY: atomic_kind_list_type  !BAJ ADDED VS
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_type
   USE kinds,                           ONLY: dp
   USE molecule_kind_list_types,        ONLY: molecule_kind_list_type
   USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                              molecule_kind_type,&
                                              vsite_constraint_type
   USE molecule_list_types,             ONLY: molecule_list_type
   USE molecule_types,                  ONLY: get_molecule,&
                                              global_constraint_type,&
                                              molecule_type
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_types,                  ONLY: particle_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: shake_vsite_int, &
             shake_vsite_ext, &
             vsite_force_control

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'constraint_vsite'

CONTAINS

! **************************************************************************************************
!> \brief control force distribution for virtual sites
!> \param force_env ...
!> \date 12.2008
!> \par History
!>      - none
!> \author Marcel Baer
! **************************************************************************************************
   SUBROUTINE vsite_force_control(force_env)
      TYPE(force_env_type), POINTER                      :: force_env

      INTEGER                                            :: i, ikind, imol, nconstraint, nkind, &
                                                            nmol_per_kind, nvsitecon, &
                                                            nparticle_kind, nparticle_local, &          !BAJ ADDED VS
                                                            iparticle, iparticle_kind, iparticle_local  !BAJ ADDED VS
      LOGICAL                                            :: do_ext_constraint
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: my_pos    !BAJ ADDED VS
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(distribution_1d_type), POINTER                :: local_molecules, local_particles
      TYPE(global_constraint_type), POINTER              :: gci
      TYPE(molecule_kind_list_type), POINTER             :: molecule_kinds
      TYPE(molecule_kind_type), DIMENSION(:), POINTER    :: molecule_kind_set
      TYPE(molecule_kind_type), POINTER                  :: molecule_kind
      TYPE(molecule_list_type), POINTER                  :: molecules
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(molecule_type), POINTER                       :: molecule
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(atomic_kind_list_type), POINTER               :: atomic_kinds

      NULLIFY (gci, subsys, local_molecules, local_particles, &
               molecule_kinds,my_pos)                                                      !BAJ ADDED VS

      CALL force_env_get(force_env=force_env, subsys=subsys)

      CALL cp_subsys_get(subsys=subsys, local_particles=local_particles, &
                         particles=particles, local_molecules=local_molecules, &
                         atomic_kinds=atomic_kinds, &                                        !BAJ ADDED VS
                         molecule_kinds=molecule_kinds, gci=gci, molecules=molecules)
      nparticle_kind = atomic_kinds%n_els

      molecule_kind_set => molecule_kinds%els
      molecule_set => molecules%els
      particle_set => particles%els
      nkind = SIZE(molecule_kind_set)
      ! Intermolecular Virtual Site Constraints
      do_ext_constraint = .FALSE.
      IF (ASSOCIATED(gci)) THEN
         do_ext_constraint = (gci%ntot /= 0)
      END IF
      IF (.NOT. ASSOCIATED(my_pos)) THEN
         ALLOCATE (my_pos(3, particles%n_els))
         my_pos = 0.0_dp
         DO iparticle_kind = 1, nparticle_kind
            nparticle_local = local_particles%n_el(iparticle_kind)
            DO iparticle_local = 1, nparticle_local
               iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
               my_pos(:, iparticle) = particles%els(iparticle)%r(:)
            END DO
         END DO
      END IF
      ! Intramolecular Virtual Site Constraints
      MOL: DO ikind = 1, nkind
         nmol_per_kind = local_molecules%n_el(ikind)
         DO imol = 1, nmol_per_kind
            i = local_molecules%list(ikind)%array(imol)
            molecule => molecule_set(i)
            molecule_kind => molecule%molecule_kind
            CALL get_molecule_kind(molecule_kind, nconstraint=nconstraint, nvsite=nvsitecon)
            IF (nconstraint == 0) CYCLE
            IF (nvsitecon /= 0) &
               !CALL force_vsite_int(molecule, particle_set, pos=my_pos)
               CALL force_vsite_int(molecule, particle_set)   !BAJ EDIT
         END DO
      END DO MOL
      ! Intermolecular Virtual Site Constraints
      IF (do_ext_constraint) THEN
         IF (gci%nvsite /= 0) &
            !CALL force_vsite_ext(gci, particle_set, pos=my_pos)
            CALL force_vsite_ext(gci, particle_set)  !BAJ EDIT
      END IF

      DEALLOCATE(my_pos)  !BAJ ADDED

   END SUBROUTINE vsite_force_control

! **************************************************************************************************
!> \brief Intramolecular virtual site
!> \param molecule ...
!> \param pos ...
!> \par History
!>      12.2008 Marcel Baer
! **************************************************************************************************
   SUBROUTINE shake_vsite_int(molecule, pos)
      TYPE(molecule_type), POINTER                       :: molecule
      REAL(KIND=dp), INTENT(INOUT)                       :: pos(:, :)

      INTEGER                                            :: first_atom, nvsite
      TYPE(molecule_kind_type), POINTER                  :: molecule_kind
      TYPE(vsite_constraint_type), POINTER               :: vsite_list(:)

      molecule_kind => molecule%molecule_kind
      CALL get_molecule_kind(molecule_kind, nvsite=nvsite, vsite_list=vsite_list)
      CALL get_molecule(molecule, first_atom=first_atom)
      ! Real Shake
      CALL shake_vsite_low(vsite_list, nvsite, first_atom, pos)

   END SUBROUTINE shake_vsite_int

! **************************************************************************************************
!> \brief Intramolecular virtual site
!> \param gci ...
!> \param pos ...
!> \par History
!>      12.2008 Marcel Baer
! **************************************************************************************************
   SUBROUTINE shake_vsite_ext(gci, pos)

      TYPE(global_constraint_type), POINTER              :: gci
      REAL(KIND=dp), INTENT(INOUT)                       :: pos(:, :)

      INTEGER                                            :: first_atom, nvsite
      TYPE(vsite_constraint_type), POINTER               :: vsite_list(:)

      first_atom = 1
      nvsite = gci%nvsite
      vsite_list => gci%vsite_list
      ! Real Shake
      CALL shake_vsite_low(vsite_list, nvsite, first_atom, pos)

   END SUBROUTINE shake_vsite_ext

! **************************************************************************************************
!> \brief ...
!> \param vsite_list ...
!> \param nvsite ...
!> \param first_atom ...
!> \param pos ...
!> \par History
!>      12.2008 Marcel Bear
! **************************************************************************************************
   SUBROUTINE shake_vsite_low(vsite_list, nvsite, first_atom, pos)
      TYPE(vsite_constraint_type)                        :: vsite_list(:)
      INTEGER, INTENT(IN)                                :: nvsite, first_atom
      REAL(KIND=dp), INTENT(INOUT)                       :: pos(:, :)
      REAL(KIND=dp)                                      :: wbc, wdc     !ADDED BAJ VS

      INTEGER                                            :: iconst, index_a, index_b, index_c, &
                                                            index_d
      REAL(KIND=dp), DIMENSION(3)                        :: r1, r2

      DO iconst = 1, nvsite
         IF (vsite_list(iconst)%restraint%active) CYCLE
         index_a = vsite_list(iconst)%a + first_atom - 1
         index_b = vsite_list(iconst)%b + first_atom - 1
         index_c = vsite_list(iconst)%c + first_atom - 1
         index_d = vsite_list(iconst)%d + first_atom - 1

         r1(:) = pos(:, index_b) - pos(:, index_c)
         r2(:) = pos(:, index_d) - pos(:, index_c)
         if (abs(vsite_list(iconst)%wdc+1)<1.0d-12) then            !!BAJ EDIT: introduction of flexible TIP4P
           if (iconst == 1) write(*,*)' this is flexible loop '
           wbc = (vsite_list(iconst)%wbc)
           wdc = (1.0d0 - wbc)/2.0d0
           pos(:, index_a) = wbc*pos(:, index_c) + wdc*(pos(:,index_b) + pos(:,index_d))    !I THINK THIS EXPRESSION IS OFF. CHECK MATH
         else
         !if (abs(vsite_list(iconst)%wdc+1)<1.0d-12) then            !!BAJ EDIT: introduction of flexible TIP4P
         !  if (iconst == 1) write(*,*)' this is flexible loop '
         !  wbc = (vsite_list(iconst)%wbc)/(dsqrt(dot_product(r1,r1))*2.0_dp)
         !  wdc = (vsite_list(iconst)%wbc)/(dsqrt(dot_product(r2,r2))*2.0_dp)
         !  pos(:, index_a) = pos(:, index_c) + wbc*r1(:) + wdc*r2(:)     !I THINK THIS EXPRESSION IS OFF. CHECK MATH
         !else
           pos(:, index_a) = pos(:, index_c) + vsite_list(iconst)%wbc*r1(:) + &
                           vsite_list(iconst)%wdc*r2(:)
         endif
      END DO
   END SUBROUTINE shake_vsite_low

! **************************************************************************************************
!> \brief Intramolecular virtual site
!> \param molecule ...
!> \param particle_set ...
!> \par History
!>      12.2008 Marcel Bear
! **************************************************************************************************
   !SUBROUTINE force_vsite_int(molecule, particle_set,pos)               !BAJ   EDIT VS
   SUBROUTINE force_vsite_int(molecule, particle_set)               !BAJ   EDIT VS
      TYPE(molecule_type), POINTER                       :: molecule
      TYPE(particle_type), POINTER                       :: particle_set(:)

      INTEGER                                            :: first_atom, iconst, index_a, index_b, &
                                                            index_c, index_d, nvsite
      !REAL(KIND=dp), INTENT(INOUT)                       :: pos(:, :)     !BAJ ADDED VS
      REAL(KIND=dp)                                      :: wb, wc, wd
      TYPE(molecule_kind_type), POINTER                  :: molecule_kind
      TYPE(vsite_constraint_type), POINTER               :: vsite_list(:)
      !REAL(KIND=dp), DIMENSION(3)                        :: r12, r13, norm_r12, norm_r13
      !REAL(KIND=dp)                                      :: norm12, norm13, fa_dot_r12, fa_dot_r13

      molecule_kind => molecule%molecule_kind
      CALL get_molecule_kind(molecule_kind, nvsite=nvsite, vsite_list=vsite_list)
      CALL get_molecule(molecule, first_atom=first_atom)

      DO iconst = 1, nvsite
         IF (vsite_list(iconst)%restraint%active) CYCLE
         !write(*,*)' first atom ', first_atom, iconst
         index_a = vsite_list(iconst)%a + first_atom - 1
         index_b = vsite_list(iconst)%b + first_atom - 1
         index_c = vsite_list(iconst)%c + first_atom - 1
         index_d = vsite_list(iconst)%d + first_atom - 1

         if (abs(vsite_list(iconst)%wdc+1)<1.0d-12) then            !!BAJ EDIT: introduction of flexible virtual site
         !  if (iconst == 1) write(*,*)' this is flexible loop - force_vsite_int'
         !    r12(:) = pos(:, index_b) - pos(:, index_c)
         !    r13(:) = pos(:, index_d) - pos(:, index_c)
         !    norm12 = dsqrt(dot_product(r12,r12))
         !    norm13 = dsqrt(dot_product(r13,r13))
         !    norm_r12(:) = r12(:)/norm12
         !    norm_r13(:) = r13(:)/norm13
         !    fa_dot_r12 = dot_product(particle_set(index_a)%f(:),norm_r12)
         !    fa_dot_r13 = dot_product(particle_set(index_a)%f(:),norm_r13)
         !    !fa_dot_r12 = 0.0d0
         !    !fa_dot_r13  = 0.0d0
         !    wb = (vsite_list(iconst)%wbc)
         !    wd = (vsite_list(iconst)%wbc)
         !    particle_set(index_b)%f(:) = particle_set(index_b)%f(:) + (wb/(2.0_dp*norm12)) *&
         !        (particle_set(index_a)%f(:) - fa_dot_r12 *norm_r12(:))
         !    particle_set(index_d)%f(:) = particle_set(index_d)%f(:) + (wd/(2.0_dp*norm13)) *&
         !        (particle_set(index_a)%f(:) - fa_dot_r13 *norm_r13(:))
         !    particle_set(index_c)%f(:) = particle_set(index_c)%f(:) - (wb/(2.0_dp*norm12)) *&
         !        (particle_set(index_a)%f(:) - fa_dot_r12 *norm_r12(:)) - &
         !        (wd/(2.0_dp*norm13)) * (particle_set(index_a)%f(:) - fa_dot_r13 *norm_r13(:)) + particle_set(index_a)%f(:)
         !    particle_set(index_a)%f(:) = 0.0_dp
         !
         !     
             wb = vsite_list(iconst)%wbc
             wd = (1.0d0-wb)/2.0d0
             particle_set(index_b)%f(:) = particle_set(index_b)%f(:) + wd*particle_set(index_a)%f(:)
             particle_set(index_c)%f(:) = particle_set(index_c)%f(:) + wb*particle_set(index_a)%f(:)
             particle_set(index_d)%f(:) = particle_set(index_d)%f(:) + wd*particle_set(index_a)%f(:)
             particle_set(index_a)%f(:) = 0.0_dp
           else
             wb = vsite_list(iconst)%wbc
             wd = vsite_list(iconst)%wdc
             wc = 1.0_dp - vsite_list(iconst)%wbc - vsite_list(iconst)%wdc
             particle_set(index_b)%f(:) = particle_set(index_b)%f(:) + wb*particle_set(index_a)%f(:)
             particle_set(index_c)%f(:) = particle_set(index_c)%f(:) + wc*particle_set(index_a)%f(:)
             particle_set(index_d)%f(:) = particle_set(index_d)%f(:) + wd*particle_set(index_a)%f(:)
             particle_set(index_a)%f(:) = 0.0_dp
           endif

      END DO

   END SUBROUTINE force_vsite_int

! **************************************************************************************************
!> \brief Intramolecular virtual site
!> \param gci ...
!> \param particle_set ...
!> \par History
!>      12.2008 Marcel Bear
! **************************************************************************************************
   !SUBROUTINE force_vsite_ext(gci, particle_set, pos)     !BAJ EDIT VS
   SUBROUTINE force_vsite_ext(gci, particle_set)     !BAJ EDIT VS
      TYPE(global_constraint_type), POINTER              :: gci
      TYPE(particle_type), POINTER                       :: particle_set(:)

      INTEGER                                            :: first_atom, iconst, index_a, index_b, &
                                                            index_c, index_d, nvsite
      !REAL(KIND=dp), INTENT(INOUT)                       :: pos(:, :)     !BAJ ADDED VS
      REAL(KIND=dp)                                      :: wb, wc, wd
      TYPE(vsite_constraint_type), POINTER               :: vsite_list(:)
      !REAL(KIND=dp), DIMENSION(3)                        :: r12, r13, norm_r12, norm_r13
      !REAL(KIND=dp)                                      :: norm12, norm13, fa_dot_r12, fa_dot_r13

      first_atom = 1
      nvsite = gci%nvsite
      vsite_list => gci%vsite_list
      ! Real Shake

      DO iconst = 1, nvsite
         IF (vsite_list(iconst)%restraint%active) CYCLE
         index_a = vsite_list(iconst)%a + first_atom - 1
         index_b = vsite_list(iconst)%b + first_atom - 1
         index_c = vsite_list(iconst)%c + first_atom - 1
         index_d = vsite_list(iconst)%d + first_atom - 1

         if (abs(vsite_list(iconst)%wdc+1)<1.0d-12) then            !!BAJ EDIT: introduction of flexible virtual
           !if (iconst == 1) write(*,*)' this is flexible loop - force_vsite_ext'
           !  r12(:) = pos(:, index_b) - pos(:, index_c)
           !  r13(:) = pos(:, index_d) - pos(:, index_c)
           !  norm12 = dsqrt(dot_product(r12,r12))
           !  norm13 = dsqrt(dot_product(r13,r13))
           !  norm_r12(:) = r12(:)/norm12
           !  norm_r13(:) = r13(:)/norm13
           !  fa_dot_r12 = dot_product(particle_set(index_a)%f(:),norm_r12)
           !  fa_dot_r13 = dot_product(particle_set(index_a)%f(:),norm_r13)
           !  !fa_dot_r12 = 0.0d0
           !  !fa_dot_r13  = 0.0d0
           !  wb = (vsite_list(iconst)%wbc)
           !  wd = (vsite_list(iconst)%wbc)
           !  particle_set(index_b)%f(:) = particle_set(index_b)%f(:) + (wb/(2.0_dp*norm12)) *&
           !      (particle_set(index_a)%f(:) - fa_dot_r12 *norm_r12(:))
           !  particle_set(index_d)%f(:) = particle_set(index_d)%f(:) + (wd/(2.0_dp*norm13)) *&
           !      (particle_set(index_a)%f(:) - fa_dot_r13 *norm_r13(:))
           !  particle_set(index_c)%f(:) = particle_set(index_c)%f(:) - (wb/(2.0_dp*norm12)) *&
           !      (particle_set(index_a)%f(:) - fa_dot_r12 *norm_r12(:)) - &
           !      (wd/(2.0_dp*norm13)) * (particle_set(index_a)%f(:) - fa_dot_r13 *norm_r13(:)) + particle_set(index_a)%f(:)
           !  particle_set(index_a)%f(:) = 0.0_dp
              
             wb = vsite_list(iconst)%wbc
             wd = (1.0d0-wb)/2.0d0
             particle_set(index_b)%f(:) = particle_set(index_b)%f(:) + wd*particle_set(index_a)%f(:)
             particle_set(index_c)%f(:) = particle_set(index_c)%f(:) + wb*particle_set(index_a)%f(:)
             particle_set(index_d)%f(:) = particle_set(index_d)%f(:) + wd*particle_set(index_a)%f(:)
             particle_set(index_a)%f(:) = 0.0_dp
           else
             wb = vsite_list(iconst)%wbc
             wd = vsite_list(iconst)%wdc
             wc = 1.0_dp - vsite_list(iconst)%wbc - vsite_list(iconst)%wdc

             particle_set(index_b)%f(:) = particle_set(index_b)%f(:) + wb*particle_set(index_a)%f(:)
             particle_set(index_c)%f(:) = particle_set(index_c)%f(:) + wc*particle_set(index_a)%f(:)
             particle_set(index_d)%f(:) = particle_set(index_d)%f(:) + wd*particle_set(index_a)%f(:)
             particle_set(index_a)%f(:) = 0.0_dp
        endif
      END DO
   END SUBROUTINE force_vsite_ext

END MODULE constraint_vsite
