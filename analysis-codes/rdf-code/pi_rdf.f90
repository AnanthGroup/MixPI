program pi_rdf
implicit none

real(kind=8),parameter      :: pi=3.141592654

integer        ::i,j,k
TYPE atom_type
  integer              ::nbead
  character(len=2)     ::atom_label
  real(kind=8),allocatable,dimension(:,:)  ::x
END TYPE
TYPE molecule_type
  character(len=160)   ::psf_input
  integer              ::natoms, total_beads, num_repeat
  TYPE(atom_type),allocatable,dimension(:)      ::atoms
  logical              ::first = .true.
END TYPE
TYPE classical_type
  character(len=2)            ::atom_label
  real(kind=8),dimension(3)   ::x
END TYPE
TYPE(molecule_type),allocatable,dimension(:)   ::molecules
TYPE(classical_type),allocatable,dimension(:)  ::classical

integer                ::num_molecules
integer                ::imolc, iatom, ibead
integer                ::jmolc, jatom, jbead
integer                ::nratio, max_beads, num_i, num_j

integer                ::nclassical, iclassical, jclassical

character(len=100)     ::trajectory_file
integer                ::num_steps, skip_steps_start, skip_steps_between, istep
integer                ::icount_steps
real(kind=8),dimension(3,3)    ::h, hinv
real(kind=8)           ::deth
character(len=2),dimension(2)  ::atom_choice
integer                ::total_particles, num_lines
real(kind=8)           ::dist

real(kind=8)             ::xmax, delx, rr, coord_num
integer                  :: nbins, ibin
integer,allocatable,dimension(:)   ::g
real(kind=8),allocatable,dimension(:)    ::g_norm

       interface
        SUBROUTINE get_inv(HMAT,HMATI,DETH)
        REAL(kind=8), intent(out), dimension(:,:) :: HMATI
        REAL(kind=8), intent(in), dimension(:,:) :: HMAT
        REAL(kind=8), intent(out) ::  DETH
        end 
       end interface

open(101,FILE='rp_molecule_info.inp')
read(101,*) num_molecules
allocate(molecules(num_molecules))
molecules(:)%natoms = 0
molecules(:)%num_repeat = 0
molecules(:)%psf_input = ''

imolc = 0
max_beads = 0
do while (imolc < num_molecules)
  imolc = imolc + 1
  read(101,*) molecules(imolc)%natoms, molecules(imolc)%num_repeat, molecules(imolc)%psf_input
  allocate(molecules(imolc)%atoms(molecules(imolc)%natoms))
  do iatom = 1, molecules(imolc)%natoms
    read(101,*) molecules(imolc)%atoms(iatom)%nbead
    allocate(molecules(imolc)%atoms(iatom)%x(molecules(imolc)%atoms(iatom)%nbead,3))
    molecules(imolc)%atoms(iatom)%x = 0.0d0
    if (molecules(imolc)%atoms(iatom)%nbead > max_beads) then 
      max_beads = molecules(imolc)%atoms(iatom)%nbead
    endif
  enddo
  molecules(imolc)%total_beads = sum(molecules(imolc)%atoms(:)%nbead)
  do i = 2, molecules(imolc)%num_repeat
    molecules(imolc+i-1)%natoms = molecules(imolc)%natoms
    allocate(molecules(imolc+i-1)%atoms(molecules(imolc)%natoms))
    molecules(imolc+i-1)%atoms(:)%nbead = molecules(imolc)%atoms(:)%nbead
    do iatom = 1, molecules(imolc+i-1)%natoms
      allocate(molecules(imolc+i-1)%atoms(iatom)%x(molecules(imolc+i-1)%atoms(iatom)%nbead,3))
      molecules(imolc+i-1)%atoms(iatom)%x = 0.0d0
    enddo
    molecules(imolc+i-1)%total_beads = sum(molecules(imolc+i-1)%atoms(:)%nbead)
    molecules(imolc+i-1)%first = .false.
  enddo
  imolc = imolc+molecules(imolc)%num_repeat-1
enddo
close(101)

open(102,file='sim_details.inp')
read(102,*) trajectory_file
read(102,*) num_steps
read(102,*) skip_steps_start
read(102,*) skip_steps_between
read(102,*) h(1,1), h(2,2), h(3,3)    !a,b,c of unit cell
read(102,*) atom_choice(1), atom_choice(2)
read(102,*) delx, xmax
close(102)

write(*,*)' delx is ', delx
write(*,*)' xmax is ', xmax
nbins = int(xmax/delx)
write(*,*)' number of bins ', nbins
allocate(g(nbins))
g = 0

call get_inv(h,hinv,deth)



open(103,file=trim(adjustl(trajectory_file)))
do istep = 1, skip_steps_start
  read(103,*) total_particles
  write(*,*)' total_particles are ', total_particles, istep
  if (istep == 1) then
    if (total_particles /= sum(molecules(:)%total_beads)) then
	write(*,*) total_particles
	write(*,*) sum(molecules(:)%total_beads)
	write(*,*) 'have classical atoms', total_particles - sum(molecules(:)%total_beads)
	nclassical = total_particles - sum(molecules(:)%total_beads)
	allocate(classical(nclassical))
    else
        nclassical = 0
    endif
  endif
  read(103,*)
  do iatom = 1, total_particles
    read(103,*)
  enddo
enddo
icount_steps = 0
skip_steps_between = skip_steps_between+1
do istep = skip_steps_start+1, num_steps
  !write(*,*)' istep ', istep
  !write(*,*)'istep-skip_steps_start+1 ', istep-skip_steps_start-1
  !write(*,*)' skip_step_between ', skip_steps_between
  !write(*,*)' mod ', mod(istep-skip_steps_start-1,skip_steps_between)
  
  if (skip_steps_between == 0 .or. mod(istep-skip_steps_start-1,skip_steps_between)==0 .or. istep == num_steps) then
    icount_steps = icount_steps + 1
    write(*,*)' istep and icount ', istep, icount_steps
    read(103,*) total_particles
  !write(*,*)' step number ', istep
  !if (total_particles /= sum(molecules(:)%total_beads)) then
!	write(*,*) total_particles
!	write(*,*) sum(molecules(:)%total_beads)
!	write(*,*)' rp_ file and trajectory file do not match '
!	STOP
!  endif
    read(103,*)
    imolc = 0
    num_lines = 0
    do imolc = 1, num_molecules
      do iatom = 1, molecules(imolc)%natoms
        do ibead = 1, molecules(imolc)%atoms(iatom)%nbead
	  num_lines = num_lines+1
          read(103,*) molecules(imolc)%atoms(iatom)%atom_label, molecules(imolc)%atoms(iatom)%x(ibead,1:3)
	  !write(104,*) molecules(imolc)%atoms(iatom)%atom_label, molecules(imolc)%atoms(iatom)%x(ibead,1:3)
        enddo
      enddo
    enddo
    !write(*,*)' num lines ', num_lines
    do iclassical = 1,  nclassical
      num_lines = num_lines + 1
      read(103,*) classical(iclassical)%atom_label, classical(iclassical)%x(1:3)
    enddo
  !write(*,*)' num lines ', num_lines

    if (istep == skip_steps_start+1) then
      num_i = 0
      num_j = 0
      do imolc = 1, num_molecules
      do iatom = 1, molecules(imolc)%natoms
        if (trim(adjustl(molecules(imolc)%atoms(iatom)%atom_label)) == trim(adjustl(atom_choice(1))) ) then
	  num_i = num_i + 1
        endif
        if (trim(adjustl(molecules(imolc)%atoms(iatom)%atom_label)) == trim(adjustl(atom_choice(2))) ) then
	  num_j = num_j + 1
        endif
      enddo
      enddo
      do iclassical = 1, nclassical
        if (trim(adjustl(classical(iclassical)%atom_label)) == trim(adjustl(atom_choice(1))) ) then
          num_i = num_i + 1
        endif
        if (trim(adjustl(classical(iclassical)%atom_label)) == trim(adjustl(atom_choice(2))) ) then
          num_j = num_j + 1
        endif
      enddo
      num_i = num_i*max_beads
      num_j = num_j!*max_beads
    !if (trim(adjustl(atom_choice(1))) == trim(adjustl(atom_choice(2)))) then
    !  num_j = num_i - max_beads
    !endif
    endif
	


    do imolc = 1, num_molecules
      do iatom = 1, molecules(imolc)%natoms
        if (trim(adjustl(molecules(imolc)%atoms(iatom)%atom_label)) == trim(adjustl(atom_choice(1))) ) then
          do jmolc = 1, num_molecules
            do jatom = 1, molecules(jmolc)%natoms
              if (trim(adjustl(molecules(jmolc)%atoms(jatom)%atom_label)) == trim(adjustl(atom_choice(2))) ) then
	        if (imolc /= jmolc .OR. (imolc == jmolc .and. iatom /= jatom)) then    !NOT SAME MOLECULE and ATOM
	            !write(*,*)' these match ', imolc, iatom, jmolc, jatom
		  !write(*,*)' nbeads are ', molecules(imolc)%atoms(iatom)%nbead, molecules(jmolc)%atoms(jatom)%nbead
                  if (molecules(imolc)%atoms(iatom)%nbead >= molecules(jmolc)%atoms(jatom)%nbead) then
                    nratio = molecules(imolc)%atoms(iatom)%nbead/molecules(jmolc)%atoms(jatom)%nbead
                    !write(*,*)' nratio is ',nratio
	            do ibead = 1, molecules(imolc)%atoms(iatom)%nbead
	              if (nratio /= 1) then
	                jbead = int(ibead/nratio-1.0d-6) + 1
		      else
	                jbead = ibead
	              endif
		      !write(*,*)' bead ', ibead, ' with ', jbead
		      call calc_distance(molecules(imolc)%atoms(iatom)%x(ibead,1:3), &
		      molecules(jmolc)%atoms(jatom)%x(jbead,1:3),hinv,h,dist)
		      if (dist > xmax) then
	                write(*,*)' atom ', imolc, iatom, jmolc, jatom, 'excluded at distance', dist
			write(*,*)' may be a problem with the defined box or xmax value '
		      else
		        ibin=int(dist/delx+.5d0) + 1
		        g(ibin) = g(ibin)+1*(max_beads/molecules(imolc)%atoms(iatom)%nbead)
		      endif
		    enddo
		  else
		    nratio = molecules(jmolc)%atoms(jatom)%nbead/molecules(imolc)%atoms(iatom)%nbead
	            do jbead = 1, molecules(jmolc)%atoms(jatom)%nbead
		      if (nratio /= 1) then
		        ibead = int(jbead/nratio-1.0d-6) + 1
		      else
		        ibead = jbead
		      endif
			!write(*,*)' bead ', jbead, ' with ', ibead
		      call calc_distance(molecules(imolc)%atoms(iatom)%x(ibead,1:3), &
		      molecules(jmolc)%atoms(jatom)%x(jbead,1:3),hinv, h, dist)
		      if (dist > xmax) then
		        write(*,*)' atom ', imolc, iatom, jmolc, jatom, 'excluded at distance', dist
		      else
		        ibin=int(dist/delx+.5d0) + 1
		        g(ibin) = g(ibin)+1*(max_beads/molecules(jmolc)%atoms(jatom)%nbead)
		      endif	  
		    enddo
		  endif
	        endif
	      endif
	    enddo
          enddo
          do jclassical = 1, nclassical
	  if (trim(adjustl(classical(jclassical)%atom_label)) == trim(adjustl(atom_choice(2))) ) then
	    !write(*,*)' these match ', imolc, iatom, jclassical
	    do ibead = 1, molecules(imolc)%atoms(iatom)%nbead
	      call calc_distance(molecules(imolc)%atoms(iatom)%x(ibead,1:3), classical(jclassical)%x(1:3), hinv,h,dist)
	      if (dist > xmax) then
    	        write(*,*)' atom ', imolc, iatom, jclassical, 'excluded at distance', dist
	      else
	        ibin=int(dist/delx+.5d0) + 1
                g(ibin) = g(ibin)+1*(max_beads/molecules(imolc)%atoms(iatom)%nbead)	
	      endif
	    enddo
	  endif
	enddo
      endif
    enddo
  enddo
  do iclassical = 1, nclassical
    if (trim(adjustl(classical(iclassical)%atom_label)) == trim(adjustl(atom_choice(1))) ) then
        do jmolc = 1, num_molecules
          do jatom = 1, molecules(jmolc)%natoms
            if (trim(adjustl(molecules(jmolc)%atoms(jatom)%atom_label)) == trim(adjustl(atom_choice(2))) ) then
	      !write(*,*)' these match ', iclassical, jmolc, jatom
	      do jbead = 1, molecules(jmolc)%atoms(jatom)%nbead
	        call calc_distance(classical(iclassical)%x(1:3), molecules(jmolc)%atoms(jatom)%x(jbead,1:3), hinv,h,dist)
		if (dist > xmax) then
		  write(*,*)' atom ', iclassical, jmolc, jatom, 'excluded at distance', dist
		else
		  ibin=int(dist/delx+.5d0) + 1
                  g(ibin) = g(ibin)+1*(max_beads/molecules(jmolc)%atoms(jatom)%nbead)
	        endif	  
	      enddo
	    endif
	  enddo
        enddo
        do jclassical = 1, nclassical
	  if (trim(adjustl(classical(jclassical)%atom_label)) == trim(adjustl(atom_choice(2))) ) then
	    !write(*,*)' these match ', iclassical, jclassical
	    if (iclassical /= jclassical) then    !NOT SAME ATOM
	      call calc_distance(classical(iclassical)%x(1:3), classical(jclassical)%x(1:3),hinv,h,dist)
		if (dist > xmax) then
		  write(*,*)' atom ', iclassical, jclassical, 'excluded at distance', dist
		else
	      ibin=int(dist/delx+.5d0) + 1
              g(ibin) = g(ibin)+1*max_beads	
	        endif
	    endif
	  endif
	enddo
    endif
  enddo
    
  else      
    read(103,*) total_particles
  !write(*,*)' step number ', istep
  !if (total_particles /= sum(molecules(:)%total_beads)) then
!	write(*,*) total_particles
!	write(*,*) sum(molecules(:)%total_beads)
!	write(*,*)' rp_ file and trajectory file do not match '
!	STOP
!  endif
    read(103,*)
    imolc = 0
    num_lines = 0
    do imolc = 1, num_molecules
      do iatom = 1, molecules(imolc)%natoms
        do ibead = 1, molecules(imolc)%atoms(iatom)%nbead
	  num_lines = num_lines+1
          read(103,*) molecules(imolc)%atoms(iatom)%atom_label, molecules(imolc)%atoms(iatom)%x(ibead,1:3)
	  !write(104,*) molecules(imolc)%atoms(iatom)%atom_label, molecules(imolc)%atoms(iatom)%x(ibead,1:3)
        enddo
      enddo
    enddo
    !write(*,*)' num lines ', num_lines
    do iclassical = 1,  nclassical
      num_lines = num_lines + 1
      read(103,*) classical(iclassical)%atom_label, classical(iclassical)%x(1:3)
    enddo
  endif
enddo



close(103)

allocate(g_norm(nbins))
g_norm = 0
write(*,*)' deth ', deth
write(*,*)' delx ', delx
write(*,*)' max_beads ', max_beads
write(*,*)' num_i and num_j ', num_i, num_j
write(*,*)' total number of steps ', icount_steps
do i = 1, nbins
  rr = ((dble(real(i))-0.50d0)*delx)**2
  g_norm(i) = deth*g(i)/(4.0d0*pi*rr*delx*dble(real(num_i*num_j)))
enddo
g_norm = g_norm/dble(real(icount_steps))

coord_num = 0.0d0
open(201,file='gofr.dat')
do i = 1, nbins
  rr = ((dble(real(i))-0.50d0)*delx)**2
  coord_num = coord_num + g(i)/(dble(real(icount_steps))*max_beads) !g_norm(i)*rr*4.0*pi*delx*num_i*num_j/deth
  write(201,*) i-1, (i-.5)*delx, g(i), g_norm(i), coord_num
enddo
close(201)

end program pi_rdf


! ----------------------------------------------------- !
      SUBROUTINE get_inv(HMAT,HMATI,DETH)
! GETS INVERSE, HMATI, OF THE IPERD DIMENSIONAL MATRIX HMAT
! (STORED AS A 3 X 3)
      REAL(kind=8), intent(out), dimension(:,:) :: HMATI
      REAL(kind=8), intent(in), dimension(:,:) :: HMAT
      REAL(kind=8), intent(out) ::  DETH
!
     deth = HMAT(1,1)*(HMAT(2,2)*HMAT(3,3)-HMAT(2,3)*HMAT(3,2))&
       + HMAT(1,2)*(HMAT(2,3)*HMAT(3,1)-HMAT(2,1)*HMAT(3,3))&
       + HMAT(1,3)*(HMAT(2,1)*HMAT(3,2)-HMAT(2,2)*HMAT(3,1))
      HMATI(1,1) = (HMAT(2,2)*HMAT(3,3)-HMAT(2,3)*HMAT(3,2))/deth
      HMATI(2,2) = (HMAT(1,1)*HMAT(3,3)-HMAT(1,3)*HMAT(3,1))/deth
      HMATI(3,3) = (HMAT(1,1)*HMAT(2,2)-HMAT(1,2)*HMAT(2,1))/deth
      HMATI(1,2) = (HMAT(1,3)*HMAT(3,2)-HMAT(1,2)*HMAT(3,3))/deth
      HMATI(2,1) = (HMAT(3,1)*HMAT(2,3)-HMAT(2,1)*HMAT(3,3))/deth
      HMATI(1,3) = (HMAT(1,2)*HMAT(2,3)-HMAT(1,3)*HMAT(2,2))/deth
      HMATI(3,1) = (HMAT(2,1)*HMAT(3,2)-HMAT(3,1)*HMAT(2,2))/deth
      HMATI(2,3) = (HMAT(1,3)*HMAT(2,1)-HMAT(2,3)*HMAT(1,1))/deth
      HMATI(3,2) = (HMAT(3,1)*HMAT(1,2)-HMAT(3,2)*HMAT(1,1))/deth
      END subroutine get_inv



      subroutine calc_distance(r1,r2,hinv,h,dist)
      implicit none
      real(kind=8),dimension(3,3)   :: hinv,h
      real(kind=8),dimension(3)     :: r1,r2, rij, sij
      real(kind=8)                  :: dist

      !distance vector
      rij = r2 - r1
      !scale
      sij = matmul(hinv,rij)
      !minimum image
      sij = sij-dnint(sij)
      !rescale
      rij = matmul(h,sij)
      dist = dsqrt(dot_product(rij,rij))
      
      end subroutine calc_distance
