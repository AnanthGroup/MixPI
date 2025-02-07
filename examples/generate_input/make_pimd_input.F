PROGRAM make_pimd_input
implicit none

TYPE atom_type
  integer                            ::nbead
  real(kind=8)                       ::mass
  real(kind=8)                       ::charge
  character(len=160)                 ::ff_type
  character(len=2)                   ::element
END TYPE atom_type
TYPE molecule_type
  character(len=160)                          ::psf_input     
  integer                                     ::natoms, total_beads, num_repeat
  TYPE(atom_type),allocatable,dimension(:)    ::atoms
  logical                                     ::first=.true.
  integer                              :: molecule_num_bonds
  integer, allocatable,dimension(:)    :: molecule_bond_a, molecule_bond_b
  integer                              :: molecule_num_bends
  integer, allocatable,dimension(:)    :: molecule_bend_a, molecule_bend_b, molecule_bend_c
END TYPE molecule_type
TYPE bonded_ff
  character(len=2)                   ::atom1, atom2
  character(len=25)		     ::k_units, r_units
  real(kind=8)                       ::k, r0
END TYPE bonded_ff
TYPE bended_ff
  character(len=2)                   ::atom1, atom2, atom3
  character(len=25)		     ::k_units, theta_units
  real(kind=8)                       ::k, theta0
END TYPE bended_ff
TYPE ljed_ff
  character(len=2)                   ::atom
  real(kind=8)                       ::eps, sigma
  real(kind=8)                       :: rcut_value
  character(len=25)                  ::eps_units, sigma_units
END TYPE ljed_ff

integer                                         ::i,j,k
integer                                         ::ibead, iatom, jatom
integer                                         ::num_molecules, imolc, jmolc, total_atoms
TYPE(molecule_type),allocatable,dimension(:)    ::molecules

logical                              :: gen_xyz, gen_psf, gen_ff
logical				     :: all_classical_ff
character(len=160)                   :: atom_label, format_string, new_label
character(len=160)                   :: format_string_ff, format_string_id
real(kind=8),dimension(3)            :: x, com_rnd
real(kind=8),allocatable,dimension(:,:)    ::bead_rnd

character(len=160)                   :: psf_format
character(len=3)                     :: char_check
integer                              :: psf_natom, nskip
integer                              :: i_dum, mol_num
real(kind=8)                         :: atom_charge, atom_mass
character(len=160)                   :: mol_name, res_id, atom_id, atom_ff_type

integer                              :: num_bonds
integer, allocatable,dimension(:)    :: bond_a, bond_b
integer                              :: ibond
integer                              :: num_bends
integer, allocatable,dimension(:)    :: bend_a, bend_b, bend_c
integer                              :: ibend
integer                              :: icount, Na, Nb, Nc, nratio, abead, bbead, aratio, bratio
integer                              :: rp_num_bonds
integer, allocatable,dimension(:)    :: rp_bond_a, rp_bond_b
integer                              :: rp_num_bends
integer, allocatable,dimension(:)    :: rp_bend_a, rp_bend_b, rp_bend_c
integer                              :: num_phi, num_imphi, num_don, num_acc, num_nb


TYPE(bonded_ff),allocatable,dimension(:)     ::bond_ff
TYPE(bended_ff),allocatable,dimension(:)     ::bend_ff
TYPE(ljed_ff), allocatable,dimension(:)      ::lj_ff
integer                                      ::bond_num, bend_num, i_ff
character(len=4)			     ::ff_type

integer                                      ::imatch, jmatch, icheck, lj_num
character(len=6),allocatable,dimension(:,:)    ::check_array
character(len=6),allocatable,dimension(:)      ::label_check
logical                                        ::check_index
integer                                        ::label_i, label_j, label_k
integer                                        ::start_index

real(kind=8),dimension(3)                    ::rnd
real(kind=8)                                 ::pi_radius

logical                                      ::print_extra

gen_xyz = .true.
gen_psf = .true.
gen_ff = .true.
all_classical_ff = .true.

print_extra = .false.

pi_radius = 0.25

write(*,*)' program to make the psf, xyz, and ff files for the relevant RP types '

open(111,file='rp_molecule_info.inp')
read(111,*) num_molecules
allocate(molecules(num_molecules))
imolc = 0
do while (imolc < num_molecules)
  imolc = imolc + 1
  read(111,*) molecules(imolc)%natoms, molecules(imolc)%num_repeat, molecules(imolc)%psf_input
  allocate(molecules(imolc)%atoms(molecules(imolc)%natoms))
  do i = 1, molecules(imolc)%natoms
    read(111,*) molecules(imolc)%atoms(i)%nbead
  enddo
  molecules(imolc)%total_beads = sum(molecules(imolc)%atoms(:)%nbead)
  do i = 2, molecules(imolc)%num_repeat
    molecules(imolc+i-1)%natoms = molecules(imolc)%natoms
    allocate(molecules(imolc+i-1)%atoms(molecules(imolc)%natoms))
    molecules(imolc+i-1)%atoms(:)%nbead = molecules(imolc)%atoms(:)%nbead
    molecules(imolc+i-1)%first = .false.
    molecules(imolc+i-1)%total_beads = molecules(imolc)%total_beads
    molecules(imolc+i-1)%psf_input = molecules(imolc)%psf_input
  enddo
  imolc = imolc+molecules(imolc)%num_repeat-1
enddo
close(111)


! --------------------------------------------------------------------------------- !
! CREATES THE XYZ FILE------------------------------------------------------------- !
call random_seed()
if (gen_xyz) then
  open(111, file='structure.xyz')
  open(333, file='RP_structure.xyz')
  read(111,*) total_atoms
  read(111,*)
  write(333,*) total_atoms - sum(molecules(:)%natoms) + sum(molecules(:)%total_beads)
  write(333,*)

  do imolc = 1, SIZE(molecules)
    do iatom = 1, molecules(imolc)%natoms
      read(111,*) atom_label, x(:) 
      if (molecules(imolc)%first) then
        molecules(imolc)%atoms(iatom)%element = trim(adjustl(atom_label))
      endif
      allocate(bead_rnd(molecules(imolc)%atoms(iatom)%nbead,3))
      do ibead = 1, molecules(imolc)%atoms(iatom)%nbead
        call random_number(rnd(1))
        call random_number(rnd(2))
        call random_number(rnd(3))
        rnd(:) = (rnd(:)-0.5d0) * 2.0d0 * pi_radius
        bead_rnd(ibead,1:3) = rnd(:)
      enddo
      com_rnd(1) = sum(bead_rnd(:,1))
      com_rnd(2) = sum(bead_rnd(:,2))
      com_rnd(3) = sum(bead_rnd(:,3))
      com_rnd(:) = com_rnd(:)/dble(real(molecules(imolc)%atoms(iatom)%nbead))
      bead_rnd(:,1) = bead_rnd(:,1)-com_rnd(1)
      bead_rnd(:,2) = bead_rnd(:,2)-com_rnd(2)
      bead_rnd(:,3) = bead_rnd(:,3)-com_rnd(3)
      do ibead = 1, molecules(imolc)%atoms(iatom)%nbead
        if (len(trim(atom_label)) == 1) then
          if (ibead < 10) format_string = "(A1,I1)"
          if (ibead < 100 .and. ibead >= 10) format_string = "(A1,I2)"
          if (ibead < 1000 .and. ibead >= 100) format_string = "(A1,I3)"
          if (ibead < 10000 .and. ibead >= 1000) format_string = "(A1,I4)"
          if (ibead >= 10000) then
            write(*,*) 'ERROR-NO INDEX LARGE ENOUGH' 
            EXIT
          endif
        else 
          if (ibead < 10) format_string = "(A2,I1)"
          if (ibead < 100 .and. ibead >= 10) format_string = "(A2,I2)"
          if (ibead < 1000 .and. ibead >= 100) format_string = "(A2,I3)"
          if (ibead < 10000 .and. ibead >= 1000) format_string = "(A2,I4)"
          if (ibead >= 10000) then 
            write(*,*) 'ERROR-NO INDEX LARGE ENOUGH' 
            EXIT
          endif
        endif
        write(new_label,format_string) atom_label, ibead
        write(333,*) trim(new_label), x(:) + bead_rnd(ibead,:)
      enddo
      deallocate(bead_rnd)
    enddo
  enddo
  do i = sum(molecules(:)%natoms)+1, total_atoms
      read(111,*) atom_label, x(:) 
      write(333,*) trim(atom_label), x(:)
  enddo

  close(111)
  close(333)
endif
! ------------------------------------------------------------------------------- !

open(111,file='RP_molset.pot')
icount = 0
do imolc = 1, SIZE(molecules)
  if (molecules(imolc)%first) then
    if (imolc < 10) write(new_label,"(A5)") 'rp-'//char(imolc+48)//'-'
    if (imolc >= 10) then
      write(*,*)' not currently supporting more then 9 psf RP at a time '
      EXIT
    endif
     
    write(111,*)'&MOLECULE'
    write(111,*)'  NMOL ', molecules(imolc)%num_repeat
    write(111,*)'  CONNECTIVITY PSF' 
    write(111,*)'  CONN_FILE_NAME ', trim(new_label)//trim(molecules(imolc)%psf_input)
    write(111,*)'&END MOLECULE'
  endif
enddo
close(111)
  

! ------------------------------------------------------------------------------- !
! GENERATES THE PSF FILE -------------------------------------------------------- !
if (gen_psf) then
  psf_format = '(I10,T12,A7,T21,I8,T30,A7,T39,A6,T47,A6,T53,F10.6,T69,F12.8,T88,I1)'
  do imolc = 1, SIZE(molecules)
    if (molecules(imolc)%first) then
      open(111, file=trim(molecules(imolc)%psf_input))
      if (imolc < 10) write(new_label,"(A5)") 'rp-'//char(imolc+48)//'-'
      if (imolc >= 10) then
        write(*,*)' not currently supporting more then 9 psf RP at a time '
        EXIT
      endif
      open(222, file=trim(new_label)//trim(molecules(imolc)%psf_input)) 
      read(111, *) char_check
      write(222,*) "PSF EXT"
      !skip line after ext
      read(111,*)
      write(222,*)
      !skip number of title lines indicated
      read(111,*) nskip
      write(222,*) nskip, "!NTITLE"
      do i = 1, nskip
        read(111,*)
        write(222,*)
      enddo
      read(111,*)
      write(222,*)
      read(111,*) psf_natom
      if (psf_natom /= molecules(imolc)%natoms) then
        write(*,*)' psf file does not match atom number of input file '
        EXIT
      endif
      write(222,*) SUM(molecules(imolc)%atoms(:)%nbead), '!NATOMS '


      if (molecules(imolc)%total_beads < 10) format_string = "(A2,I1)"
      if (molecules(imolc)%total_beads < 100 .and. molecules(imolc)%total_beads >= 10) format_string = "(A2,I2)"
      if (molecules(imolc)%total_beads < 1000 .and. molecules(imolc)%total_beads >= 100) format_string = "(A2,I3)"
      if (molecules(imolc)%total_beads < 10000 .and. molecules(imolc)%total_beads >= 1000) format_string = "(A2,I4)"
      if (molecules(imolc)%total_beads >= 10000) then 
        write(*,*) 'ERROR-NO INDEX LARGE ENOUGH' 
        EXIT
      endif
      do iatom = 1, molecules(imolc)%natoms
        read(111,psf_format) i_dum, mol_name, mol_num, res_id, atom_id, atom_ff_type, atom_charge, atom_mass
	write(mol_name,format_string) 'RP',molecules(imolc)%total_beads

	molecules(imolc)%atoms(iatom)%mass = atom_mass/real(molecules(imolc)%atoms(iatom)%nbead)
	molecules(imolc)%atoms(iatom)%charge = atom_charge/real(molecules(imolc)%atoms(iatom)%nbead)

        if (molecules(imolc)%atoms(iatom)%nbead < 10) format_string_ff = "(I1,A2)"
        if (molecules(imolc)%atoms(iatom)%nbead < 100 .and. molecules(imolc)%atoms(iatom)%nbead >= 10) &
		format_string_ff = "(I2,A2)"
        if (molecules(imolc)%atoms(iatom)%nbead < 1000 .and. molecules(imolc)%atoms(iatom)%nbead >= 100) &
		format_string_ff = "(I3,A2)"
        if (molecules(imolc)%atoms(iatom)%nbead < 10000 .and. molecules(imolc)%atoms(iatom)%nbead >= 1000) &
		format_string_ff = "(I4,A2)"
        if (molecules(imolc)%atoms(iatom)%nbead >= 10000) then 
          write(*,*) 'ERROR-NO INDEX LARGE ENOUGH' 
          EXIT
        endif

	atom_ff_type = trim(adjustl(atom_ff_type))
	write(molecules(imolc)%atoms(iatom)%ff_type,format_string_ff) molecules(imolc)%atoms(iatom)%nbead, atom_ff_type(1:2)
	do i = 1, molecules(imolc)%atoms(iatom)%nbead
          if (i < 10) format_string_id = "(A2,I1)"
          if (i < 100 .and. i >= 10) format_string_id = "(A2,I2)"
          if (i < 1000 .and. i >= 100) format_string_id = "(A2,I3)"
          if (i < 10000 .and. i >= 1000) format_string_id = "(A2,I4)"
          if (i >= 10000) then 
            write(*,*) 'ERROR-NO INDEX LARGE ENOUGH' 
            EXIT
          endif
	  write(atom_id,format_string_id) atom_ff_type(1:2),i
          write(222,psf_format) sum(molecules(imolc)%atoms(1:iatom-1)%nbead)+i, &
		  mol_name, imolc, 'RP', atom_id, molecules(imolc)%atoms(iatom)%ff_type, &
		  molecules(imolc)%atoms(iatom)%charge, molecules(imolc)%atoms(iatom)%mass, 1
        enddo

      enddo

      !assume 2 blank lines between sections
      read(111,*)
      read(111,*)
      write(222,*)
      write(222,*)

      read(111,*) molecules(imolc)%molecule_num_bonds
      num_bonds = molecules(imolc)%molecule_num_bonds
      if (num_bonds > 0) then
	allocate(bond_a(num_bonds), bond_b(num_bonds))
	allocate(molecules(imolc)%molecule_bond_a(num_bonds))
	allocate(molecules(imolc)%molecule_bond_b(num_bonds))
        do i = 1, num_bonds, 4
	  read(111,*) (bond_a(i + j-1), bond_b(i+j-1), j=1, min(num_bonds-i+1, 4)) 
	enddo

	molecules(imolc)%molecule_bond_a = bond_a
	molecules(imolc)%molecule_bond_b = bond_b

	!allocate(bond_list(molecules(imolc)%natoms))
	!do i = 1, molecules(imolc)%natoms
	!  allocate(bond_list(i)%array(0))
	!enddo
	!CALL reorder_bonds(bond_list, bond_a, bond_b, num_bonds) 

	write(*,*)' imolc is ', imolc
	if (print_extra) then 
          write(*,*)' with bond list '
          do ibond = 1, num_bonds
            write(*,*) bond_a(ibond), bond_b(ibond)
          enddo
        endif

	write(*,*)' num_bonds', num_bonds
	icount = 0
	do ibond = 1, num_bonds
	  Na = molecules(imolc)%atoms(bond_a(ibond))%nbead
	  Nb = molecules(imolc)%atoms(bond_b(ibond))%nbead
	  if (print_extra) write(*,*)' na and nb ', Na, Nb
	  if (Na < Nb) then
            nratio = Nb/Na
	    if (print_extra) write(*,*)' Na is smaller so ratio is ', nratio
            do abead = 1, Na
              do bratio = 1, nratio
	        icount = icount+1
	      enddo
	    enddo
	  else
            nratio = Na/Nb
	    if (print_extra) write(*,*)' Nb is smaller so ratio is ', nratio
	    do bbead = 1, Nb
              do aratio = 1, nratio
	        icount = icount+1
	      enddo
            enddo
          endif
        enddo

	if (print_extra) write(*,*)' number of RP new bonds are ', icount

        rp_num_bonds = icount
        allocate(rp_bond_a(rp_num_bonds), rp_bond_b(rp_num_bonds))

	icount = 0
	do ibond = 1, num_bonds
	  Na = molecules(imolc)%atoms(bond_a(ibond))%nbead
	  Nb = molecules(imolc)%atoms(bond_b(ibond))%nbead
	  if (Na < Nb) then
            nratio = Nb/Na
            do abead = 1, Na
              do bratio = 1, nratio
	        icount = icount+1
		if (bond_a(ibond) == 1) then
		  rp_bond_a(icount) = abead
		else
		  rp_bond_a(icount) = abead+sum(molecules(imolc)%atoms(1:bond_a(ibond)-1)%nbead)
		endif
		if (bond_b(ibond) == 1) then
		  rp_bond_b(icount) = ((abead-1)*nratio+bratio)
		else
		  rp_bond_b(icount) = ((abead-1)*nratio+bratio)+sum(molecules(imolc)%atoms(1:bond_b(ibond)-1)%nbead)
		endif
	      enddo
	    enddo
	  else
            nratio = Na/Nb
	    do bbead = 1, Nb
              do aratio = 1, nratio
	        icount = icount+1
		if (bond_a(ibond) == 1) then
		  rp_bond_a(icount) = ((bbead-1)*nratio+aratio)
		else
		  rp_bond_a(icount) = ((bbead-1)*nratio+aratio)+sum(molecules(imolc)%atoms(1:bond_a(ibond)-1)%nbead)
		endif
		if (bond_b(ibond) == 1) then
		  rp_bond_b(icount) = bbead
		else
		  rp_bond_b(icount) = bbead+sum(molecules(imolc)%atoms(1:bond_b(ibond)-1)%nbead)
		endif
	      enddo
            enddo
          endif
        enddo
        if (print_extra) then  
	  write(*,*)' with rp_ bond list '
	  do ibond = 1, rp_num_bonds
	    write(*,*) rp_bond_a(ibond), rp_bond_b(ibond)
	  enddo
        endif

	write(222,*) rp_num_bonds, '!NBONDS'
	do i = 1, rp_num_bonds, 4
	  write(222,'(8I10)') (rp_bond_a(i+j-1), rp_bond_b(i+j-1),j=1,min(rp_num_bonds-i+1,4))
	enddo

	!allocate(bond_list(molecules(imolc)%natoms))

	deallocate(bond_a, bond_b)
	deallocate(rp_bond_a, rp_bond_b)

      else
        write(222,*) 0, '!NBOND'

      endif	!END OF BOND IF STATEMENT

      !assume 2 blank lines between sections
      read(111,*)
      read(111,*)
      write(222,*)
      write(222,*)

      read(111,*) molecules(imolc)%molecule_num_bends
      num_bends = molecules(imolc)%molecule_num_bends
      if (num_bends > 0) then

	allocate(bend_a(num_bends), bend_b(num_bends), bend_c(num_bends))
	allocate(molecules(imolc)%molecule_bend_a(num_bends))
	allocate(molecules(imolc)%molecule_bend_b(num_bends))
	allocate(molecules(imolc)%molecule_bend_c(num_bends))
        do i = 1, num_bonds, 3
	  read(111,*) (bend_a(i+j-1), bend_b(i+j-1), bend_c(i+j-1), j=1, min(num_bends-i+1, 3)) 
	enddo
	molecules(imolc)%molecule_bend_a = bend_a
	molecules(imolc)%molecule_bend_b = bend_b
	molecules(imolc)%molecule_bend_c = bend_c
	if (print_extra) then
          write(*,*)' with bend list '
	  do ibend = 1, num_bends
	    write(*,*) bend_a(ibend), bend_b(ibend), bend_c(ibend)
	  enddo
        endif

	icount = 0
	do ibend = 1, num_bends
	  Na = molecules(imolc)%atoms(bend_a(ibend))%nbead
	  Nb = molecules(imolc)%atoms(bend_b(ibend))%nbead
	  Nc = molecules(imolc)%atoms(bend_c(ibend))%nbead

	  if (Na >= Nb .and. Na >= Nc) then
              do abead = 1, Na
	        icount = icount+1
		
	      enddo
	  else if (Nb >= Na .and. Nb >= Nc) then
              do abead = 1, Nb
	        icount = icount+1
	      enddo
      else if (Nc >= Na .and. Nc >= Nb) then
              do abead = 1, Nc
	        icount = icount+1
	      enddo
	  endif

        enddo	  

	if (print_extra) write(*,*)' num of rp_bends is ', icount
	rp_num_bends = icount
	allocate(rp_bend_a(rp_num_bends), rp_bend_b(rp_num_bends), rp_bend_c(rp_num_bends))

	icount = 0
	do ibend = 1, num_bends
	  Na = molecules(imolc)%atoms(bend_a(ibend))%nbead
	  Nb = molecules(imolc)%atoms(bend_b(ibend))%nbead
	  Nc = molecules(imolc)%atoms(bend_c(ibend))%nbead

	  if (Na >= Nb .and. Na >= Nc) then
              do abead = 1, Na
		icount = icount+1
		if (bend_a(ibend) == 1) then
		  rp_bend_a(icount) = abead
		else
		  rp_bend_a(icount) = abead+sum(molecules(imolc)%atoms(1:bend_a(ibend)-1)%nbead)
		endif
		aratio = ((abead-1)*Nb)/Na+1
		if (bend_b(ibend) == 1) then
	          rp_bend_b(icount) = aratio
		else
	          rp_bend_b(icount) = aratio+sum(molecules(imolc)%atoms(1:bend_b(ibend)-1)%nbead)
		endif
		bratio = ((abead-1)*Nc)/Na+1
		if (bend_c(ibend) == 1) then
	          rp_bend_c(icount) = bratio
		else
	          rp_bend_c(icount) = bratio+sum(molecules(imolc)%atoms(1:bend_c(ibend)-1)%nbead)
		endif
	      enddo
          else if (Nb >= Na .and. Nb >= Nc) then
              do abead = 1, Nb
		icount = icount+1
		if (bend_b(ibend) == 1) then
		  rp_bend_b(icount) = abead
		else
		  rp_bend_b(icount) = abead+sum(molecules(imolc)%atoms(1:bend_b(ibend)-1)%nbead)
		endif
		aratio = ((abead-1)*Na)/Nb+1
		if (bend_a(ibend) == 1) then
	          rp_bend_a(icount) = aratio
		else
	          rp_bend_a(icount) = aratio+sum(molecules(imolc)%atoms(1:bend_a(ibend)-1)%nbead)
		endif
		bratio = ((abead-1)*Nc)/Nb+1
		if (bend_c(ibend) == 1) then
	          rp_bend_c(icount) = bratio
		else
	          rp_bend_c(icount) = bratio+sum(molecules(imolc)%atoms(1:bend_c(ibend)-1)%nbead)
		endif
	      enddo
	   else if (Nc >= Na .and. Nc >= Nb) then
              do abead = 1, Nc
		icount = icount+1
		if (bend_c(ibend) == 1) then
		  rp_bend_c(icount) = abead
		else
		  rp_bend_c(icount) = abead+sum(molecules(imolc)%atoms(1:bend_c(ibend)-1)%nbead)
		endif
		aratio = ((abead-1)*Na)/Nc+1
		if (bend_a(ibend) == 1) then
	          rp_bend_a(icount) = aratio
		else
	          rp_bend_a(icount) = aratio+sum(molecules(imolc)%atoms(1:bend_a(ibend)-1)%nbead)
		endif
		bratio = ((abead-1)*Nb)/Nc+1
		if (bend_b(ibend) == 1) then
	          rp_bend_b(icount) = bratio
		else
	          rp_bend_b(icount) = bratio+sum(molecules(imolc)%atoms(1:bend_b(ibend)-1)%nbead)
		endif
	      enddo
	   else
	      WRITE(*,*)' SOMETHING HAS GONE WRONG IN BEAD LOOP FOR BENDS '
	      EXIT
           endif
        enddo
        if (print_extra) then
	  write(*,*)' with rp_ bend list '
	  do ibend = 1, rp_num_bends
	    write(*,*) rp_bend_a(ibend), rp_bend_b(ibend), rp_bend_c(ibend)
	  enddo
        endif

	write(222,*) rp_num_bends, '!NTHETA'
	do i = 1, rp_num_bends, 3
	  write(222,'(9I10)') (rp_bend_a(i+j-1), rp_bend_b(i+j-1), rp_bend_c(i+j-1),j=1,min(rp_num_bends-i+1,3))
	enddo

	deallocate(bend_a, bend_b, bend_c)
	deallocate(rp_bend_a, rp_bend_b, rp_bend_c)
      else
        write(222,*) 0, '!NTHETA'


      endif

      !assume 2 blank lines between sections
      read(111,*)
      read(111,*)
      write(222,*)
      write(222,*)

      read(111,*) num_phi
      write(222,*) num_phi, '!NPHI'
      if (num_phi > 0) WRITE(*,*)' NPHI (and beyond) is not currently supported by this code '

      !assume 2 blank lines between sections
      read(111,*)
      read(111,*)
      write(222,*)
      write(222,*)

      read(111,*) num_imphi
      write(222,*) num_imphi, '!NIMPHI'

      !assume 2 blank lines between sections
      read(111,*)
      read(111,*)
      write(222,*)
      write(222,*)

      read(111,*) num_don
      write(222,*) num_don, '!NDON'

      !assume 2 blank lines between sections
      read(111,*)
      read(111,*)
      write(222,*)
      write(222,*)

      read(111,*) num_acc
      write(222,*) num_acc, '!NACC'

      !assume 2 blank lines between sections
      read(111,*)
      read(111,*)
      write(222,*)
      write(222,*)

      read(111,*) num_nb
      write(222,*) num_nb, '!NNB'

      !assume 2 blank lines between sections
      read(111,*)
      read(111,*)
      write(222,*)
      write(222,*)



      close(222)
      close(111)
      
    endif
  enddo
endif
! ------------------------------------------------------------------ !

open(111,file='RP_type.pot')
icount = 0
allocate(label_check(total_atoms))
do imolc = 1, SIZE(molecules)
  if (molecules(imolc)%first) then
    do iatom = 1, molecules(imolc)%natoms
      icount = 0
      do i = 1, icount
        if (label_check(i) == trim(molecules(imolc)%atoms(iatom)%ff_type)) then
          icount = -1
	  EXIT
        endif
      enddo
      if (icount == 0) then
        write(111,*)'&KIND ', trim(molecules(imolc)%atoms(iatom)%ff_type)
	write(111,*)'  ELEMENT ', trim(molecules(imolc)%atoms(iatom)%element)
	write(111,*)'  MASS ', molecules(imolc)%atoms(iatom)%mass
	write(111,*)'&END KIND'
      endif
    enddo
  endif
enddo
close(111)
      
    





! ------------------------------------------------------------------ !
! GENERATES THE FF INFORMATION (CURRENT VERSION ASSUMES A CP2K ----- !
! STYLE INPUT FILE CALLED FF_FILE.POT THAT CONTAINS BOND, BEND, ---- !
! TORSION, AND NONBONDED (LJ) STYLE INTERACTIONS. MORE CAN BE ADDED- !
if (gen_ff .and. gen_psf) then
  open(222, file='RP_FF_bond.pot')
  open(111,file='FF_info.inp')

  if (print_extra) write(*,*) ' inside loop of gen_ff and gen_psf'
  read(111,*)
  read(111,*) ff_type, bond_num
  allocate(bond_ff(bond_num))
  do i = 1, bond_num
    read(111,*) bond_ff(i)%atom1, bond_ff(i)%atom2, bond_ff(i)%k, bond_ff(i)%k_units, bond_ff(i)%r0, bond_ff(i)%r_units
  enddo

  do imolc = 1, SIZE(molecules) 
    if (print_extra) write(*,*)' imolc is ', imolc
    if (molecules(imolc)%first .and. molecules(imolc)%molecule_num_bonds > 0) then
    allocate(bond_a(molecules(imolc)%molecule_num_bonds), bond_b(molecules(imolc)%molecule_num_bonds))
    bond_a = molecules(imolc)%molecule_bond_a
    bond_b = molecules(imolc)%molecule_bond_b
    do ibond = 1, SIZE(bond_a)
      if (print_extra) write(*,*)' ibond is ', ibond
      do i_ff = 1, bond_num   !loop over number of bond types defined in FF file
	      if (print_extra) write(*,*)' molecueles ', trim(molecules(imolc)%atoms(bond_a(ibond))%ff_type), bond_ff(i_ff)%atom1

	label_i = LEN(trim(molecules(imolc)%atoms(bond_a(ibond))%ff_type))
	label_j = LEN(trim(molecules(imolc)%atoms(bond_b(ibond))%ff_type))
        if (trim(molecules(imolc)%atoms(bond_a(ibond))%ff_type(label_i-1:)) == bond_ff(i_ff)%atom1 .and. &
	      trim(molecules(imolc)%atoms(bond_b(ibond))%ff_type(label_j-1:)) == bond_ff(i_ff)%atom2) then
                write(222,*)'&BOND'
                write(222,*)'  ATOMS ', trim(molecules(imolc)%atoms(bond_a(ibond))%ff_type), '      ', &
			trim(molecules(imolc)%atoms(bond_b(ibond))%ff_type)
                write(222,*)'  K     ', trim(bond_ff(i_ff)%k_units), bond_ff(i_ff)%k/&
			real(max(molecules(imolc)%atoms(bond_a(ibond))%nbead, &
			molecules(imolc)%atoms(bond_b(ibond))%nbead))
                write(222,*)'  R0    ', bond_ff(i_ff)%r_units, bond_ff(i_ff)%r0
		write(222,*)'  KIND G87 '
                write(222,*)'&END BOND'
		EXIT
	else if (trim(molecules(imolc)%atoms(bond_b(ibond))%ff_type(label_i-1:)) == bond_ff(i_ff)%atom1 .and. &
	      trim(molecules(imolc)%atoms(bond_a(ibond))%ff_type(label_j-1:)) == bond_ff(i_ff)%atom2) then
                write(222,*)'&BOND'
                write(222,*)'  ATOMS ', trim(molecules(imolc)%atoms(bond_a(ibond))%ff_type), '      ', &
			trim(molecules(imolc)%atoms(bond_b(ibond))%ff_type)
                write(222,*)'  K     ', trim(bond_ff(i_ff)%k_units), bond_ff(i_ff)%k/&
			real(max(molecules(imolc)%atoms(bond_a(ibond))%nbead, &
			molecules(imolc)%atoms(bond_b(ibond))%nbead))
                write(222,*)'  R0    ', bond_ff(i_ff)%r_units, bond_ff(i_ff)%r0
		write(222,*)'  KIND G87 '
                write(222,*)'&END BOND'
		EXIT
	endif
      enddo
    enddo
    deallocate(bond_a, bond_b)
    endif
  enddo

  read(111,*)
  read(111,*)
  read(111,*) ff_type, bend_num
  allocate(bend_ff(bend_num))
  do i = 1, bend_num
    read(111,*) bend_ff(i)%atom1, bend_ff(i)%atom2, bend_ff(i)%atom3, bend_ff(i)%k, bend_ff(i)%k_units, bend_ff(i)%theta0, &
	    bend_ff(i)%theta_units
  enddo

  do imolc = 1, SIZE(molecules)
    if (molecules(imolc)%first .and. molecules(imolc)%molecule_num_bends > 0) then
      allocate(bend_a(molecules(imolc)%molecule_num_bends))
      allocate(bend_b(molecules(imolc)%molecule_num_bends))
      allocate(bend_c(molecules(imolc)%molecule_num_bends))
      bend_a = molecules(imolc)%molecule_bend_a
      bend_b = molecules(imolc)%molecule_bend_b
      bend_c = molecules(imolc)%molecule_bend_c
      do ibend = 1, SIZE(bend_a)
        do i_ff = 1, bend_num
          label_i = len(trim(molecules(imolc)%atoms(bend_a(ibend))%ff_type))
          label_j = len(trim(molecules(imolc)%atoms(bend_b(ibend))%ff_type))
          label_k = len(trim(molecules(imolc)%atoms(bend_c(ibend))%ff_type))
	  if (trim(molecules(imolc)%atoms(bend_a(ibend))%ff_type(label_i-1:)) == bend_ff(i_ff)%atom1 .or. &
		  trim(molecules(imolc)%atoms(bend_a(ibend))%ff_type(label_i-1:)) == bend_ff(i_ff)%atom2 .or. &
		  trim(molecules(imolc)%atoms(bend_a(ibend))%ff_type(label_i-1:)) == bend_ff(i_ff)%atom3 ) then
	    if (trim(molecules(imolc)%atoms(bend_b(ibend))%ff_type(label_j-1:)) == bend_ff(i_ff)%atom1 .or. &
		  trim(molecules(imolc)%atoms(bend_b(ibend))%ff_type(label_j-1:)) == bend_ff(i_ff)%atom2 .or. &
		  trim(molecules(imolc)%atoms(bend_b(ibend))%ff_type(label_j-1:)) == bend_ff(i_ff)%atom3 ) then
	      if (trim(molecules(imolc)%atoms(bend_c(ibend))%ff_type(label_k-1:)) == bend_ff(i_ff)%atom1 .or. &
		  trim(molecules(imolc)%atoms(bend_c(ibend))%ff_type(label_k-1:)) == bend_ff(i_ff)%atom2 .or. &
		  trim(molecules(imolc)%atoms(bend_c(ibend))%ff_type(label_k-1:)) == bend_ff(i_ff)%atom3 ) then
                write(222,*)'&BEND'
                write(222,*)'  ATOMS ', trim(molecules(imolc)%atoms(bend_a(ibend))%ff_type), '      ', &
			trim(molecules(imolc)%atoms(bend_b(ibend))%ff_type), '     ', &
			trim(molecules(imolc)%atoms(bend_c(ibend))%ff_type)
                write(222,*)'  K     ', trim(bend_ff(i_ff)%k_units), bend_ff(i_ff)%k/&
			real(max(molecules(imolc)%atoms(bend_a(ibend))%nbead, &
			molecules(imolc)%atoms(bend_b(ibend))%nbead, &
			molecules(imolc)%atoms(bend_c(ibend))%nbead))
                write(222,*)'  THETA0    ', bend_ff(i_ff)%theta_units, bend_ff(i_ff)%theta0
		write(222,*)'  KIND G87 '
                write(222,*)'&END BEND'
		EXIT
	      endif
	    endif
	  endif
	enddo
      enddo
    endif
  enddo

  close(222)
  open(222, file='RP_FF_nonbond.pot')

  read(111,*)
  read(111,*)
  read(111,*) ff_type, lj_num
  if (print_extra) write(*,*)'ff_type and lj_num ', ff_type, lj_num
  allocate(lj_ff(lj_num))
  do i = 1, lj_num
    read(111,*) lj_ff(i)%atom, lj_ff(i)%eps, lj_ff(i)%eps_units, lj_ff(i)%sigma, &
	    lj_ff(i)%sigma_units, lj_ff(i)%rcut_value
  enddo 

  allocate(check_array(100,2))
  check_array(:,:) = ' '

  icount = 0
  do imolc = 1, SIZE(molecules)
    if(molecules(imolc)%first) then
    do iatom = 1, molecules(imolc)%natoms
      imatch = 0
      do i_ff = 1, lj_num
	label_i = LEN(trim(molecules(imolc)%atoms(iatom)%ff_type))
        if (lj_ff(i_ff)%atom == molecules(imolc)%atoms(iatom)%ff_type(label_i-1:)) then
	  imatch = i_ff
	  if (print_extra) write(*,*)' imatch '
	  EXIT
	endif
      enddo
      do jatom = 1, molecules(imolc)%natoms   !LOOP FOR BEADS/ATOMS IN SAME MOLECULE
	jmatch = 0
        do i_ff = 1, lj_num
	  label_j = LEN(trim(molecules(imolc)%atoms(jatom)%ff_type))
	  if (lj_ff(i_ff)%atom == molecules(imolc)%atoms(jatom)%ff_type(label_j-1:)) then
	    jmatch = i_ff
	    if (print_extra) write(*,*)' jmatch '
	    EXIT
	  endif
	enddo
	if (print_extra) write(*,*)'outside array with ', imatch, jmatch
	if (imatch == 0 .or. jmatch == 0) then 
          write(*,*) 'combination does not have LJ input'
	else
          check_index = .true.
          do icheck = 1, icount
	    !write(*,*)' check_array '
            if (trim(molecules(imolc)%atoms(iatom)%ff_type) == trim(molecules(imolc)%atoms(jatom)%ff_type)) then
              if (trim(molecules(imolc)%atoms(iatom)%ff_type) == trim(check_array(icheck,1)) .and. &
			  trim(molecules(imolc)%atoms(jatom)%ff_type) == trim(check_array(icheck,2))) then
	        check_index = .false.
		EXIT
	      endif
            else if ((trim(molecules(imolc)%atoms(iatom)%ff_type) == check_array(icheck,1) .or. &
		    trim(molecules(imolc)%atoms(iatom)%ff_type) == check_array(icheck,2)) .and. &
		    (trim(molecules(imolc)%atoms(jatom)%ff_type) == check_array(icheck,1) .or. & 
		      trim(molecules(imolc)%atoms(jatom)%ff_type) == check_array(icheck,2))) then 
		  if (print_extra) write(*,*)' already defined this pair'
		  check_index = .false.
		  EXIT
	    endif
          enddo
	  if (check_index) then
	      icount = icount + 1
	      if (print_extra) write(*,*)' icount is ', icount
	      check_array(icount,1) = trim(molecules(imolc)%atoms(iatom)%ff_type) 
	      check_array(icount,2) = trim(molecules(imolc)%atoms(jatom)%ff_type)

	      write(222,*)'&LENNARD-JONES'
	      write(222,*)'  ATOMS ', trim(molecules(imolc)%atoms(iatom)%ff_type), '     ',&
		      trim(molecules(imolc)%atoms(jatom)%ff_type)
	      write(222,*)'  EPSILON ', lj_ff(imatch)%eps_units, dsqrt( (lj_ff(imatch)%eps * lj_ff(jmatch)%eps) )/&
		      real(max(molecules(imolc)%atoms(iatom)%nbead, molecules(imolc)%atoms(jatom)%nbead))
	      write(222,*)'  SIGMA ', lj_ff(imatch)%sigma_units, (lj_ff(imatch)%sigma + lj_ff(jmatch)%sigma)/2.0d0
	      write(222,*)'  RCUT ', lj_ff(imatch)%sigma_units, min(lj_ff(imatch)%rcut_value, lj_ff(jmatch)%rcut_value)
	      write(222,*)'&END LENNARD-JONES'
	  endif
	endif
      enddo




      if (molecules(imolc)%num_repeat > 0) then 
        start_index = imolc
      else
        start_index = imolc+1
      endif
      do jmolc = start_index, SIZE(molecules)
	if (molecules(jmolc)%first .or. jmolc == imolc) then
          do jatom = 1, molecules(jmolc)%natoms
            jmatch = 0
            do i_ff = 1, lj_num
	      label_j = len(trim(molecules(jmolc)%atoms(jatom)%ff_type))
	      if (lj_ff(i_ff)%atom == molecules(jmolc)%atoms(jatom)%ff_type(label_j-1:)) then
	        jmatch = i_ff
		EXIT
	      endif
	    enddo
	    if (imatch == 0 .or. jmatch == 0) then 
              write(*,*) 'combination does not have LJ input'
	    else
              check_index = .true.
              do icheck = 1, icount
		if (trim(molecules(imolc)%atoms(iatom)%ff_type) == trim(molecules(jmolc)%atoms(jatom)%ff_type)) then
		  if (trim(molecules(imolc)%atoms(iatom)%ff_type) == trim(check_array(icheck,1)) .and. &
			  trim(molecules(jmolc)%atoms(jatom)%ff_type) == trim(check_array(icheck,2))) then
		    check_index = .false.
		    EXIT
		  endif
		else if (trim(molecules(imolc)%atoms(iatom)%ff_type) == trim(check_array(icheck,1)) &
			.or.trim(molecules(imolc)%atoms(iatom)%ff_type) == trim(check_array(icheck,2))) then
		  if (trim(molecules(jmolc)%atoms(jatom)%ff_type) == trim(check_array(icheck,1)) &
			  .or. trim(molecules(jmolc)%atoms(jatom)%ff_type) == trim(check_array(icheck,2))) then
		    check_index = .false.
		    EXIT
		  endif
	        endif
              enddo
	      if (check_index) then
	        icount = icount + 1
	        if (print_extra) write(*,*)' icount is ', icount
	        check_array(icount,1) = trim(molecules(imolc)%atoms(iatom)%ff_type) 
	        check_array(icount,2) = trim(molecules(jmolc)%atoms(jatom)%ff_type)

	        write(222,*)'&LENNARD-JONES'
	        write(222,*)'  ATOMS ', trim(molecules(imolc)%atoms(iatom)%ff_type), '     ',&
		      trim(molecules(jmolc)%atoms(jatom)%ff_type)
	        write(222,*)'  EPSILON ', lj_ff(imatch)%eps_units, dsqrt( (lj_ff(imatch)%eps * lj_ff(jmatch)%eps) )/&
		      real(max(molecules(imolc)%atoms(iatom)%nbead, molecules(jmolc)%atoms(jatom)%nbead))
	        write(222,*)'  SIGMA ', lj_ff(imatch)%sigma_units, (lj_ff(imatch)%sigma + lj_ff(jmatch)%sigma)/2.0d0
	        write(222,*)'  RCUT ', lj_ff(imatch)%sigma_units, min(lj_ff(imatch)%rcut_value, lj_ff(jmatch)%rcut_value)
	        write(222,*)'&END LENNARD-JONES'
	      endif
	      check_index = .true.
	    endif
	  enddo
	endif
      enddo



      if (all_classical_ff) then
        do i_ff = 1, lj_num
	    if (imatch == 0) then 
              write(*,*) 'combination does not have LJ input'
	    else
              check_index = .true.
              do icheck = 1, icount
	        !write(*,*)' check_array '
	        if ((trim(molecules(imolc)%atoms(iatom)%ff_type) == check_array(icheck,1) .or. &
		    trim(molecules(imolc)%atoms(iatom)%ff_type) == check_array(icheck,2)) .and. &
		    (trim(lj_ff(i_ff)%atom) == check_array(icheck,1) .or. & 
		      trim(lj_ff(i_ff)%atom) == check_array(icheck,2))) then 
		  write(*,*)' already defined this pair'
		  check_index = .false.
		  EXIT
	        endif
              enddo
	      if (check_index) then
	        icount = icount + 1
	        if (print_extra) write(*,*)' icount is ', icount
	        check_array(icount,1) = trim(molecules(imolc)%atoms(iatom)%ff_type) 
	        check_array(icount,2) = trim(lj_ff(i_ff)%atom)

	        write(222,*)'&LENNARD-JONES'
	        write(222,*)'  ATOMS ', trim(molecules(imolc)%atoms(iatom)%ff_type), '     ',&
		      trim(lj_ff(i_ff)%atom)
	        write(222,*)'  EPSILON ', lj_ff(imatch)%eps_units, dsqrt( (lj_ff(imatch)%eps * lj_ff(i_ff)%eps) )/&
		      real(molecules(imolc)%atoms(iatom)%nbead)
	        write(222,*)'  SIGMA ', lj_ff(imatch)%sigma_units, (lj_ff(imatch)%sigma + lj_ff(i_ff)%sigma)/2.0d0
	        write(222,*)'  RCUT ', lj_ff(imatch)%sigma_units, lj_ff(imatch)%rcut_value
	        write(222,*)'&END LENNARD-JONES'
	      endif
	    endif

        enddo
      endif


    enddo
    endif
  enddo

	       

	
       
	        

endif
      








if (gen_ff .eqv. .true. .and. gen_psf .eqv. .false.) then 
  open(222, file='FF_RP_file.pot')
  open(111,file='FF_info.inp')

  read(111,*)
  read(111,*) ff_type, bond_num
  allocate(bond_ff(bond_num))
  do i = 1, bond_num
    read(111,*) bond_ff(i)%atom1, bond_ff(i)%atom2, bond_ff(i)%k, bond_ff(i)%k_units, bond_ff(i)%r0, bond_ff(i)%r_units
  enddo

  do imolc = 1, SIZE(molecules)
    do iatom = 1, molecules(imolc)%natoms
      do ibond = 1, bond_num
	      if (print_extra) write(*,*)' imolc is ', trim(molecules(imolc)%atoms(iatom)%ff_type)
	      if (print_extra) write(*,*)' bond_ff is ', bond_ff(ibond)%atom1
        if (trim(molecules(imolc)%atoms(iatom)%ff_type(1:2)) == bond_ff(ibond)%atom1) then
		write(*,*)' iatom matches '
          do jatom = 1, molecules(imolc)%natoms
	    if (iatom /= jatom) then
	      if (trim(molecules(imolc)%atoms(jatom)%ff_type(1:2)) == bond_ff(ibond)%atom2) then
		      if (print_extra) write(*,*)' jatom matches '

                write(222,*)'&BOND'
                write(222,*)'  ATOMS ', trim(molecules(imolc)%atoms(iatom)%ff_type), '      ', &
			trim(molecules(imolc)%atoms(jatom)%ff_type)
                write(222,*)'  K     ', trim(bond_ff(ibond)%k_units), bond_ff(ibond)%k/&
			real(max(molecules(imolc)%atoms(iatom)%nbead, &
			molecules(imolc)%atoms(jatom)%nbead))
                write(222,*)'  R0    ', bond_ff(ibond)%r_units, bond_ff(ibond)%r0
		write(222,*)'  KIND G87 '
                write(222,*)'&END BOND'
	      endif
	    endif
          enddo
	endif
      enddo
    enddo
  enddo
         
  close(111)
  close(222)
endif 
      
    

END PROGRAM make_pimd_input
