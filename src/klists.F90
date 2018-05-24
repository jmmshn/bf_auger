program klists

  use base,   only: r8, freeunit, bohr2cm, kelvin2au, ivbm, icbm, &
      nr1, nr2, nr3, nb, cEIGfile
  use kgrids, only: setup_kgrid, find_knum, firstbz, &
       nktot, equiv_ee, wk_e, kgbz, &
       nktot_h, equiv_hh, wk_h, map_he
  use energies, only: setup_energies, fill_electrons, fill_holes, &
       elda, T_Kelvin, kT, density_au, density_cm3, evbm, ecbm
  
  implicit none

  !-- part 1
  integer :: nk, fu, ios

  !-- part 2
  integer :: nkh, nke, nkeep, log, index
  real(r8) :: efermi, ecutoff

  !-- calculation of k4
  logical, allocatable :: isinlist(:)
  integer, allocatable :: klist_e_ful(:)
  integer, allocatable :: klist_h_ful(:)
  integer, allocatable :: klist_h_irr(:)
  integer :: nocc_e_ful
  integer :: nocc_h_ful
  integer :: nocc_h_irr
  integer :: k1, k2, k3
  real(r8) :: kdiff(3)

  integer :: G_umklapp(3)
  real(r8) :: kbz(3)

! ============================================================
! Get the user inputs
! Using NAMELIST and routines taken from QE
!
! ============================================================
  NAMELIST / auger / nb, nr1, nr2, nr3, ivbm, &
      cEIGfile
  NAMELIST / klist / T_Kelvin, density_cm3
! ========== Defaults
nb = 0 ; nr1 = 0 ; nr2 = 0 ; nr3 = 0; 
! ========== Read in 
  CALL input_from_file ( )
! 
  READ (5, auger, iostat = ios)
  READ (5, klist, iostat = ios)
  close(5)
! 
  write(6,*) '======== Parameters ME_INP ========='
  write(6,*) 'nb           =', nb

  write(6,'(A20,3i4)') 'nr1,  nr2,  nr3  = ', nr1, nr2, nr3
  write(6,*) 'T_Kelvin     = ', T_Kelvin
  write(6,*) 'density_cm3  = ', density_cm3
  write(6,*) 'ivbm         = ', ivbm
  write(6,*) 'cEIGfile     = ', cEIGfile
  write(6,*) '===================================='
  

! ========== Conversions
  kT = T_Kelvin * kelvin2au
  density_au = density_cm3 * (bohr2cm)**3
  icbm = ivbm + 1
  
!-- part 1

  call setup_kgrid

  write (*,"(a)") 'writing irreducible k-point list'
  call freeunit(fu)
  open (unit=fu, file='klist.elec.irr')
  nkeep = 0
  do nk = 1, nktot
     if (nk .eq. equiv_ee(nk)) then
        write(fu, "(i7,3f13.9,i3)") nk, kgbz(1:3,nk), wk_e(nk)
        nkeep = nkeep + 1
     end if
  end do
  close(fu)
  write(*,"(a,1x,i0,1x,a)") '... found', nkeep, 'irreducible k-points'

  !-- part 2

  call freeunit(log)
  open (unit=log, file='klist.dat')
  write (log, "(i7,3x,a)") nkeep, 'irreducible k-points'

  call setup_energies(log)

  allocate(isinlist(nktot))

  call fill_electrons(efermi, ecutoff)

  !OPEN(UNIT=fu, FILE="density_efermi_elec.dat")
  !WRITE(fu,*) "# density(a.u.) efermi(Ha)"
  !WRITE(fu,*) density_au, efermi
  !CLOSE(fu)

  !OPEN(UNIT=fu, FILE='nkT_cutoff_elec.dat')
  !WRITE(fu,*) iT
  !CLOSE(fu)

  write (*,"(a)") 'writing occupied electron k-points'
  call freeunit(fu)
  open (unit=fu, file='klist.elec.occ.ful')
  isinlist = .false.
  nkeep = 0
  do nk = 1, nktot
     if (elda(icbm,equiv_ee(nk)).lt.ecutoff) then
        write(fu, "(i7,3f13.9,i3)") nk, kgbz(1:3,nk), wk_e(nk)
        isinlist(nk) = .true.
        nkeep = nkeep + 1
     end if
  end do
  close(fu)
  write(*,"(a,1x,i0,1x,a)") '... found', nkeep, 'k-points'

  nocc_e_ful = nkeep
  allocate(klist_e_ful(nocc_e_ful))
  index = 1
  do nk = 1, nktot
     if(isinlist(nk)) then
        klist_e_ful(index) = nk
        index = index + 1
     end if
  end do

  write (log, "(2i7,3x,a)") nocc_e_ful, nktot, 'e_occ_ful, e_total'
  write (log, "(2es30.15,3x,a)") efermi-ecbm, ecutoff-ecbm, '   elec E_fermi, E_cutoff (Ha)'

  call fill_holes(efermi, ecutoff)

  !OPEN(UNIT=fu, FILE='density_efermi_hole.dat')
  !WRITE(fu,*) "# density(a.u.) efermi(Ha)"
  !WRITE(fu,*) density_au, efermi
  !CLOSE(fu)

  !OPEN(UNIT=fu, FILE='nkT_cutoff_hole.dat')
  !WRITE(fu,*) iT
  !CLOSE(fu)

  ! NEED TO WRITE
  ! INFO: nktot, kT, density
  ! efermi, ecutoff
  ! nk_e_occ

  write (*,"(a)") 'writing occupied hole k-points'
  call freeunit(fu)
  open (unit=fu, file='klist.hole.occ.ful')
  isinlist = .false.
  nkeep = 0
  do nkh = 1, nktot_h ; nke = map_he(nkh)
     if ( elda(ivbm,equiv_ee(nke)) .gt. ecutoff ) then
        write(fu, "(i7,3f13.9,i3)") nke, kgbz(1:3,nke), wk_h(nkh)
        isinlist(nke) = .true.
        nkeep = nkeep + 1
     end if
  end do
  close(fu)
  write(*,"(a,1x,i0,1x,a)") '... found', nkeep, 'k-points'

  nocc_h_ful = nkeep
  allocate(klist_h_ful(nocc_h_ful))
  index = 1
  do nk = 1, nktot
     if(isinlist(nk)) then
        klist_h_ful(index) = nk
        index = index + 1
     end if
  end do

  open (unit=fu, file='klist.hole.occ.irr')
  isinlist = .false.
  nkeep = 0
  do nkh = 1, nktot_h ; nke = map_he(nkh)
     if ( elda(ivbm,equiv_ee(nke)) .gt. ecutoff .and. nkh .eq. equiv_hh(nkh)) then
        write(fu, "(i7,3f13.9,i3)") nke, kgbz(1:3,nke), wk_h(nkh) !-- cannot do weight afterwards
        isinlist(nke) = .true.
        nkeep = nkeep + 1
     end if
  end do
  close(fu)
  write(*,"(a,1x,i0,1x,a)") '... found', nkeep, 'irreducible k-points'

  nocc_h_irr = nkeep
  allocate(klist_h_irr(nocc_h_irr))
  index = 1
  do nk = 1, nktot
     if(isinlist(nk)) then
        klist_h_irr(index) = nk
        index = index + 1
     end if
  end do

  write (log, "(3i7,3x,a)") nocc_h_irr, nocc_h_ful, nktot_h, 'h_occ_irr, h_occ_ful, h_total'
  write (log, "(2es30.15,3x,a)") efermi-evbm, ecutoff-evbm, '   hole E_fermi, E_cutoff (Ha)'

  !-- eeh k4
  write (*,"(a)") 'writing eeh k4'
  isinlist = .false.
  do k1 = 1, nocc_e_ful
  do k2 = 1, nocc_e_ful
  do k3 = 1, nocc_h_irr
     kdiff(:) = kgbz(:, klist_e_ful(k1)) + kgbz(:, klist_e_ful(k2)) - kgbz(:, klist_h_irr(k3))
     call find_knum (nk, kdiff)
     isinlist(nk) = .true.
     call firstbz (kdiff, kbz, G_umklapp)
     if (.not. all(G_umklapp .eq. 0) ) then
        write (*, '(a18,i8,i8,i8,i8,i8,i8)' ) &
            'found an Umklapp:', &
            k1,nocc_e_ful, &
            k2,nocc_e_ful, &
            k3,nocc_h_irr
     end if
  end do
  end do
  end do

  open (unit=fu, file='klist.eeh4')
  nkeep = 0
  do nk = 1, nktot
     if(isinlist(nk)) then
        write(fu, "(i7,3f13.9)") nk, kgbz(1:3,nk)
        nkeep = nkeep + 1
     end if
  end do
  close(fu)
  write(*,"(a,1x,i0,1x,a)") '... found', nkeep, 'eeh k4-points'

  write (log, "(i7,3x,a)") nkeep, 'k4 eeh'

  !-- hhe k4
  write (*,"(a)") 'writing hhe k4'
  isinlist = .false.
  do k1 = 1, nocc_h_ful
  do k2 = 1, nocc_h_irr
  do k3 = 1, nocc_e_ful
     kdiff(:) = kgbz(:, klist_h_ful(k1)) + kgbz(:, klist_h_irr(k2)) - kgbz(:, klist_e_ful(k3))
     call find_knum (nk, kdiff)
     isinlist(nk) = .true.
     call firstbz (kdiff, kbz, G_umklapp)
     if (.not. all(G_umklapp .eq. 0) ) then
        write (*, '(a18,i8,i8,i8,i8,i8,i8)' ) &
            'found an Umklapp:', &
            k1,nocc_h_ful, &
            k2,nocc_h_irr, &
            k3,nocc_e_ful
     end if
  end do
  end do
  end do

  open (unit=fu, file='klist.hhe4')
  nkeep = 0
  do nk = 1, nktot
     if(isinlist(nk)) then
        write(fu, "(i7,3f13.9)") nk, kgbz(1:3,nk)
        nkeep = nkeep + 1
     end if
  end do
  close(fu)
  write(*,"(a,1x,i0,1x,a)") '... found', nkeep, 'eeh k4-points'

  write (log, "(i7,3x,a)") nkeep, 'k4 hhe'

  close (log)
  ! check what the shift does
  !call freeunit(fu)
  !open (unit=fu, file='elda.dat')
  !do nk = 1 ,20
      !write(fu,"(i,10f7.3)")  nk, elda(1:10,nk)*Ha2eV
  !end do

  close(fu)

end program klists
