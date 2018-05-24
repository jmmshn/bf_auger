module energies
  use base, only: r8, bohr2cm, kelvin2au, ivbm, icbm, nb
  implicit none
  save

  !include 'includes.f90'
  real(r8) :: density_cm3 
  real(r8) :: T_Kelvin 

  !real(r8) :: density_au, kT
  real(r8) :: kT
  real(r8) :: density_au 
  !real(r8), parameter :: n_free = density_au
  
  real(r8), allocatable :: elda(:,:)
  real(r8) :: evbm, ecbm, egap_lda
  integer  :: evbm_idx, ecbm_idx
  
  private
  public :: setup_energies, fill_electrons, fill_holes, &
       elda, T_Kelvin, density_cm3, kT, density_au, evbm, ecbm

contains

  subroutine setup_energies(log)
    use base,   only: freeunit, read_elda, ha2eV
    use kgrids, only: nktot, equiv_ee
    implicit none
    integer :: log

    integer :: fu, nk, idummy

    logical :: file_exists
    inquire (file="kpoints/k_1/UNK00001.1", exist=file_exists)
    if (.not.file_exists) then
       stop
    end if

    call freeunit(fu)

    open (unit=fu, file="kpoints/k_1/UNK00001.1", form='UNFORMATTED')
    read (fu) idummy, idummy, idummy, idummy, nb
    close (fu)

    write (*,"(a)") 'reading LDA eigenvalues'
    allocate(elda(nb,nktot))
    elda = -500.0_r8

    call read_elda (1, nb, elda(1:nb,1), fu)
    evbm = elda(ivbm,1)
    ecbm = elda(icbm,1)
    do nk = 1, nktot
       if (nk .eq. equiv_ee(nk)) then ! only read first occurance
          call read_elda (nk, nb, elda(1:nb, nk), fu)
       end if
      if (elda(ivbm,nk)>evbm) then
        evbm_idx = nk
        evbm = elda(ivbm,nk) 
      end if
      if (elda(icbm,nk)<ecbm) then
        ecbm_idx = nk
        ecbm = elda(icbm,nk)
      end if

    end do

    !-- density and temperature
    !OPEN(UNIT=fu, FILE='DENSITY', action='READ')
    !READ(fu,*) density_cm3
    !PRINT*, 'Density = ', density_cm3
    !CLOSE(fu)
    !density_au = density_cm3 * (bohr2cm)**3
    !OPEN(UNIT=fu, FILE='TEMPERATURE', action='READ')
    !READ(fu,*) kT_kelvin
    !PRINT*, 'T = ',kT_kelvin
    !CLOSE(fu)
    !kT = kT_kelvin * kelvin2au

    write (log, "(2i3,3x,a)") ivbm, nb, 'ivbm, nb'
    write (log, "(i8,f19.13,3x,a)") evbm_idx, evbm, 'evbm_idx, evbm (Ha)'
    write (log, "(i8,f19.13,3x,a)") ecbm_idx, ecbm, 'ecbm_idx, ecbm (Ha)'
    write (log, *) density_au, density_cm3, '   density (au, cm3)'
    write (log, *) kT, T_Kelvin, '   temperature (Ha, K)'
  end subroutine setup_energies

  subroutine fill_electrons(efermi, ecutoff)
    use base,   only: Ha2eV
    use kgrids, only: Vcell => omega, &
         nktot, equiv_ee, wk_e
    implicit none

    integer :: nk, ib, step, iT

    real(r8) :: efermi1,  efermi2,  efermi
    real(r8) :: density1, density2, density
    real(r8) :: ecutoff

    write (*,"(a)") 'calculate Efermi for electrons'

    !-- Reference energies to CBM:
    write(*,*)  '========================'
    ecbm=100.
    do nk = 1, nktot
      if (nk .ne. equiv_ee(nk)) cycle
      if (elda(icbm,nk) < ecbm) ecbm = elda(icbm,nk)
    end do
    print *, 'ecbm',  ecbm
    !elda = elda - ecbm


    efermi1 = -20.0_r8/ha2eV
    efermi2 =  20.0_r8/ha2eV

    density1 = 0.0_r8
    density2 = 0.0_r8
  
    do nk = 1, nktot
       if (nk .eq. equiv_ee(nk)) then
          do ib = icbm, nb
             density1 = density1 + real(wk_e(nk),r8)/(1.0_r8 + EXP( (-efermi1+elda(ib,nk))/kT) )
             density2 = density2 + real(wk_e(nk),r8)/(1.0_r8 + EXP( (-efermi2+elda(ib,nk))/kT) )
          end do
       end if
    end do
    !-- prefactors
    density1 = density1 * 2.0_r8/(Vcell*nktot)
    density2 = density2 * 2.0_r8/(Vcell*nktot)

    if ( (density1-density_au)*(density2-density_au) .gt. 0.0_r8 ) then
       print *, 'wrong energy window'
       stop
    end if

    bisect: do step = 1, 100

       efermi = (efermi1 + efermi2)/2.0_r8
       density = 0.0_r8
       do nk = 1, nktot
          if (nk .eq. equiv_ee(nk)) then
             do ib = icbm, nb
                density = density + real(wk_e(nk),r8)/(1.0_r8 + EXP( (-efermi+elda(ib,nk))/kT) )
             end do
          end if
       end do
       density = density * 2.0_r8/(Vcell*nktot)
       
       write (*, "(i4,f14.9,e18.9)") step, efermi*Ha2eV, density/(bohr2cm)**3
       
       if ( abs(density - density_au)/density_au .lt. 1.0E-6 ) then
          exit bisect
       end if
       
       if ( density .gt. density_au ) then
          !-- Root in [1,m]
          efermi2 = efermi
          density2 = density
       else
          !-- Root in [m,2]
          efermi1 = efermi
          density1 = density
       end if
       
    end do bisect

    if (step .ge. 99) then
       stop 'bisecting failed'
    end if
    
    write(*,"(a)") 'calculate cutoff that contains 99% of charge'
    cutoff: do iT = 0, 20
       ecutoff = efermi + iT*kT
       
       density = 0.0_r8
       do nk = 1, nktot
          if (nk .eq. equiv_ee(nk)) then
             do ib = icbm, nb
                if (elda(ib,nk).lt.ecutoff) then
                   density = density + real(wk_e(nk),r8)/(1.0_r8 + EXP( (-efermi+elda(ib,nk))/kT) )
                end if
             end do
          end if
       end do
       density = density * 2.0_r8/(Vcell*nktot)
       
       write(*,"(i3,a,e18.9)") iT, ' kT, density = ', density/(bohr2cm)**3
       if ( abs(density-density_au)/density_au .lt. 1.0E-2_r8 ) then
          exit cutoff
       end if
    end do cutoff
    
  end subroutine fill_electrons

  subroutine fill_holes(efermi, ecutoff)
    use base,   only: Ha2eV
    use kgrids, only: Vcell => omega, &
         equiv_ee, &
         nktot_h, equiv_hh, wk_h, map_he  
    implicit none

    integer :: nkh, nke, ib, step, iT

    real(r8) :: efermi1,  efermi2,  efermi
    real(r8) :: density1, density2, density
    real(r8) :: ecutoff

    write (*,"(a)") 'calculate Efermi for holes'
  
    !-- Reference energies to VBM:
    write(*,*)  '========================'
    evbm=-100.
    do nkh = 1, nktot_h
      if (nkh .ne. equiv_ee(nkh)) cycle
      if (elda(ivbm,nkh) > evbm) evbm = elda(ivbm,nkh)
    end do
    print *, 'evbm',  evbm
    ! print the band edges
    print *, 'BANDs (MIN, MAX)'
    do ib = 1, ivbm
      print *, ib, (MINVAL(elda(ib,:)),MAXVAL(elda(ib,:)))
    end do

    !elda = elda - evbm

    ! print the band edges
    print *, 'BANDs-EVBM (MIN, MAX)'
    do ib = 1, ivbm
      print *, ib, (MINVAL(elda(ib,:)),MAXVAL(elda(ib,:)))
    end do
    
    
    
    efermi1 = -20.0_r8/ha2eV
    efermi2 =  20.0_r8/ha2eV

    density1 = 0.0_r8
    density2 = 0.0_r8
  
    do nkh = 1, nktot_h ; nke = map_he(nkh)
       if (nkh .eq. equiv_hh(nkh)) then
          do ib = 1, ivbm
             density1 = density1 + real(wk_h(nkh),r8)/(1.0_r8 + EXP( (efermi1-elda(ib,equiv_ee(nke)))/kT) )
             density2 = density2 + real(wk_h(nkh),r8)/(1.0_r8 + EXP( (efermi2-elda(ib,equiv_ee(nke)))/kT) )
          end do
       end if
    end do
    !-- prefactors
    density1 = density1 * 2.0_r8/(Vcell*nktot_h)
    density2 = density2 * 2.0_r8/(Vcell*nktot_h)

    if ( (density1-density_au)*(density2-density_au) .gt. 0.0_r8 ) then
       print *, 'wrong energy window'
       stop
    end if

    bisect: do step = 1, 100

       efermi = (efermi1 + efermi2)/2.0_r8
       density = 0.0_r8
       do nkh = 1, nktot_h ; nke = map_he(nkh)
          if (nkh .eq. equiv_hh(nkh)) then
             do ib = 1, ivbm
                density = density + real(wk_h(nkh),r8)/(1.0_r8 + EXP( (efermi-elda(ib,equiv_ee(nke)))/kT) )
             end do
          end if
       end do
       density = density * 2.0_r8/(Vcell*nktot_h)

       write (*, "(i4,f14.9,e18.9)") step, efermi*Ha2eV, density/(bohr2cm)**3

       if ( abs(density - density_au)/density_au .lt. 1.0E-6 ) then
          exit bisect
       end if

       if ( density .lt. density_au ) then
          !-- Root in [1,m]
          efermi2 = efermi
          density2 = density
       else
          !-- Root in [m,2]
          efermi1 = efermi
          density1 = density
       end if

    end do bisect

    if (step .ge. 99) then
       stop 'bisecting failed'
    end if

    write(*,"(a)") 'calculate cutoff that contains 99% of charge'
    cutoff: do iT = 0, 20
       ecutoff = efermi - iT*kT
       
       density = 0.0_r8
       do nkh = 1, nktot_h ; nke = map_he(nkh)
          if (nkh .eq. equiv_hh(nkh)) then
             do ib = 1, ivbm
                if ( elda(ib,equiv_ee(nke)) .gt. ecutoff ) then
                   density = density + real(wk_h(nkh),r8)/(1.0_r8 + EXP( (efermi-elda(ib,equiv_ee(nke)))/kT) )
                end if
             end do
          end if
       end do
       !-- prefactors
       density = density * 2.0_r8/(Vcell*nktot_h)
       
       write(*,"(i3,a,e18.9)") iT, ' kT, density = ', density/(bohr2cm)**3
       if ( abs(density-density_au)/density_au .lt. 1.0E-2_r8 ) then
          exit cutoff
       end if
    end do cutoff

  end subroutine fill_holes

end module energies
