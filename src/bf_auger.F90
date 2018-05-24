program bf_auger
  !-- E. Kioupakis, 01/18/2010
  !-- Calculate the Auger coefficient in a brute-force way, eeh
  !-- J. Shen, 01/05/2016  
  !-- modified the program to take spinor wavefunctions
  !-- J. Shen, 01/05/2017  
  !-- modified the program to take namelist inputs

  use base, only: dp, pi, Ha2eV, bohr2cm, time_au2sec, &
       freeunit, read_elda, read_wfn, nr1, nr2, nr3, ivbm, icbm,&
       nb, cEIGfile, lsoc
  use kgrids, only: nktot, nktot_h, kgbz, &
       Vcell => omega, &
       setup_kgrid, find_knum, firstbz
  use me_fft_direct, only: setup_fft, setup_epsilon, coulomb_me_fft
  implicit none

  !-- width of gaussian functions approximating the delta function
  integer, parameter :: neta = 8
  integer :: ieta
  real(dp) :: eta(neta) = (/ 0.01_dp, 0.02_dp, 0.03_dp, 0.05_dp, 0.1_dp, 0.2_dp, 0.3_dp, 0.5_dp /) ! in eV

  integer, parameter:: ngap = 401 ! sample size of gaps
  integer :: U_OUT
  integer :: igap
  real(dp) :: Egap_min  ! energy in hartree
  real(dp) :: Egap_step ! energy step in hartree
  real(dp) :: Egap(ngap)
  real(dp) :: Delta_scissors(ngap)

  integer :: nx, ny, nz, np

  complex(dp), allocatable :: wfn1(:,:), wfn2(:,:), wfn3(:,:), wfn4(:,:)

  !integer :: ivbm, icbm, nb
  real(dp) :: n_free
  real(dp) :: kT

  !-- occupied k-points and symmetry weights
  integer :: nk_e_ful, nk_irr
  integer :: nk_h_ful
  integer :: nk_h_irr
  integer, allocatable :: klist_irr(:)
  integer, allocatable :: klist_e_ful(:)
  integer, allocatable :: klist_h_ful(:)
  integer, allocatable :: klist_h_irr(:)
  real(dp) :: weight_e_ful
  real(dp) :: weight_h_ful
  real(dp), allocatable :: weight_h_irr(:)
  integer, allocatable :: ik1(:), ik2(:), ik3(:)
  integer :: i_loop, k4num, k1, k2, k3, k4

  integer :: n1, n2, n3, n4
  real(dp) :: f1, f2, f3, f4

  integer :: ik, i1, i2, i3

  real(dp) :: efermi_e, efermi_h, ecutoff_e, ecutoff_h
  real(dp) :: evbm, ecbm, egap_lda
  integer  :: evbm_idx, ecbm_idx
  real(dp), allocatable :: elda1(:), elda2(:), elda3(:), elda4(:)

  real(dp) :: auger_coef(ngap,neta), auger_coef_global(ngap,neta)

  real(dp) :: kdiff(3), kbz(3), mtrans(3)
  integer :: G_Umklapp(3)

  complex(dp) :: auger_d, auger_x
  real(dp) :: auger_me
  real(dp) :: DeltaE

  CHARACTER(LEN=10) :: c
  integer :: fu
  real(dp) :: rdummy
  integer :: idummy, ios

  integer nprocs, ident
  real(dp) :: t_start, t1, t2 ! timing
  !integer  :: t_counter 
#ifdef MPI
  include 'mpif.h'
  integer :: status(MPI_STATUS_SIZE)
  integer :: ierr
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, ident, ierr)
#else
  nprocs = 1
  ident = 0
#endif

! ============================================================
! Get the user inputs
! Using NAMELIST and routines taken from QE
!
! ============================================================
! ========== Read in 
  NAMELIST / auger / nb, nr1, nr2, nr3, cEIGfile, lsoc, &
      Egap_min, Egap_step, ivbm
!call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  if (ident .eq. 0) then
! ========== Defaults
    nb = 0 ; nr1 = 0 ; nr2 = 0 ; nr3 = 0; 
    cEIGfile='eigen.val'; lsoc=.False.
    Egap_min=0.1; Egap_step=0.02
! ==========
    CALL input_from_file()
! 
    READ (5, auger, iostat = ios)
    close(5)
    if (ios.ne.0) then
      write(6,*) 'Auger input ERRO IOS = ' , ios
      STOP
    end if
! 
    call freeunit (U_OUT)
    OPEN (U_OUT,FILE="auger_eeh.log", ACTION='WRITE',STATUS='REPLACE')
    write(U_OUT,*) '======== Parameters ME_INP ========='
    write(U_OUT,*) 'nb=', nb

    write(U_OUT,'(A20,3i4)') 'nr1,  nr2,  nr3  = ', nr1, nr2, nr3
    write(U_OUT,*) 'cEIGfile  = ', cEIGfile
    write(U_OUT,*) 'lsoc  = ', lsoc
    write(U_OUT,*) 'ivbm  = ', ivbm
    write(U_OUT,*) 'Egap_min  = ', Egap_min, Egap_min/ha2eV
    write(U_OUT,*) 'Egap_step  = ', Egap_step, Egap_step/ha2eV
    write(U_OUT,*) '===================================='
    
    Egap_min = Egap_min/ha2eV
    Egap_step = Egap_step/ha2eV
  end if
! ========== Broadcast
  !NAMELIST / auger / nb, nr1, nr2, nr3, cEIGfile, lsoc, &
      !Egap_min, Egap_step, ivbm
#ifdef MPI
  CALL mpi_bcast( nb ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr )
  CALL mpi_bcast( nr1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr )
  CALL mpi_bcast( nr2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr )
  CALL mpi_bcast( nr3,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr )
  CALL mpi_bcast( ivbm,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr )
  CALL mpi_bcast( lsoc ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr )

  CALL mpi_bcast( Egap_min ,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr )
  CALL mpi_bcast( Egap_step ,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr )

  CALL mpi_bcast( cEIGfile,50,MPI_CHAR,0,MPI_COMM_WORLD,ierr )
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
  !print *, 'ident ', ident,', val ', cEIGfile
  
  !CALL cpu_time(t_start)
  
  eta = eta / Ha2eV
  call freeunit (fu)

  OPEN (UNIT=fu, FILE="kpoints/k_1/UNK00001.1", FORM='UNFORMATTED')
  READ (fu) nx, ny, nz, idummy, nb
  CLOSE (fu)
  np = nx * ny * nz

  ALLOCATE(wfn1(np,nb))
  ALLOCATE(wfn2(np,nb))
  ALLOCATE(wfn3(np,nb))
  ALLOCATE(wfn4(np,nb))

  allocate (elda1(nb))
  allocate (elda2(nb))
  allocate (elda3(nb))
  allocate (elda4(nb))

  open (unit=fu, file='klist.dat', action='read')
  read (fu,*) nk_irr
  read (fu,*) ivbm
  read (fu,*) evbm_idx, evbm
  read (fu,*) ecbm_idx, ecbm
  read (fu,*) n_free
  read (fu,*) kT
  read (fu,*) nk_e_ful
  read (fu,*) efermi_e, ecutoff_e
  read (fu,*) nk_h_irr, nk_h_ful
  read (fu,*) efermi_h, ecutoff_h
  close (fu)
  if (ident.eq.0) write(*,*)  'Number of irreducible kpts = ', nk_irr

  icbm = ivbm + 1

  ! CONT_ create a continuation file (unformated) RESTART
  ! CONT_ copy existing file over to a tmp file


  call setup_kgrid
  call setup_fft (nx, ny, nz)
  call setup_epsilon (n_free, kT)


  !! find the vbm and cbm energies
  !allocate (klist_irr(nk_irr))
  !open (unit=fu, file='klist.elec.irr', action='read')
  !do ik = 1, nk_irr
  !   read (fu,*) klist_irr(ik)
  !end do
  !close (fu)

  !call read_elda (1, nb, elda1, fu)
  egap_lda = ecbm - evbm



  !print *, evbm*ha2eV, ecbm*ha2eV, egap_lda*ha2eV
  DO igap = 1, ngap
     Egap(igap) = Egap_min + (igap-1)*Egap_step
     Delta_scissors(igap) = Egap(igap) - egap_lda !-- The correction needed to the LDA eigenvalues to achieve the gap Egap(i)
     !write(*,*) egap(igap)*ha2eV, delta_scissors(igap)*ha2eV
  END DO

  allocate (klist_e_ful(nk_e_ful))
  open (unit=fu, file='klist.elec.occ.ful', action='read')
  do ik = 1, nk_e_ful
     read (fu,*) klist_e_ful(ik)
  end do
  close (fu)
  weight_e_ful = 1.0_dp / real(nktot, dp)

  allocate (klist_h_ful(nk_h_ful))
  open (unit=fu, file='klist.hole.occ.ful', action='read')
  do ik = 1, nk_h_ful
     read (fu,*) klist_h_ful(ik)
  end do
  close (fu)
  weight_h_ful = 1.0_dp / real(nktot_h, dp)

  allocate (klist_h_irr(nk_h_irr))
  allocate (weight_h_irr(nk_h_irr))
  open (unit=fu, file='klist.hole.occ.irr', action='read')
  do ik = 1, nk_h_irr
     read (fu,*) klist_h_irr(ik), rdummy, rdummy, rdummy, idummy
     weight_h_irr(ik) = real(idummy, dp) / real(nktot_h, dp)
  end do
  close (fu)
  
#ifdef EEH
  k4num = nk_e_ful*nk_e_ful*nk_h_irr
#else
  k4num = nk_h_ful*nk_h_irr*nk_e_ful
#endif

  allocate (ik1(k4num))
  allocate (ik2(k4num))
  allocate (ik3(k4num))
  i_loop = 0
#ifdef EEH
  do i1 = 1, nk_e_ful
  do i2 = 1, nk_e_ful
  do i3 = 1, nk_h_irr
#else
  do i1 = 1, nk_h_ful
  do i2 = 1, nk_h_irr
  do i3 = 1, nk_e_ful
#endif

     i_loop = i_loop + 1
     ik1(i_loop) = i1
     ik2(i_loop) = i2
     ik3(i_loop) = i3
  end do
  end do
  end do

  auger_coef = 0.0_dp
  
  iloop: do i_loop = ident + 1, k4num, nprocs
     if (ident.eq.0) write (*,*) i_loop, ' of ', k4num

#ifdef EEH
     k1 = klist_e_ful(ik1(i_loop))
     k2 = klist_e_ful(ik2(i_loop))
     k3 = klist_h_irr(ik3(i_loop))
#else
     k1 = klist_h_ful(ik1(i_loop))
     k2 = klist_h_irr(ik2(i_loop))
     k3 = klist_e_ful(ik3(i_loop))
#endif
     kdiff(:) = kgbz(:, k1) + kgbz(:, k2) - kgbz(:, k3)
     call find_knum (k4, kdiff)
     call firstbz (kdiff, kbz, G_umklapp)

     call read_elda (k1, nb, elda1, fu)
     call read_elda (k2, nb, elda2, fu)
     call read_elda (k3, nb, elda3, fu)
     call read_elda (k4, nb, elda4, fu)

     call read_wfn (k1, np, nb, wfn1, fu)
     call read_wfn (k2, np, nb, wfn2, fu)
     call read_wfn (k3, np, nb, wfn3, fu)
     call read_wfn (k4, np, nb, wfn4, fu)
     

#ifdef EEH
     n1_loop: DO n1 = icbm, nb
        IF ( (elda1(n1)-ecbm) .GT.  ecutoff_e ) cycle n1_loop
        
     n2_loop: DO n2 = icbm, nb
        IF ( (elda2(n2)-ecbm) .GT.  ecutoff_e ) cycle n2_loop
           
     n3_loop: DO n3 = 1, ivbm
        IF (-(elda3(n3)-evbm) .GT. -ecutoff_h ) cycle n3_loop
              
     n4_loop: DO n4 = icbm, nb

       !print *, 'e1 =' ,elda1(n1)*ha2eV,', e2 = ',elda2(n2)*ha2eV,', e3, ',elda3(n3)*ha2eV
        
        !write(*,*) 'n1', n1, 'n2', n2, 'n3', n3
        f1 = 1.0_dp / (1.0_dp + EXP( ((elda1(n1)-ecbm) - efermi_e) / kT) )
        f2 = 1.0_dp / (1.0_dp + EXP( ((elda2(n2)-ecbm) - efermi_e) / kT) )
        f3 = 1.0_dp / (1.0_dp + EXP( (efermi_h - (elda3(n3)-evbm)) / kT) )
        f4 = 1.0_dp / (1.0_dp + EXP( ((elda4(n4)-ecbm) - efermi_e) / kT) )
#else
     n1_loop: DO n1 = 1, ivbm
        IF (-(elda1(n1)-evbm) .GT. -ecutoff_h ) cycle n1_loop
        
     n2_loop: DO n2 = 1, ivbm
        IF (-(elda2(n2)-evbm) .GT. -ecutoff_h ) cycle n2_loop

     n3_loop: DO n3 = icbm, nb
        IF ( (elda3(n3)-ecbm) .GT.  ecutoff_e ) cycle n3_loop

     n4_loop: DO n4 = 1, ivbm
        ! if (ident.eq.0)  write(*,*) 'e1 =' ,elda1(n1)*ha2eV,', e2 = ',elda2(n2)*ha2eV,', e3, ',elda3(n3)*ha2eV
        
        !write(*,*) 'n1', n1, 'n2', n2, 'n3', n3
        f1 = 1.0_dp / (1.0_dp + EXP( (efermi_h - (elda1(n1)-evbm)) / kT) )
        f2 = 1.0_dp / (1.0_dp + EXP( (efermi_h - (elda2(n2)-evbm)) / kT) )
        f3 = 1.0_dp / (1.0_dp + EXP( ((elda3(n3)-ecbm) - efermi_e) / kT) )
        f4 = 1.0_dp / (1.0_dp + EXP( (efermi_h - (elda4(n4)-evbm)) / kT) )
#endif

        ! CONT_ READ auger_me here if the ME is already here, then skip to the
        ! next CONT_
        mtrans = kgbz(:,k1) - kgbz(:,k3)
        auger_d = coulomb_me_fft (wfn1(:,n1), wfn2(:,n2), wfn3(:,n3), wfn4(:,n4), G_umklapp, mtrans)

        mtrans = kgbz(:,k2) - kgbz(:,k3)
        auger_x = coulomb_me_fft (wfn2(:,n2), wfn1(:,n1), wfn3(:,n3), wfn4(:,n4), G_umklapp, mtrans)

        auger_me = ABS(auger_d-auger_x)**2 + ABS(auger_d)**2 + ABS(auger_x)**2
     
        ! CONT_ set auger_me to the value in the file
      !t_counter
! Since we cannot finish the calculation at once 
! --- write the auger_coef calculated on each node
! --- open each file and delete that last line
#ifdef EEH
#else
#endif

     !if (ident.eq.0) CALL cpu_time(t1)
        
        !write(*,*) 'Auger_me', auger_me
        !write(*,*) 'With fermi_fact', weight_e_ful * weight_e_ful * weight_h_irr(ik3(i_loop)) * f1 * f2 * f3 * (1 - f4) * auger_me
        DO ieta = 1, neta
           DO igap = 1, ngap
#ifdef EEH
              DeltaE = elda1(n1) + elda2(n2) + Delta_scissors(igap) - elda3(n3) - elda4(n4)
              auger_coef(igap,ieta) = auger_coef(igap,ieta) + &
                   weight_e_ful * weight_e_ful * weight_h_irr(ik3(i_loop)) * &
                   f1 * f2 * f3 * (1 - f4) * auger_me * &
                   EXP( -DeltaE**2 / eta(ieta)**2 ) / SQRT(pi*eta(ieta)**2)
#else
              DeltaE = elda1(n1) + elda2(n2) - Delta_scissors(igap) - elda3(n3) - elda4(n4)
              auger_coef(igap,ieta) = auger_coef(igap,ieta) + &
                   weight_h_ful * weight_h_irr(ik2(i_loop)) * weight_e_ful * &
                   f1 * f2 * f3 * (1 - f4) * auger_me * &
                   EXP( -DeltaE**2 / eta(ieta)**2 ) / SQRT(pi*eta(ieta)**2)
               ! write the auger_coef calculated on each node
#endif
           END DO
        END DO
     !if (ident.eq.0) CALL cpu_time(t2)
     !if (ident.eq.0) write (*,*) i_loop, ' of ', k4num, 'Dt:', t2-t1

     end DO n4_loop
     end DO n3_loop
     end DO n2_loop
     end DO n1_loop
  end do iloop

  !-- Prefactor
  auger_coef = auger_coef * 2.0_dp * (2.0_dp/pi)*(1.0_dp/Vcell**3) / n_free**3
  !-- Convert to cm6sec-1
  auger_coef = auger_coef *(bohr2cm)**6/time_au2sec

  auger_coef_global = 0.0_dp
#ifdef MPI
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  call MPI_REDUCE( auger_coef, auger_coef_global, ngap*neta, MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#else
  auger_coef_global = auger_coef
#endif
  
  IF ( ident .EQ. 0 ) THEN
     DO ieta = 1, neta
        WRITE(c,"(i4.4)") ieta
#ifdef EEH
        OPEN (UNIT=fu, FILE="auger_coef_eeh_vs_gap_" // TRIM(c) // ".dat")
#else
        OPEN (UNIT=fu, FILE="auger_coef_hhe_vs_gap_" // TRIM(c) // ".dat")
#endif
        WRITE(fu,"(a,f7.3)") "# Egap(eV), Auger_coef(cm6s-1), eta (eV) = ", eta(ieta)*ha2eV
        DO igap = 1, ngap
           if (auger_coef_global(igap,ieta).lt.1.d-200) auger_coef_global(igap,ieta) = 0.d0
           WRITE(fu,"(f7.3,e18.9)") Egap(igap)*ha2eV, auger_coef_global(igap,ieta)
        END DO
        CLOSE(fu)
     END DO
  END IF

#ifdef MPI
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  call MPI_FINALIZE(ierr)
#endif
  
end program bf_auger
