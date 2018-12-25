program bf_auger_hhe
  !-- Calculate the Auger coefficient in a brute-force way, hhe
  !-- E. Kioupakis, 01/18/2010

  use base, only: dp, pi, Ha2eV, bohr2cm, time_au2sec, &
       freeunit, read_elda, read_wfn
  use kgrids, only: nktot, nktot_h, kgbz, &
       Vcell => omega, &
       setup_kgrid, find_knum, firstbz
  use me_fft_direct, only: setup_fft, setup_epsilon, coulomb_me_fft
  implicit none

  !-- width of gaussian functions approximating the delta function
  integer, parameter :: neta = 8
  integer :: ieta
  real(dp) :: eta(neta) = (/ 0.01_dp, 0.02_dp, 0.03_dp, 0.05_dp, 0.1_dp, 0.2_dp, 0.3_dp, 0.5_dp /) ! in eV

#ifdef GAAS
  integer, parameter :: ngap = 201 !GaAs
#endif
#ifdef INAS
  integer, parameter :: ngap = 201 !InAs
#endif
#if defined (INP) || defined (ALN4) || defined (GAN4) || defined (ALN8) || defined (GAN8) || defined (ALGAN8)
  integer, parameter :: ngap = 601
#endif
  integer :: igap
  real(dp), parameter :: Egap_min = 0.1_dp/ha2eV
  real(dp), parameter :: Egap_step = 0.01_dp/ha2eV
  real(dp) :: Egap(ngap)
  real(dp) :: Delta_scissors(ngap)

  integer :: nx, ny, nz, np

  complex(dp), allocatable :: wfn1(:,:), wfn2(:,:), wfn3(:,:), wfn4(:,:)

  integer :: ivbm, icbm, nb
  real(dp) :: n_free
  real(dp) :: kT

  !-- occupied k-points and symmetry weights
  integer :: nk_e_ful
  integer :: nk_h_ful
  integer :: nk_h_irr
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
  integer :: idummy

  integer nprocs, ident
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
  read (fu,*) ivbm
  read (fu,*) n_free
  read (fu,*) kT
  read (fu,*) nk_e_ful
  read (fu,*) efermi_e, ecutoff_e
  read (fu,*) nk_h_irr, nk_h_ful
  read (fu,*) efermi_h, ecutoff_h
  close (fu)

  icbm = ivbm + 1

  call setup_kgrid
  call setup_fft (nx, ny, nz)
  call setup_epsilon (n_free, kT)

  call read_elda (1, nb, elda1, fu)
  evbm = elda1(ivbm)
  ecbm = elda1(icbm)
  egap_lda = ecbm - evbm
  !print *, evbm*ha2eV, ecbm*ha2eV, egap_lda*ha2eV
  DO igap = 1, ngap
     Egap(igap) = Egap_min + (igap-1)*Egap_step
     Delta_scissors(igap) = Egap(igap) - egap_lda !-- The correction needed to the LDA eigenvalues to achieve the gap Egap(i)
     !print *, egap(i)*ha2eV, delta_scissors(i)*ha2eV
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

  k4num = nk_h_ful*nk_h_irr*nk_e_ful

  !-- k4_vec = k_hole_ful(:,k1) + k_hole_irr(:,k2) - k_elec_ful(:,k3)
  allocate (ik1(k4num))
  allocate (ik2(k4num))
  allocate (ik3(k4num))
  i_loop = 0
  do i1 = 1, nk_h_ful
  do i2 = 1, nk_h_irr
  do i3 = 1, nk_e_ful
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

     k1 = klist_h_ful(ik1(i_loop))
     k2 = klist_h_irr(ik2(i_loop))
     k3 = klist_e_ful(ik3(i_loop))
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

     n1_loop: DO n1 = 1, ivbm
        IF (-(elda1(n1)-evbm) .GT. -ecutoff_h ) cycle n1_loop
        
     n2_loop: DO n2 = 1, ivbm
        IF (-(elda2(n2)-evbm) .GT. -ecutoff_h ) cycle n2_loop
           
     n3_loop: DO n3 = icbm, nb
        IF ( (elda3(n3)-ecbm) .GT.  ecutoff_e ) cycle n3_loop
              
     n4_loop: DO n4 = 1, ivbm
        
        f1 = 1.0_dp / (1.0_dp + EXP( (efermi_h - (elda1(n1)-evbm)) / kT) )
        f2 = 1.0_dp / (1.0_dp + EXP( (efermi_h - (elda2(n2)-evbm)) / kT) )
        f3 = 1.0_dp / (1.0_dp + EXP( ((elda3(n3)-ecbm) - efermi_e) / kT) )
        f4 = 1.0_dp / (1.0_dp + EXP( (efermi_h - (elda4(n4)-evbm)) / kT) )

        mtrans = kgbz(:,k1) - kgbz(:,k3)
        auger_d = coulomb_me_fft (wfn1(:,n1), wfn2(:,n2), wfn3(:,n3), wfn4(:,n4), G_umklapp, mtrans)
        
        mtrans = kgbz(:,k2) - kgbz(:,k3)
        auger_x = coulomb_me_fft (wfn2(:,n2), wfn1(:,n1), wfn3(:,n3), wfn4(:,n4), G_umklapp, mtrans)
        
        auger_me = ABS(auger_d-auger_x)**2 + ABS(auger_d)**2 + ABS(auger_x)**2

        DO ieta = 1, neta
           DO igap = 1, ngap
              DeltaE = elda1(n1) + elda2(n2) - Delta_scissors(igap) - elda3(n3) - elda4(n4)
              auger_coef(igap,ieta) = auger_coef(igap,ieta) + &
                   weight_h_ful * weight_h_irr(ik2(i_loop)) * weight_e_ful * &
                   f1 * f2 * f3 * (1 - f4) * auger_me * &
                   EXP( -DeltaE**2 / eta(ieta)**2 ) / SQRT(pi*eta(ieta)**2)
           END DO
        END DO

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
        OPEN (UNIT=fu, FILE="auger_coef_hhe_vs_gap_" // TRIM(c) // ".dat")
        WRITE(fu,"(a,f7.3)") "# Egap(eV), Auger_coef(cm6s-1), eta (eV) = ", eta(ieta)*ha2eV
        DO igap = 1, ngap
           if (auger_coef_global(igap,ieta).lt.1.d-90) auger_coef_global(igap,ieta) = 0.d0
           WRITE(fu,"(f7.3,g18.9)") Egap(igap)*ha2eV, auger_coef_global(igap,ieta)
        END DO
        CLOSE(fu)
     END DO
  END IF

#ifdef MPI
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  call MPI_FINALIZE(ierr)
#endif

end program bf_auger_hhe
