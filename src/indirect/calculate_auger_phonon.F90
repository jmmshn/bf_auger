program phononauger

  use base
  use tableio
  use me_fft
  implicit none

#ifdef MPI
  include 'mpif.h'
  integer :: ierr
#endif
  integer, parameter :: dp = selected_real_kind(15,300)

  real(dp), parameter :: pi = 3.14159265358979323846264338327950288_dp
  real(dp), parameter :: ryd = 13.605698066_dp, ha2eV = 2.0_dp*ryd
  real(dp), parameter :: ha2cmm1 = 219474.63_dp  !-- multiply by this to convert from Hartree to cm-1
  real(dp), parameter :: K2Ha = 3.16682968067E-6_dp
  real(dp), parameter :: bohr2cm = 5.29177249E-9_dp
  real(dp), parameter :: time_au2sec = 2.418884326505E-17_dp  
  real(dp), parameter :: cmm3_2_au = 1.48184743477E-25_dp

  integer :: U_OUT
  integer :: nskip
  integer :: nb
  integer :: ivbm
  integer :: icbm
  integer :: nphonon
  integer :: vbmdeg
  real(dp) :: Egap_min
  real(dp) :: Egap_max
  integer :: ngap
  real(dp) :: kT
  integer :: ios
  character(LEN=50) :: eigFile
  character(LEN=50) :: elphFile
  character(LEN=50) :: phfreqFile

  real(dp) :: T_Kelvin
  real(dp) :: n_free
  real(dp) :: epsilon_infty
  real(dp) :: Vcell
  real(dp) :: n_val
  real(dp) :: a0
  real(dp) :: latVec11, latVec12, latVec13
  real(dp) :: latVec21, latVec22, latVec23
  real(dp) :: latVec31, latVec32, latVec33
 

  integer, parameter :: neta = 8
  real(dp) :: eta(neta)

  integer :: nk, iq, ieta, igap, inu, ideg, ib
  real(dp), allocatable :: kpoint(:,:) , w(:)
  real(dp) :: sumw, rdummy

  !-- electron energies and gaps
  real(dp), allocatable :: elda0(:,:)
  real(dp), allocatable :: elda(:,:)
  real(dp), allocatable  :: ebndG(:), ebndQ(:)
  real(dp) :: egap_lda, Delta_scissors
  real(dp), allocatable  :: Egap(:)

  integer :: nx, ny, nz, index
  integer :: idummy

  integer :: i, j !-- loop variables

  !-- wave functions
  complex(dp), allocatable :: wfnG(:,:), wfnq(:,:), wfnmq(:,:)

  !-- rate coefficients
#if  (defined EEH)
  integer :: maxdeg
  real(dp), parameter :: ehsign =  1.0_dp
#elif (defined HHE)
  integer :: maxdeg
  real(dp), parameter :: ehsign = -1.0_dp
#endif

  real(dp), allocatable :: auger_coef_ep_abs(:,:,:,:), auger_coef_ep_abs_global(:,:,:,:)
  real(dp), allocatable :: auger_coef_ep_emi(:,:,:,:), auger_coef_ep_emi_global(:,:,:,:)


  character(LEN=10) :: c
  character(LEN=100) :: filename
  character(LEN=100) :: sdummy

  real(dp), allocatable :: omega_q(:)
  real(dp) :: omega, odummy, x, y
  complex(dp), allocatable :: g_elphon(:,:,:)

  integer :: nb1, nb2, nb3, nb4, m

  complex(dp), allocatable :: auger_me_d_1(:),auger_me_d_2(:),auger_me_d_3(:),auger_me_d_4(:) !-- Auger ME, direct  , 4 terms
  complex(dp), allocatable  :: auger_me_x_1(:),auger_me_x_2(:),auger_me_x_3(:),auger_me_x_4(:) !-- Auger ME, exchange, 4 terms

  complex(dp) :: idelta
  complex(dp) :: ep_me_1, ep_me_2, ep_me_3, ep_me_4
  complex(dp) :: auger_gme_ep_d, auger_gme_ep_x
  real(dp) :: denom1, denom2, denom3, denom4

  !-- matrix elements squared and rate coefficients
  real(dp) :: auger_me2_ep_total, rate
  real(dp), allocatable :: ratebymode_abs(:)
  real(dp), allocatable :: ratebymode_emi(:)

  real :: time1, time2

  integer :: ident, nprocs
#ifdef MPI
  integer :: status(MPI_STATUS_SIZE)
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, ident,  ierr)
#else
  nprocs = 1
  ident = 0
#endif
 
    NAMELIST / auger / nskip, nb, ivbm, nphonon, vbmdeg, &
      Egap_min, Egap_max, ngap, T_Kelvin, eigFile, elphFile, phfreqFile, &
      Vcell, a0, n_free, n_val, epsilon_infty, &
      latVec11, latVec12, latVec13, latVec21, latVec22, latVec23, latVec31, latVec32, latVec33
    if (ident .eq. 0 ) then
     nskip = 0
     nb = 0
     ivbm = 0
     nphonon = 0
     vbmdeg = 1
     Egap_min = 0.1
     Egap_max = 1.0
     ngap = 91
     eigFile = 'eig.dat'
     elphFile = 'epmat1.ik01.dat'
     phfreqFile = 'omega'
     T_Kelvin = 300.0
     Vcell = 1.0
     a0 = 1.0
     latVec11 = 1.0
     latVec12 = 0.0
     latVec13 = 0.0
     latVec21 = 0.0
     latVec22 = 1.0
     latVec23 = 0.0
     latVec31 = 0.0
     latVec32 = 0.0
     latVec33 = 1.0
     n_val = 1.0
     n_free = 1E18
     epsilon_infty = 1.0
     
     CALL  input_from_file()

     READ (5, auger, iostat = ios)
     CLOSE(5)
     if  (ios .ne. 0) then
       write(6,*) 'Auger input ERRO IOS = ', ios
       STOP
     end if

     CALL freeunit (U_OUT)
     OPEN (U_OUT, FILE="auger.log", ACTION='WRITE',STATUS='REPLACE')
     write(U_OUT,*) '======PARAMETERS SETUP======'
     write(U_OUT,*) 'nskip = ', nskip
     write(U_OUT,*) 'nb = ', nb
     write(U_OUT,*) 'ivbm = ', ivbm
     write(U_OUT,*) 'nphonon = ', nphonon
     write(U_OUT,*) 'vbmdeg = ', vbmdeg
     write(U_OUT,*) 'Egap_min = ', Egap_min
     write(U_OUT,*) 'Egap_max = ', Egap_max
     write(U_OUT,*) 'ngap = ', ngap
     write(U_OUT,*) 'eigFile = ', eigFile
     write(U_OUT,*) 'elphFile = ', elphFile
     write(U_OUT,*) 'phfreqFile = ', phfreqFile
     write(U_OUT,*) 'T_Kelvin = ', T_Kelvin
     write(U_OUT,*) 'n_free = ', n_free
     write(U_OUT,*) 'n_val = ', n_val
     write(U_OUT,*) 'epsilon_infty = ', epsilon_infty
     write(U_OUT,*) 'Vcell = ', Vcell
     write(U_OUT,*) 'a0 = ', a0
     write(U_OUT,*) 'latVec11 = ', latVec11
     write(U_OUT,*) 'latVec12 = ', latVec12
     write(U_OUT,*) 'latVec13 = ', latVec13
     write(U_OUT,*) 'latVec21 = ', latVec21
     write(U_OUT,*) 'latVec22 = ', latVec22
     write(U_OUT,*) 'latVec23 = ', latVec23
     write(U_OUT,*) 'latVec31 = ', latVec31
     write(U_OUT,*) 'latVec32 = ', latVec32
     write(U_OUT,*) 'latVec33 = ', latVec33
     write(U_OUT,*) '============================'
     
     Egap_min = Egap_min/ha2eV
     Egap_max = Egap_max/ha2eV
     kT = T_Kelvin*K2Ha
     ivbm = ivbm - nskip
     icbm = ivbm + 1
     n_free = n_free * cmm3_2_au
     n_val = n_val / Vcell

#if (defined EEH)
     maxdeg = vbmdeg
#elif (defined HHE)
     maxdeg = vbmdeg ** 2
#endif
  
end if

#ifdef MPI
  CALL mpi_bcast( nskip ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr )
  CALL mpi_bcast( nb ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr )
  CALL mpi_bcast( ivbm ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr )
  CALL mpi_bcast( icbm ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr )
  CALL mpi_bcast( nphonon ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr )
  CALL mpi_bcast( vbmdeg ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr )
  CALL mpi_bcast( Egap_min ,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr )
  CALL mpi_bcast( Egap_max ,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr )
  CALL mpi_bcast( ngap ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr )
  CALL mpi_bcast( eigFile ,50,MPI_CHAR,0,MPI_COMM_WORLD,ierr )
  CALL mpi_bcast( elphFile ,50,MPI_CHAR,0,MPI_COMM_WORLD,ierr )
  CALL mpi_bcast( phfreqFile ,50,MPI_CHAR,0,MPI_COMM_WORLD,ierr )
  CALL mpi_bcast( kT ,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr )
  CALL mpi_bcast( maxdeg ,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr )
  CALL mpi_bcast( T_Kelvin ,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr )
  CALL mpi_bcast( a0 ,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr )
  CALL mpi_bcast( latVec11 ,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr )
  CALL mpi_bcast( latVec12 ,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr )
  CALL mpi_bcast( latVec13 ,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr )
  CALL mpi_bcast( latVec21 ,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr )
  CALL mpi_bcast( latVec22 ,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr )
  CALL mpi_bcast( latVec23 ,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr )
  CALL mpi_bcast( latVec31 ,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr )
  CALL mpi_bcast( latVec32 ,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr )
  CALL mpi_bcast( latVec33 ,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr )
  CALL mpi_bcast( n_free ,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr )
  CALL mpi_bcast( n_val ,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr )
  CALL mpi_bcast( epsilon_infty ,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr )
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

#endif

  allocate(ebndG(nb),ebndQ(nb))
  allocate(Egap(ngap))
  allocate(auger_coef_ep_abs(maxdeg,neta,ngap,nphonon),auger_coef_ep_abs_global(maxdeg,neta,ngap,nphonon)) 
  allocate(auger_coef_ep_emi(maxdeg,neta,ngap,nphonon),auger_coef_ep_emi_global(maxdeg,neta,ngap,nphonon))
  allocate(omega_q(nphonon))
  allocate(g_elphon(nphonon,nb,nb))
  allocate(auger_me_d_1(nb),auger_me_d_2(nb),auger_me_d_3(nb),auger_me_d_4(nb))
  allocate(auger_me_x_1(nb),auger_me_x_2(nb),auger_me_x_3(nb),auger_me_x_4(nb))
  allocate(ratebymode_abs(nphonon),ratebymode_emi(nphonon))

  eta = (/ 0.01_dp, 0.02_dp, 0.03_dp, 0.05_dp, &
       &   0.1_dp , 0.2_dp , 0.3_dp , 0.5_dp /)
  eta = eta/ha2eV
  
  do i = 1, ngap
     Egap(i) = Egap_min + (i-1) * (Egap_max-Egap_min) / real(ngap-1,dp)
  end do

  !-- q-kpoints: number, coordinates, and electron energies

  nk = file_row_count ( "klist_weights.dat" )

  allocate ( kpoint(3,nk), w(nk) )
  sumw = 0.0_dp
  open(UNIT=1,FILE="klist_weights.dat",ACTION='READ')
  do iq = 1, nk
     read(1,*) kpoint(:,iq), w(iq)
     sumw = sumw + w(iq)
  end do
  close(1)
  w = w / sumw  !-- Renormalize weights so that they sum to 1

  !-- this could be done differently and without this file and without allocate
  !-- if the wannier band energies are read into eldaG and elgaQ
  allocate ( elda0(nb,nk) )
  allocate ( elda(nb,nk) )

  do iq = 1, nk
   c = 'X'
   write(c,"(i0)") iq
   filename = "phonon_" // trim(c) // "/" // trim(eigFile)
   open(UNIT=1,FILE=filename,ACTION='READ')
   !- skip the first nb at gamma (in each directory)
   !- we can go back to reading a big list to port over VASP eigenvalues
   do ib = 1,nskip
    read(1,*) idummy, idummy, rdummy
   end do
   do ib = 1,nb-nskip
    read(1,*) idummy, idummy, elda0(ib,iq)
   end do
   !- then skip the nskip at q
   do ib = 1,nskip 
    read(1,*) idummy, idummy, rdummy
   end do
   do ib = 1,nb-nskip 
    read(1,*) idummy, idummy, elda(ib,iq) 
   end do
   close(1)
  end do
  elda0 = elda0/ha2eV
  elda = elda/ha2eV
  egap_lda = elda0(icbm,1) - elda0(ivbm,1)

  !-- wave functions: read (reduced) grid dimensions, allocate and set up FFT

  open(UNIT=1,FILE="phonon_1/UNK00001.1",FORM='UNFORMATTED',ACTION='READ')
  read(1) nx,ny,nz,idummy,idummy
  close(1)
  
  allocate(wfnG(nx*ny*nz,nb))
  allocate(wfnq(nx*ny*nz,nb))
  allocate(wfnmq(nx*ny*nz,nb))

  call initialize_me_fft(nx,ny,nz,a0,latVec11,latVec12,latVec13,latVec21,latVec22,latVec23,latVec31,latVec32,latVec33)

  !print *,'calculating matrix elements'
  auger_coef_ep_abs = 0.0_dp
  auger_coef_ep_emi = 0.0_dp
  
  do iq = 1+ident, nk, nprocs

     if ( mod(iq,1) .eq. 0 ) print *, ident, iq
     c = 'X'
     write(c,"(i0)") iq

     !PRINT*, 'Reading wfn at Gamma:'
     filename = "phonon_" // trim(c) // "/UNK00001.1"
     open(UNIT=1,FILE=filename,FORM='UNFORMATTED',ACTION='READ')
     read(1) idummy,idummy,idummy,idummy,idummy
     do m = 1, nb
        !PRINT*, 'band=', m
        read(1) (wfnG(index,m),index=1,nx*ny*nz)
     end do
     close(1)
     
     !PRINT*, 'Reading wfn at q:'
     filename = "phonon_" // trim(c) // "/UNK00002.1"
     open(UNIT=1,FILE=filename,FORM='UNFORMATTED',ACTION='READ')
     read(1) idummy,idummy,idummy,idummy,idummy
     do m = 1, nb
        read(1) (wfnq(index,m),index=1,nx*ny*nz)
     end do
     close(1)
     !PRINT*, 'DONE: Reading wfn at q'
     
     !PRINT*, 'Reading wfn at -q:'
     wfnmq = conjg(wfnq)

     !-- Reading phonon frequencies:
     filename = "phonon_" // trim(c) //"/"// trim(phfreqFile)
     !PRINT*, filename
     open(UNIT=1,FILE=filename,ACTION='READ')
     do inu = 1, nphonon
        read(1,*) omega_q(inu)
        !PRINT*, omega_q(inu,iq)
        if ( omega_q(inu) .LT. 0.0_dp ) omega_q(inu) = 0.0_dp
     end do
     close(1)
     omega_q = omega_q / ha2cmm1  !-- Convert from cm-1 to Ha

     !-- electron-phonon
     filename = "phonon_" // trim(c) //"/"// trim(elphFile)
     !PRINT*, filename
     open(UNIT=1, FILE=filename,ACTION='READ')
     read(1,*) ! to skip the 1st line
     do i = 1, nb+nskip
     do j = 1, nb+nskip
        do inu = 1, nphonon
           read(1,"(4I6,1F15.9,2E20.10)") idummy, idummy, idummy, idummy, odummy, x, y
           !-- Don't forget to convert to Ha and divide by SQRT(omega_q)!
           if ( i .GT. nskip .AND. j .GT. nskip ) then
              if ( omega_q(inu) .LT. 1.0E-10_dp ) then
                 g_elphon(inu,i-nskip,j-nskip) = CMPLX(0.0_dp, 0.0_dp)
              else
                 g_elphon(inu,i-nskip,j-nskip) = (CMPLX(x,y)/SQRT(omega_q(inu))/2.0_dp) 
              end if
           end if

        end do
     end do
     end do
     close(1)

     ebndG(1:ivbm) = ehsign * elda0(1:ivbm, 1)
     ebndQ(1:ivbm) = ehsign * elda(1:ivbm,iq)

     call cpu_time(time1)

#if (defined EEH)
     nb1 = icbm
     nb2 = icbm
     ideg = 0
#ifdef SO
    do nb3 = ivbm - 2, ivbm, 2
#else
     do nb3 = ivbm - vbmdeg + 1, ivbm
#endif
     ideg = ideg + 1
     do nb4 = icbm, nb
#elif (defined HHE)
     ideg = 0
#ifdef SO
     do nb1 = ivbm - 2, ivbm, 2
     do nb2 = ivbm - 2, ivbm, 2
#else
     do nb1 = ivbm - vbmdeg + 1, ivbm
     do nb2 = ivbm - vbmdeg + 1, ivbm
#endif
     ideg = ideg + 1
     nb3 = icbm
     do nb4 = 1, ivbm
#endif

        !-- Calculate Auger MEs only once, not for every phonon mode iteration:
        do m = 1, nb
           auger_me_d_1(m) = coulomb_me_fft( wfnq(:, m ), wfnG(:,nb2), wfnG(:,nb3), wfnq(:,nb4), kpoint(:,iq),T_Kelvin,n_val,n_free,epsilon_infty )
           auger_me_d_2(m) = coulomb_me_fft( wfnG(:,nb1), wfnq(:,m  ), wfnG(:,nb3), wfnq(:,nb4), kpoint(:, 1),T_Kelvin,n_val,n_free,epsilon_infty )
           auger_me_d_3(m) = coulomb_me_fft( wfnG(:,nb1), wfnG(:,nb2), wfnmq(:, m), wfnq(:,nb4), kpoint(:,iq),T_Kelvin,n_val,n_free,epsilon_infty )
           auger_me_d_4(m) = coulomb_me_fft( wfnG(:,nb1), wfnG(:,nb2), wfnG(:,nb3), wfnG(:, m ), kpoint(:, 1),T_Kelvin,n_val,n_free,epsilon_infty )

           auger_me_x_1(m) = coulomb_me_fft( wfnG(:,nb2), wfnq(:, m ), wfnG(:,nb3), wfnq(:,nb4), kpoint(:, 1),T_Kelvin,n_val,n_free,epsilon_infty )
           auger_me_x_2(m) = coulomb_me_fft( wfnq(:, m ), wfnG(:,nb1), wfnG(:,nb3), wfnq(:,nb4), kpoint(:,iq),T_Kelvin,n_val,n_free,epsilon_infty )
           auger_me_x_3(m) = coulomb_me_fft( wfnG(:,nb2), wfnG(:,nb1), wfnmq(:, m), wfnq(:,nb4), kpoint(:,iq),T_Kelvin,n_val,n_free,epsilon_infty )
           auger_me_x_4(m) = coulomb_me_fft( wfnG(:,nb2), wfnG(:,nb1), wfnG(:,nb3), wfnG(:, m ), kpoint(:, 1),T_Kelvin,n_val,n_free,epsilon_infty )

        end do !m = 1, nb
        
        do inu = 1, nphonon
           if ( omega_q(inu) .lt. 1.0E-10_dp ) cycle !-- pathological phonons do not contribute (already excluded above)
           omega = omega_q(inu)

           do igap = 1, ngap
              Delta_scissors = Egap(igap) - egap_lda !-- The correction needed to the LDA eigenvalues to achieve the gap Egap(igap)
              ebndG(icbm:nb) = ehsign * ( elda0(icbm:nb, 1) + Delta_scissors )
              ebndQ(icbm:nb) = ehsign * ( elda(icbm:nb,iq) + Delta_scissors )

              do ieta = 1, neta
                 idelta = ehsign * CMPLX(0.0_dp,eta(ieta))

                 !-- caveat: the following is written for given ieta, inu, igap. These indeces are dropped
                 !-- m over all bands has to be the inner loop!

                 !-- absorption of +q phonon
                 auger_gme_ep_d = CMPLX(0.0_dp,0.0_dp)
                 auger_gme_ep_x = CMPLX(0.0_dp,0.0_dp)
                 
                 do m = 1, nb
                    !-- Electron-phonon scattering
                    ep_me_1 = g_elphon(inu,nb1,m)
                    ep_me_2 = g_elphon(inu,nb2,m)
                    ep_me_3 = conjg(g_elphon(inu,nb3,m))
                    ep_me_4 = g_elphon(inu,m,nb4)
                    
                    denom1 = ebndQ(m) - ebndG(nb1) - omega
                    denom2 = ebndQ(m) - ebndG(nb2) - omega
                    denom3 = ebndQ(m) - ebndG(nb3) + omega
                    denom4 = ebndG(m) - ebndQ(nb4) + omega
                    
                    auger_gme_ep_d = auger_gme_ep_d + &
                         & ep_me_1*auger_me_d_1(m) / (denom1+idelta) + &
                         & ep_me_2*auger_me_d_2(m) / (denom2+idelta) + &
                         & auger_me_d_3(m)*ep_me_3 / (denom3+idelta) + &
                         & auger_me_d_4(m)*ep_me_4 / (denom4+idelta)
                    auger_gme_ep_x = auger_gme_ep_x + &
                         & ep_me_1*auger_me_x_1(m) / (denom1+idelta) + &
                         & ep_me_2*auger_me_x_2(m) / (denom2+idelta) + &
                         & auger_me_x_3(m)*ep_me_3 / (denom3+idelta) + &
                         & auger_me_x_4(m)*ep_me_4 / (denom4+idelta)
                 end do
                 
                 auger_me2_ep_total = ABS(auger_gme_ep_d - auger_gme_ep_x)**2 + ABS(auger_gme_ep_d)**2 + ABS(auger_gme_ep_x)**2 
                 
                 rate =  w(iq) * auger_me2_ep_total * ( 1.0_dp / (exp(omega/kT)-1.0_dp) ) * &
                      & EXP( -(ebndG(nb1)+ebndG(nb2)+omega-ebndG(nb3)-ebndQ(nb4))**2/eta(ieta)**2 )/SQRT(pi*eta(ieta)**2)

                 auger_coef_ep_abs(ideg,ieta,igap,inu) = auger_coef_ep_abs(ideg,ieta,igap,inu) + rate

                 !-- emission of -q phonon
                 !-- omega -> -omega, bose -> bose + 1
                 auger_gme_ep_d = CMPLX(0.0_dp,0.0_dp)
                 auger_gme_ep_x = CMPLX(0.0_dp,0.0_dp)
                 
                 do m = 1, nb
                    !-- Electron-phonon scattering
                    ep_me_1 = g_elphon(inu,nb1,m)
                    ep_me_2 = g_elphon(inu,nb2,m)
                    ep_me_3 = conjg(g_elphon(inu,nb3,m))
                    ep_me_4 = g_elphon(inu,m,nb4)
                    
                    denom1 = ebndQ(m) - ebndG(nb1) + omega
                    denom2 = ebndQ(m) - ebndG(nb2) + omega
                    denom3 = ebndQ(m) - ebndG(nb3) - omega
                    denom4 = ebndG(m) - ebndQ(nb4) - omega

                    auger_gme_ep_d = auger_gme_ep_d + &
                         & ep_me_1*auger_me_d_1(m) / (denom1+idelta) + &
                         & ep_me_2*auger_me_d_2(m) / (denom2+idelta) + &
                         & auger_me_d_3(m)*ep_me_3 / (denom3+idelta) + &
                         & auger_me_d_4(m)*ep_me_4 / (denom4+idelta)
                    auger_gme_ep_x = auger_gme_ep_x + &
                         & ep_me_1*auger_me_x_1(m) / (denom1+idelta) + &
                         & ep_me_2*auger_me_x_2(m) / (denom2+idelta) + &
                         & auger_me_x_3(m)*ep_me_3 / (denom3+idelta) + &
                         & auger_me_x_4(m)*ep_me_4 / (denom4+idelta)
                 end do
                 
                 auger_me2_ep_total = ABS(auger_gme_ep_d - auger_gme_ep_x)**2 + ABS(auger_gme_ep_d)**2 + ABS(auger_gme_ep_x)**2 

                 rate = w(iq) * auger_me2_ep_total * ( 1.0_dp / (exp(omega/kT)-1.0_dp) + 1.0_dp) * &
                      & EXP( -(ebndG(nb1)+ebndG(nb2)-omega-ebndG(nb3)-ebndQ(nb4))**2/eta(ieta)**2 )/SQRT(pi*eta(ieta)**2)
                 
                 auger_coef_ep_emi(ideg,ieta,igap,inu) = auger_coef_ep_emi(ideg,ieta,igap,inu) + rate

              end do ! ieta = 1, neta
           end do ! igap = 1, ngap
        end do ! inu = 1, nphonon

        ! now we have the rate contribution for all
        ! phonon modes, smearing widths, band gaps
        ! for given bands and q-point

     end do !bands
     end do !bands
#ifdef HHE
     end do !bands
#endif

     call cpu_time(time2)
     write (*,*) "timing", time2 - time1

  end do !iq nk phonons

  auger_coef_ep_abs_global = 0.0_dp
  auger_coef_ep_emi_global = 0.0_dp
#ifdef MPI
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  call MPI_REDUCE( auger_coef_ep_abs, auger_coef_ep_abs_global, maxdeg*neta*ngap*nphonon, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  call MPI_REDUCE( auger_coef_ep_emi, auger_coef_ep_emi_global, maxdeg*neta*ngap*nphonon, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#else
  auger_coef_ep_abs_global = auger_coef_ep_abs
  auger_coef_ep_emi_global = auger_coef_ep_emi
#endif
  
  auger_coef_ep_abs_global = auger_coef_ep_abs_global * (pi / 2.0_dp) * (bohr2cm)**6 / time_au2sec
  auger_coef_ep_emi_global = auger_coef_ep_emi_global * (pi / 2.0_dp) * (bohr2cm)**6 / time_au2sec

  if ( ident .EQ. 0 ) then
     do ieta = 1, neta
        write(c,"(i4.4)") ieta

#if (defined EEH)
        sdummy = 'eeh'
#elif (defined HHE)
        sdummy = 'hhe'
#endif
        filename = "auger_coef_ep_" // trim(sdummy) // "_absorption_" // trim(c) // ".dat"
        open(UNIT=11,FILE=filename)
        filename = "auger_coef_ep_" // trim(sdummy) // "_emission_"   // trim(c) // ".dat"
        open(UNIT=12,FILE=filename)
        write(11,*) "# Egap(eV), Auger_coef(cm6s-1), eta (eV) = ", eta(ieta)*ha2eV
        write(12,*) "# Egap(eV), Auger_coef(cm6s-1), eta (eV) = ", eta(ieta)*ha2eV
        do igap = 1, ngap
           ratebymode_abs = 0._dp
           ratebymode_emi = 0._dp
           do inu = 1, nphonon
              do ideg = 1, maxdeg
                 ratebymode_abs(inu) = ratebymode_abs(inu) + auger_coef_ep_abs_global(ideg,ieta,igap,inu)
                 ratebymode_emi(inu) = ratebymode_emi(inu) + auger_coef_ep_emi_global(ideg,ieta,igap,inu)
              end do
           end do
           write(11,"(f7.3,24g14.6)") Egap(igap)*ha2eV, ratebymode_abs(1:nphonon) / real(maxdeg,dp)
           write(12,"(f7.3,24g14.6)") Egap(igap)*ha2eV, ratebymode_emi(1:nphonon) / real(maxdeg,dp)
        end do
        close(11)
        close(12)

     end do

  end if

#ifdef MPI
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  call MPI_FINALIZE(ierr)
#endif
  stop
end program phononauger
