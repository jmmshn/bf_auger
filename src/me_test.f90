program me_test
  USE kinds, ONLY : DP
  USE base

  implicit none
  save
  integer :: nk, nb, nr1, nr2, nr3, ios
  integer :: np 
  real :: start1, finish1, start2, finish2

  integer :: ii, fu
  real(dp),allocatable :: elda(:,:), elda_k(:)
  complex(dp),allocatable :: wfn(:,:), wfn2(:,:)

! ============================================================
! Get the user inputs
! Using NAMELIST and routines taken from QE
!
! ============================================================
  NAMELIST / me_inp / nk, nb, nr1, nr2, nr3
! ========== Defaults
  nk = 0 ; nb = 0 ; nr1 = 0 ; nr2 = 0 ; nr3 = 0; 
! ========== Read in 
  CALL input_from_file ( )
! 
  READ (5, me_inp, iostat = ios)
! 
  write(6,*) '======== Parameters ME_INP =========', nk
  write(6,*) 'nk=', nk
  write(6,*) 'nb=', nb

  write(6,'(A20,3i4)') 'nr1,  nr2,  nr3  = ', nr1, nr2, nr3
  close(5)


! ========== Declarations
  np = nr1*nr2*nr3
  allocate(elda_k(nb))
  allocate(wfn(np,nb))
  allocate(wfn2(np/8,nb))
! ============================================================
! Compare the UNKs and WFN_Rs
! ============================================================

! ========== Defaults
  call freeunit (fu)
  write(6,*) 'START elda'
  call read_elda(1,nb,elda_k,fu)
  write(6,*) 'eigs  = '
  write(6,*)(elda_k(ii), ii=1,nb)
 
  write(6,*) 'UNK'
  call cpu_time(start1)
  !do ii = 1, 10000
  call read_wfn ( 1, np/8, nb, wfn2, fu)
  !end do
  call cpu_time(finish1)
  
  write(6,*) 'WFN'
  call cpu_time(start2)
  !do ii = 1, 10000
  call read_wfn_fft ( 1, np, nb, wfn, fu)
  !end do
  call cpu_time(finish2)
  
  write(6,*) 'UNK times:', finish1-start1
  write(6,*) 'WFN times:', finish2-start2
  !call setup_fft (nx, ny, nz)

end program me_test
