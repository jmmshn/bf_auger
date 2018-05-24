module base
  implicit none
  save

 !integer, parameter :: dp = SELECTED_REAL_KIND(12,256) ! Manos
 !integer, parameter :: dp = SELECTED_REAL_KIND(14,200) ! QE
  integer, parameter :: dp = SELECTED_REAL_KIND(15,300) ! W90
  integer, parameter :: r8 = dp
  logical :: lsoc  ! are we doing spin-orbit

  REAL(dp), PARAMETER :: pi    = 3.14159265358979323846264338327950288_dp
  real(r8), parameter :: ryd   = 13.605698066_r8
  real(r8), parameter :: ha2eV = 2.0_r8 * ryd
 !REAL(dp), PARAMETER :: ha2eV = 27.2113961318_dp

  ! these are very confusing but I think everything is right in the code
  real(r8), parameter :: bohr2cm = 5.29177249E-9_r8 ! this is actually bohr_per_cm and all calculations should be right
  REAL(dp), PARAMETER :: cmm3tobohrm3 = 1.48184743477E-25_dp ! this is "bohr2cm" ** 3

  real(r8), parameter :: kelvin2au = 3.16682968067E-6_r8
  real(r8), parameter :: K2Ha      = Kelvin2au
  REAL(dp), PARAMETER :: bohr2Ang = 0.5291772108_dp
  real(dp), parameter :: time_au2sec = 2.418884326505E-17_dp
  integer :: ivbm , icbm , nr1, nr2, nr3, nb
  character(len=50) :: cEIGfile
  
  public

contains
  
  subroutine freeunit(ifu)
    implicit none

    integer :: ifu, ios
    logical :: isopen

    ios = 1

    ifu = 10
    do while ( isopen .or. (ios.ne.0) )
       ifu = ifu + 1
       inquire (unit=ifu, opened=isopen, iostat=ios)
    end do

    !print *, 'found free unit ', ifu
    
  end subroutine freeunit

  subroutine read_elda (nk, nb, elda, fu)
    implicit none

    integer, intent(in) :: nk, nb
    real(r8) :: elda(nb)
    integer, optional, intent(in) :: fu

    character(len=20) :: c
    integer :: ib, idummy
    
    if (.not.present(fu)) then
       call freeunit(fu)
    end if
  
    c = 'X'
    write(c,"(i0)") nk
    open(unit=fu, file='kpoints/k_'//trim(c)//'/'//trim(cEIGfile), action='READ')
    do ib = 1, nb
       read(fu,*) idummy, idummy, elda(ib)
    end do
    close(fu)

    elda = elda / Ha2eV

  end subroutine read_elda

  subroutine read_wfn ( nk, np, nb, wfn, fu)
    implicit none
    
    integer, intent(in) :: nk, nb, np
    complex(dp) :: wfn(np, nb)
    integer, optional, intent(in) :: fu

    character(len=20) :: c
    integer :: ib, j, nx, ny, nz, idummy, nb1
    
    if (.not.present(fu)) then
       call freeunit(fu)
    end if

    c = 'X'
    write(c,"(i0)") nk
    if (lsoc .eq. .False.) then
      OPEN (UNIT=fu, FILE="kpoints/k_"//TRIM(c)//"/UNK00001.1", FORM='UNFORMATTED', ACTION='READ')
    else
      OPEN (UNIT=fu, FILE="kpoints/k_"//TRIM(c)//"/UNK00001.SOC.1", FORM='UNFORMATTED', ACTION='READ')
    endif 

    READ (fu) nx, ny, nz, idummy, nb1
    if (np.ne.nx*ny*nz) then
      write (*,*) 'ERROR: grid in read_wfn'
      stop
    end if
    if (nb.ne.nb1) then
       write (*,*) 'ERROR: band in read_wfn'
       stop
    end if
    DO ib = 1, nb
       READ (fu) (wfn(j,ib), j = 1, nx*ny*nz)
    END DO
    CLOSE(fu)

  end subroutine read_wfn

end module base
