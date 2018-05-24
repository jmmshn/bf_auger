module kgrids
  use base, only: dp
  implicit none
  save

  integer :: nk1, nk2, nk3, nktot, nktot_h
  real(dp) :: alat, omega, at(3,3), bg(3,3)

  integer, allocatable :: equiv_ee(:), wk_e(:)
  integer, allocatable :: equiv_hh(:), wk_h(:)
  integer, allocatable :: map_he(:)
  real(dp), allocatable :: kgbz(:,:)

  private
  public :: nktot, equiv_ee, wk_e, kgbz, &
       nktot_h, equiv_hh, wk_h, map_he, &
       setup_kgrid, create_coord_anew, find_knum, firstbz, &
       omega, bg

contains

  subroutine setup_kgrid
    use base, only: freeunit, pi
    implicit none

    integer :: nk
    integer :: n, i, j, k, ncheck
    integer :: nk1h, nk2h, nk3h
    real(dp) :: xkg(3), xkc(3), xkgh(3)
    integer :: G_Umklapp(3)

    integer :: fu

    write (*,"(a)") 'setting up k-point grids'

    call freeunit(fu)
    open(unit=fu, file='kpoints.elec')

    read(fu,*) nk1, nk2, nk3
    nktot = nk1 * nk2 * nk3
    write(*,"(a,3i4)") 'elec grid: ', nk1, nk2, nk3

    read(fu,*) alat, omega
    do i = 1, 3
       read(fu,*) at(1:3,i)
    end do
    do i = 1, 3
       read(fu,*) bg(1:3,i)
    end do

    at = at * alat
    bg = bg * 2._dp * pi / alat
    
    allocate(equiv_ee(nktot))
    allocate(wk_e(nktot))
    allocate(kgbz(3,nktot))

    do nk = 1, nktot
       read(fu,*) n, equiv_ee(n), i, j, k, xkg(1:3), xkc(1:3), wk_e(n)
       call firstbz(xkg(1:3), kgbz(1:3,n), G_Umklapp(1:3))
       !-- paranoid check
       call find_knum(ncheck, kgbz(1:3,n))
       if (ncheck .ne. n) then
          stop
       end if
    end do

    close(fu)

    open(unit=fu, file='kpoints.hole')
    read(fu,*) nk1h, nk2h, nk3h
    nktot_h = nk1h * nk2h * nk3h
    write(*,"(a,3i4)") 'hole grid: ', nk1h, nk2h, nk3h

    read(fu,*) alat, omega
    do i = 1, 3
       read(fu,*) at(1:3,i)
    end do
    do i = 1, 3
       read(fu,*) bg(1:3,i)
    end do
 
    at = at * alat
    bg = bg * 2._dp * pi / alat

    allocate(equiv_hh(nktot_h))
    allocate(wk_h(nktot_h))
    allocate(map_he(nktot_h))

    do nk = 1, nktot_h
       read(fu,*) n, equiv_hh(n), i, j, k, xkgh(1:3), xkc(1:3), wk_h(n)
       call find_knum(ncheck, xkgh(1:3))
       map_he(n) = ncheck !-- map the kpt_hole index onto the kpt_elec index
    end do

    close(fu)

  end subroutine setup_kgrid

  subroutine create_coord_anew
    implicit none
    integer :: i, j, k, l, nk
    real(dp) :: xkg(3,nktot)
    integer :: G_Umklapp(3)

    !-- this is verbatim from QE
    do i = 1, nk1
    do j = 1, nk2
    do k = 1, nk3
       !  this is nothing but consecutive ordering
       nk = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
       !  xkg are the components of the complete grid in crystal axis
       xkg(1,nk) = DBLE(i-1)/nk1
       xkg(2,nk) = DBLE(j-1)/nk2
       xkg(3,nk) = DBLE(k-1)/nk3
       !  bring back into to the first BZ
       do l = 1, 3
          xkg(l,nk) = xkg(l,nk) - nint(xkg(l,nk))
       end do
       call firstbz(xkg(1:3,nk), kgbz(1:3,nk), G_Umklapp(1:3))
    end do
    end do
    end do

  end subroutine create_coord_anew

  subroutine find_knum(nk, krecp)
    implicit none

    integer :: nk
    real(dp) :: krecp(3)

    integer :: i, j, k
    
    i = mod(9 * nk1 + nint(krecp(1) * nk1), nk1)
    j = mod(9 * nk2 + nint(krecp(2) * nk2), nk2)
    k = mod(9 * nk3 + nint(krecp(3) * nk3), nk3)

    nk = k + j * nk3 + i * nk2 * nk3 + 1
           
  end subroutine find_knum

  subroutine firstbz(krecp, kgbz, G_Umklapp)
    implicit none
    
    real(dp), intent(in) :: krecp(3)
    real(dp) :: kgbz(3)
    integer :: G_Umklapp(3)

    real(dp) :: qmin2, ktry(3), qtry2
    integer :: j1, j2, j3

    kgbz = krecp
    call len2_krecp(qmin2, krecp)
    G_Umklapp = (/0,0,0/)

    do j1 = -3, 3
    do j2 = -3, 3
    do j3 = -3, 3
       ktry(:) = (/ j1, j2, j3 /) + krecp(:)
       call len2_krecp(qtry2, ktry)
       if ( qtry2 .LT. qmin2 ) then
          qmin2 = qtry2
          kgbz = ktry
          G_Umklapp = (/ j1, j2, j3 /)
       end if
    end do
    end do
    end do

  end subroutine firstbz

  subroutine len2_krecp(len2, krecp)
    implicit none
    real(dp) :: len2, krecp(3), kcart(3)
    kcart = MATMUL(bg,krecp) !-change from reciprocal to cartesian
    len2 = DOT_PRODUCT(kcart, kcart)    !-returns the square of the length
                                        !-as the first argument
  end subroutine len2_krecp
  
end module kgrids
