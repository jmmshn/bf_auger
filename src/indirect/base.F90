module base
  implicit none
  save 
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
 
  end subroutine freeunit

end module base
