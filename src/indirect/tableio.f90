module tableio
  implicit none
  private
  public :: file_column_count, file_row_count

contains
  !*****************************************************************************80
  !
  !! FILE_COLUMN_COUNT counts the number of columns in the first line of a file.
  !
  !  Discussion:
  !
  !    The file is assumed to be a simple text file.
  !
  !    Most lines of the file is presumed to consist of COLUMN_NUM words,
  !    separated by spaces.  There may also be some blank lines, and some
  !    comment lines,
  !    which have a "#" in column 1.
  !
  !    The routine tries to find the first non-comment non-blank line and
  !    counts the number of words in that line.
  !
  !    If all lines are blanks or comments, it goes back and tries to analyze
  !    a comment line.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    21 June 2001
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) INPUT_FILENAME, the name of the file.
  !
  !    Output, integer :: COLUMN_NUM, the number of columns in the file.
  !
  integer function file_column_count ( input_filename )
    
    implicit none
    
    integer :: column_num
    logical got_one
    character ( len = * ) input_filename
    integer :: input_status
    integer :: input_unit
    character ( len = 255 ) line
    !
    !  Open the file.
    !
    call get_unit ( input_unit )
    
    open ( unit = input_unit, file = input_filename, status = 'old', &
         form = 'formatted', access = 'sequential', iostat = input_status )
    
    if ( input_status /= 0 ) then
       column_num = -1
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Fatal error!'
       write ( *, '(a,i8)' ) '  Could not open the input file "' &
            // trim ( input_filename ) // '" on unit ', input_unit
       return
    end if
    !
    !  Read one line, but skip blank lines and comment lines.
    !
    got_one = .false.
    
    do
       
       read ( input_unit, '(a)', iostat = input_status ) line
       
       if ( input_status /= 0 ) then
          exit
       end if
       
       if ( len_trim ( line ) == 0 ) then
          cycle
       end if
       
       if ( line(1:1) == '#' ) then
          cycle
       end if
       
       got_one = .true.
       exit
       
    end do
    
    if ( .not. got_one ) then
       
       rewind ( input_unit )
       
       do
          
          read ( input_unit, '(a)', iostat = input_status ) line
          
          if ( input_status /= 0 ) then
             exit
          end if
          
          if ( len_trim ( line ) == 0 ) then
             cycle
          end if
          
          got_one = .true.
          exit
          
       end do
       
    end if
    
    close ( unit = input_unit )
    
    if ( .not. got_one ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Warning!'
       write ( *, '(a)' ) '  The file does not seem to contain any data.'
       column_num = -1
       return
    end if
    
    call s_word_count ( line, column_num )
    
    file_column_count = column_num
    return
  end function file_column_count
  
  !*****************************************************************************80
  !
  !! FILE_ROW_COUNT counts the number of row records in a file.
  !
  !  Discussion:
  !
  !    It does not count lines that are blank, or that begin with a
  !    comment symbol '#'.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    06 March 2003
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
  !
  !    Output, integer :: ROW_NUM, the number of rows found.
  !
  integer function file_row_count ( input_filename )
    
    implicit none
    
    character ( len = * ) input_filename
    character ( len = 255 ) line
    integer :: input_status
    integer :: input_unit
    integer :: record_num
    integer :: comment_num
    integer :: row_num
    
    call get_unit ( input_unit )
    
    open ( unit = input_unit, file = input_filename, status = 'old', &
         iostat = input_status )
    
    if ( input_status /= 0 ) then
       row_num = -1;
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'FILE_ROW_COUNT - Fatal error!'
       write ( *, '(a,i8)' ) '  Could not open the input file "' // &
            trim ( input_filename ) // '" on unit ', input_unit
       stop
    end if
    
    comment_num = 0
    row_num = 0
    record_num = 0
    
    do
       
       read ( input_unit, '(a)', iostat = input_status ) line
       
       if ( input_status /= 0 ) then
          exit
       end if
       
       record_num = record_num + 1
       
       if ( line(1:1) == '#' ) then
          comment_num = comment_num + 1
          cycle
       end if
       
       if ( len_trim ( line ) == 0 ) then
          comment_num = comment_num + 1
          cycle
       end if
       
       row_num = row_num + 1
       
    end do
    
    close ( unit = input_unit )
    
    file_row_count = row_num
    return
  end function file_row_count

  !*****************************************************************************80
  !
  !! GET_UNIT returns a free FORTRAN unit number.
  !
  !  Discussion:
  !
  !    A "free" FORTRAN unit number is a value between 1 and 99 which
  !    is not currently associated with an I/O device.  A free FORTRAN unit
  !    number is needed in order to open a file with the OPEN command.
  !
  !    If IUNIT = 0, then no free FORTRAN unit could be found, although
  !    all 99 units were checked (except for units 5, 6 and 9, which
  !    are commonly reserved for console I/O).
  !
  !    Otherwise, IUNIT is a value between 10 and 99, representing a
  !    free FORTRAN unit.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    26 October 2008
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Output, integer :: IUNIT, the free unit number.
  !
  subroutine get_unit ( iunit )
    
    implicit none
    
    integer :: i
    integer :: ios
    integer :: iunit
    logical lopen
    
    iunit = 0
    
    do i = 10, 99
       
       inquire ( unit = i, opened = lopen, iostat = ios )
       
       if ( ios == 0 ) then
          if ( .not. lopen ) then
             iunit = i
             return
          end if
       end if
       
    end do
    
    return
  end subroutine get_unit

  !*****************************************************************************80
  !
  !! S_WORD_COUNT counts the number of "words" in a string.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string to be examined.
  !
  !    Output, integer :: NWORD, the number of "words" in the string.
  !    Words are presumed to be separated by one or more blanks.
  !
  subroutine s_word_count ( s, nword )
    
    implicit none
    
    logical blank
    integer :: i
    integer :: lens
    integer :: nword
    character ( len = * )  s
    
    nword = 0
    lens = len ( s )
    
    if ( lens <= 0 ) then
       return
    end if
    
    blank = .true.
    
    do i = 1, lens
       
       if ( s(i:i) == ' ' ) then
          blank = .true.
       else if ( blank ) then
          nword = nword + 1
          blank = .false.
       end if
       
    end do
    
    return
  end subroutine s_word_count

end module tableio
