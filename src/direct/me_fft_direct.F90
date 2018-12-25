module me_fft_direct
  use base, only: dp

  implicit none
  save

  include 'fftw3.f'
  integer(8) :: plan_b, plan_f

  complex(dp), allocatable :: fftw_in_b(:,:,:), fftw_out_b(:,:,:)
  complex(dp), allocatable :: fftw_in_f(:,:,:), fftw_out_f(:,:,:)
  complex(dp), allocatable :: I13(:,:,:), I24(:,:,:)

  integer :: nx, ny, nz, np
  integer, allocatable :: Gmin_index(:,:,:,:)

  !-- Dielectric function parameters:
  real(dp) :: q2_TF
  real(dp) :: omega2_p
  real(dp) :: lambda2
  
  private
  public :: setup_fft, setup_epsilon, coulomb_me_fft

contains

  subroutine setup_fft (nx1, ny1, nz1)
    use base,   only: pi
    use kgrids, only: B => bg
    implicit none

    integer, intent(in) :: nx1, ny1, nz1

    integer :: i, j, k, i1, i2, i3
    real(dp) :: G(3), Gmin_norm2, Gnorm2

    nx = nx1
    ny = ny1
    nz = nz1

    ALLOCATE (fftw_in_f(nx,ny,nz),fftw_out_f(nx,ny,nz))
    ALLOCATE (fftw_in_b(nx,ny,nz),fftw_out_b(nx,ny,nz))
    call dfftw_plan_dft_3d (plan_b,nx,ny,nz,fftw_in_b,fftw_out_b,FFTW_BACKWARD,FFTW_MEASURE)
    call dfftw_plan_dft_3d (plan_f,nx,ny,nz,fftw_in_f,fftw_out_f,FFTW_FORWARD, FFTW_MEASURE)
    
    ALLOCATE (I13(nx,ny,nz),I24(nx,ny,nz))

    ALLOCATE (Gmin_index(3,nx,ny,nz))
    !IF ( ident .EQ. 0 ) PRINT *, 'Generating G min indices'
    do k = 1, nz ; do j = 1, ny ; do i = 1, nx
       !-- Map G vector to Wigner-Seitz supercell
       G(:) = (i-1)*b(:,1) + (j-1)*b(:,2) + (k-1)*b(:,3)
       Gnorm2 = DOT_PRODUCT(G,G)
       Gmin_index(:,i,j,k) = (/ i-1, j-1, k-1 /)
       Gmin_norm2 = Gnorm2
       do i1 = -1, 1
          do i2 = -1, 1
             do i3 = -1, 1
                G(:) = (i-1+i1*nx)*b(:,1) + (j-1+i2*ny)*b(:,2) + (k-1+i3*nz)*b(:,3)
                Gnorm2 = DOT_PRODUCT(G,G)
                IF ( Gnorm2 .LT. Gmin_norm2 ) THEN
                   Gmin_index(:,i,j,k) = (/ i-1+i1*nx, j-1+i2*ny, k-1+i3*nz /)
                   Gmin_norm2 = Gnorm2
                END IF
             end do
          end do
       end do
       !WRITE(*,"6I4") i, j, k, Gmin_index(:,i,j,k)
    end do; end do ; end do
    !IF ( ident .EQ. 0 ) PRINT *, 'Generating G min indices: DONE'
    
  end subroutine setup_fft

  subroutine setup_epsilon (n_free, kT)
    use base,   only: pi
    use kgrids, only: Vcell => omega
    implicit none

    real(dp), intent(in) :: n_free, kT

    real(dp) :: n_val

!#if defined (GAAS) || defined (INAS) || defined (INP)
    n_val = real(3+5, dp)/Vcell !GaAs
!#endif

    q2_TF = 4.0_dp * ( 3.0_dp/pi * n_val ) ** (1.0_dp/3.0_dp)
    omega2_p = 4.0_dp * pi * n_val
    lambda2 = 4.0_dp * pi * n_free / kT

  end subroutine setup_epsilon

  complex(dp) function coulomb_me_fft (wfn1, wfn2, wfn3, wfn4, G_Umklapp, mtrans)
    use base,   only: pi
    use kgrids, only: B => bg
    implicit none
    complex(dp) :: wfn1(:), wfn2(:), wfn3(:), wfn4(:)
    real(dp)    :: mtrans(:) !-- Coulomb momentum transfer
    integer     :: G_Umklapp(:)

    complex(dp), parameter :: iimag = cmplx(0.d0, 1.d0)

    real(dp) :: epsilon
    !-- Dielectric function parameters:
    real(dp), parameter :: alpha = 1.563_dp
!#ifdef GAAS
    real(dp), parameter :: epsilon_0 = 10.89_dp !GaAs experimental, RPA 8.9 (see Bechstedt's paper)
!#endif

    integer :: index, i, j, k
    real(dp) :: Gtemp(3), Gnorm2
    complex(dp) :: factor, me

    !-- Auger ME:
    index = 0
    do k = 1, nz ; do j = 1, ny ; do i = 1, nx
       index = index + 1
       factor = EXP( 2.0_dp * pi * iimag * &
            DOT_PRODUCT( real(G_Umklapp,dp), &
            (/ real(i-1,dp)/nx, real(j-1,dp)/ny, real(k-1,dp)/nz /) ) )
       I13(i,j,k) = CONJG(wfn1(index)) * wfn3(index) / (nx*ny*nz)
       I24(i,j,k) = CONJG(wfn2(index)) * wfn4(index) * factor / (nx*ny*nz)
    end do; end do ; end do

    fftw_in_b = I13 ; call dfftw_execute(plan_b) ; I13 = fftw_out_b
    fftw_in_f = I24 ; call dfftw_execute(plan_f) ; I24 = fftw_out_f
        
    me = CMPLX(0.0_dp,0.0_dp)
    do k = 1, nz ; do j = 1, ny ; do i = 1, nx
       Gtemp = MATMUL( b, Gmin_index(:,i,j,k) + mtrans(:) )
       Gnorm2 = DOT_PRODUCT(Gtemp,Gtemp)
       epsilon = 1.0_dp + 1.0_dp/( 1.0_dp/(epsilon_0-1.0_dp) + alpha*Gnorm2/q2_TF + Gnorm2**2/(4.0_dp*omega2_p) ) 
       me = me + (1.0_dp/epsilon) * I13(i,j,k) * I24(i,j,k) / (Gnorm2 + lambda2)
    end do; end do ; end do

    coulomb_me_fft = me * 4.0_dp * pi

  end function coulomb_me_fft

end module me_fft_direct
