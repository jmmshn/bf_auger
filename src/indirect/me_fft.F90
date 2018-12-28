module me_fft
  
  implicit none
  save
  
  include 'fftw3.f'
  integer(8) :: plan_f, plan_b
  
  integer , parameter :: dp = selected_real_kind(15,300)
  real(dp), parameter :: pi = 3.14159265358979323846264338327950288_dp
  real(dp), parameter :: cmm3_2_au = 1.48184743477E-25_dp
  real(dp), parameter :: K2Ha = 3.16682968067E-6_dp
  
  complex(dp), allocatable :: fftw_in_b(:,:,:), fftw_out_b(:,:,:)
  complex(dp), allocatable :: fftw_in_f(:,:,:), fftw_out_f(:,:,:)
  complex(dp), allocatable :: I13(:,:,:), I24(:,:,:)
  
  integer, allocatable :: Gmin_index(:,:,:,:)
  integer  :: nx, ny, nz
  
  real(dp) :: b(3,3)
  
  real(dp) :: kT

  
  !-- Dielectric function parameters:
  real(dp), parameter :: alpha    = 1.563_dp
  real(dp) :: q2_TF
  real(dp) :: omega2_p
  real(dp) :: lambda2
  
  private
  public :: initialize_me_fft, coulomb_me_fft
  
contains
  
  subroutine initialize_me_fft(nx1,ny1,nz1,a0,latVec11,latVec12,latVec13,latVec21,latVec22,latVec23,latVec31,latVec32,latVec33)
    implicit none
    integer, intent(in) :: nx1, ny1, nz1
    real(dp) :: a0
    real(dp) :: latVec11
    real(dp) :: latVec12
    real(dp) :: latVec13
    real(dp) :: latVec21
    real(dp) :: latVec22
    real(dp) :: latVec23
    real(dp) :: latVec31
    real(dp) :: latVec32
    real(dp) :: latVec33
    integer :: i, j, k, i1, i2, i3
    real(dp) :: Gtemp(3), Gnorm2, Gmin_norm2

    b(:,1) = (/ latVec11 , latVec12, latVec13 /)
    b(:,2) = (/ latVec21 , latVec22, latVec23 /)
    b(:,3) = (/ latVec31 , latVec32, latVec33 /)
    
    b = b*2.0_dp*pi/a0
    
    nx = nx1
    ny = ny1
    nz = nz1
    
    allocate( fftw_in_b(nx,ny,nz), fftw_out_b(nx,ny,nz) )
    allocate( fftw_in_f(nx,ny,nz), fftw_out_f(nx,ny,nz) )
    CALL dfftw_plan_dft_3d(plan_b,nx,ny,nz,fftw_in_b,fftw_out_b,FFTW_BACKWARD,FFTW_MEASURE)
    CALL dfftw_plan_dft_3d(plan_f,nx,ny,nz,fftw_in_f,fftw_out_f,FFTW_FORWARD, FFTW_MEASURE)

    allocate( I13(nx,ny,nz), I24(nx,ny,nz) )

    allocate(Gmin_index(3,nx,ny,nz))
    do k = 1, nz ; do j = 1, ny ; do i = 1, nx
       !-- Map G vector to Wigner-Seitz supercell
       Gtemp(:) = (i-1)*b(:,1) + (j-1)*b(:,2) + (k-1)*b(:,3)
       Gnorm2 = dot_product(Gtemp,Gtemp)
       Gmin_index(:,i,j,k) = (/ i-1, j-1, k-1 /)
       Gmin_norm2 = Gnorm2
       do i1 = -1, 1 ; do i2 = -1, 1 ; do i3 = -1, 1
          Gtemp(:) = (i-1+i1*nx)*b(:,1) + (j-1+i2*ny)*b(:,2) + (k-1+i3*nz)*b(:,3)
          Gnorm2 = dot_product(Gtemp,Gtemp)
          if ( Gnorm2 .lt. Gmin_norm2 ) then
             Gmin_index(:,i,j,k) = (/ i-1+i1*nx, j-1+i2*ny, k-1+i3*nz /)
             Gmin_norm2 = Gnorm2
          end if
       end do ; end do ; end do
    end do ; end do ; end do
    
  end subroutine initialize_me_fft


  complex(dp) function coulomb_me_fft (wfn1, wfn2, wfn3, wfn4, mtrans,T_Kelvin,n_val,n_free,epsilon_infty)
    implicit none
    complex(dp) :: wfn1(:), wfn2(:), wfn3(:), wfn4(:)
    real(dp)    :: mtrans(:) !-- Coulomb momentum transfer
    real(dp)    :: T_Kelvin
    real(dp)    :: n_val
    real(dp)    :: n_free
    real(dp)    :: epsilon_infty
    integer  :: index, i, j, k

    complex(dp) :: me
    real(dp)    :: Gtemp(3), Gnorm2, epsilon
    
    kT = T_Kelvin * K2Ha
    q2_TF    = 4.0_dp * (3.0_dp/pi * n_val)**(1.0_dp/3.0_dp)
    omega2_p = 4.0_dp * pi * n_val
    lambda2  = 4.0_dp * pi * n_free / kT
    !-- FFT
    index = 0
    do k = 1, nz; do j = 1, ny; do i = 1, nx
       index = index + 1
       I13(i,j,k) = conjg(wfn1(index)) * wfn3(index) / (nx*ny*nz)
       I24(i,j,k) = conjg(wfn2(index)) * wfn4(index) / (nx*ny*nz)
    end do;       end do;       end do

    fftw_in_b = I13; call dfftw_execute(plan_b); I13 = fftw_out_b
    fftw_in_f = I24; call dfftw_execute(plan_f); I24 = fftw_out_f

    !-- Calculate ME
    me = CMPLX(0.0_dp,0.0_dp)
    do k = 1, nz; do j = 1, ny; do i = 1, nx
       Gtemp  = MATMUL( b, Gmin_index(:,i,j,k) + mtrans(:) )
       Gnorm2 = dot_product(Gtemp,Gtemp)
       epsilon = 1.0_dp + 1.0_dp/( 1.0_dp/(epsilon_infty-1.0_dp) + alpha*Gnorm2/q2_TF + Gnorm2**2/(4.0_dp*omega2_p) ) 
       me = me + (1.0_dp/epsilon) * I13(i,j,k) * I24(i,j,k) / ( Gnorm2 + lambda2 )
    end do;       end do;       end do
    coulomb_me_fft = me * 4.0_dp * pi

  end function coulomb_me_fft
  
end module me_fft
