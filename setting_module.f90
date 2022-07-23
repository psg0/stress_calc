module setting_module
  implicit none
  private

  !==================const=================================
  real(8), parameter, public :: pi=4.0_8*atan(1.0_8)
  real(8), parameter, public :: pi2=pi**2, pi3=pi**3
  
  !==================parameter=============================
  real(8), parameter, public  :: alpha=1.0_8, beta=alpha/sqrt(3.0_8)
  integer(8), parameter, public :: Lx=40, Ly=5, Lz=2
  real(8), parameter, public :: Dx=1.0_8, Ds=Dx, Dt=0.25_8
  integer(8), parameter, public :: Nx=int((2*Lx)/Ds)-1, Ny=int(Ly/(sqrt(3.0_8)/2*Ds)), Nz=int(Lz/Ds)+1
  integer(8), parameter, public :: Nmax=Nx*Ny*Nz
  integer(8), parameter, public :: Mmax=80, Nqmax=2, Mst=1 !Mmax: Maximum time index, Nqmax: Number of quantization, Mst: index of the start of time
  
  !=========================================================
  real(8), public :: xa(Nmax), ya(Nmax), za(Nmax)
  real(8), public :: xb(Nmax), yb(Nmax), zb(Nmax)
  real(8), public :: xc(Nmax), yc(Nmax), zc(Nmax)
  real(8), public :: xg(Nmax), yg(Nmax), zg(Nmax) !Coordinates of the center of gravity
  real(8) :: sigma11(Mmax,Nmax), sigma12(Mmax,Nmax), sigma13(Mmax,Nmax)
  real(8) :: sigma21(Mmax,Nmax), sigma22(Mmax,Nmax), sigma23(Mmax,Nmax)
  real(8) :: sigma31(Mmax,Nmax), sigma32(Mmax,Nmax), sigma33(Mmax,Nmax)
  real(8) :: D1(Mmax,Nmax),D2(Mmax,Nmax),D3(Mmax,Nmax) !slip
  real(8) :: D1hatFp(Mmax,Nmax), D2hatFp(Mmax,Nmax), D3hatFp(Mmax,Nmax) !amplitude
  real(8) :: D1hatFs(Mmax,Nmax), D2hatFs(Mmax,Nmax), D3hatFs(Mmax,Nmax)
  real(8) :: Ds1(Mmax,Nmax), Ds2(Mmax,Nmax), Ds3(Mmax,Nmax) !cumulative slip
    
  public triangular_mesh!, make_slip_vel
contains

  subroutine triangular_mesh()
    implicit none
    integer(8) :: i,j,k,ind,ind1,ind2
    real(8) :: x,y,z
    !real(8) ::
    do k=1,Nz
       ind1=Nx*Ny*(k-1)
       z=Ds*(k-1)
       do j=1,Ny
          ind2=Nx*(j-1)
          y=sqrt(3.0_8)*0.50_8*j
          do i=1,Nx,2
             ind=ind1+ind2+i
             x=Ds*0.50_8*i
             xa(ind)=x
             ya(ind)=y
             za(ind)=z
             ga(ind)=(x+y+z)/3
          end do
       end do
    end do
       

  end subroutine triangular_mesh
end module setting_module
