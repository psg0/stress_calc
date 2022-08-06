module setting_module
  implicit none
  private

  !==================const=================================
  real(8), parameter, public :: pi=4.0_8*atan(1.0_8)
  real(8), parameter, public :: pi2=pi**2, pi3=pi**3,hpi=pi*0.50_8
  
  !==================parameter=============================
  !real(8), parameter, public :: alpha=1.0_8, beta=alpha/sqrt(3.0_8)
  real(8), parameter, public :: alpha=sqrt(3.0_8), beta=1.0_8
  real(8), parameter, public :: alpha2=alpha**2, beta2=beta**2, alpha3=alpha2*alpha, beta3=beta2*beta
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
  real(8) :: D1(Mmax,Nmax),D2(Mmax,Nmax),D3(Mmax,Nmax) !slip velocity
  real(8) :: D1hatFp(Mmax,Nmax), D2hatFp(Mmax,Nmax), D3hatFp(Mmax,Nmax) !amplitude
  real(8) :: D1hatFs(Mmax,Nmax), D2hatFs(Mmax,Nmax), D3hatFs(Mmax,Nmax)
  real(8) :: Ds1(Mmax,Nmax), Ds2(Mmax,Nmax), Ds3(Mmax,Nmax) !cumulative slip vel=> slip
    
  public triangular_mesh, make_slip_vel
contains
  subroutine make_slip_vel()
    implicit none
    integer(8) :: m,i,j,k
    real(8) :: x,y,z,r

    do i=1,Nmax
       x=xg(i)-floor(Lx*0.50)
       y=yg(i)-floor(Ly*0.50)
       z=zg(i)-floor(Lz*0.50)
       r=sqrt(x*x+y*y+z*z)
       if (r<=1) then
          write(*,*) "source index",i,r
          do m=1,Mmax
             D1(m,i)=1.0_8
             D2(m,i)=0.0_8
             D3(m,i)=0.0_8
          end do
       else
          do m=1,Mmax
             D1(m,i)=0.0_8
             D2(m,i)=0.0_8
             D3(m,i)=0.0_8
          end do
       end if
    end do
    
  end subroutine make_slip_vel
  
  subroutine triangular_mesh()
    implicit none
    integer(8) :: i,j,k,ind,ind1,ind2
    real(8) :: x,y,z,xpre,xnex,ypre

    do k=1,Nz
       ind1=Nx*Ny*(k-1)
       z=Ds*(k-1)
       do j=1,Ny
          ind2=Nx*(j-1)
          y=Ds*sqrt(3.0_8)*0.50_8*j
          ypre=Ds*sqrt(3.0_8)*0.50_8*(j-1)
          do i=1,Nx,2 !odd
             ind=ind1+ind2+i
             x=Ds*0.50_8*i
             xpre=Ds*0.50_8*(i-1)
             xnex=xpre+Ds
             !A
             xa(ind)=x
             ya(ind)=y
             za(ind)=z
             !B
             xb(ind)=xpre
             yb(ind)=ypre
             zb(ind)=z
             !C
             xc(ind)=xnex
             yc(ind)=ypre
             zc(ind)=z
             !G
             xg(ind)=x !(x+xpre+xnex)/3
             yg(ind)=(y+ypre+ypre)/3
             zg(ind)=z
          end do
          do i=2,Nx,2 !even
             ind=ind1+ind2+i
             x=Ds*0.50_8*i
             xpre=Ds*0.50_8*(i-1)
             xnex=xpre+Ds
             !A
             xa(ind)=x
             ya(ind)=ypre
             za(ind)=z
             !B
             xb(ind)=xnex
             yb(ind)=y
             zb(ind)=z
             !C
             xc(ind)=xpre
             yc(ind)=y
             zc(ind)=z
             !G
             xg(ind)=x !(x+xnex+xpre)/3
             yg(ind)=(ypre+y+y)/3
             zg(ind)=z
          end do
       end do
    end do

  end subroutine triangular_mesh
  
end module setting_module
