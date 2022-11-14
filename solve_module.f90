module solve_module
  use setting_module
  use kernel_module
  use func_module
  use test_module
  implicit none
  private
  public exact_solve,fdph_solve,error
  ! public exact_new, fdph_new,error
contains
  ! subroutine exact_solve(s11,s12,s13,s21,s22,s23,s31,s32,s33) !old
  !   implicit none
  !   integer :: i,j,m,dm
  !   real(8), dimension(Mmax,Nmax), intent(inout) :: s11,s12,s13,s21,s22,s23,s31,s32,s33
  !   real(8) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t
  !   real(8) :: dd1,dd2,dd3,keisu
  !   open(12,file="calc_process_exact.dat",status="replace")
  !   !===initialization===
  !   call sigma_ini(s11,s12,s13,s21,s22,s23,s31,s32,s33)
  !   !===time marching=== honest convolution
  !   do i=1,Nmax !receiver   
  !      if (mod(i,50)==0) then
  !         write(*,*) "now:", i
  !         write(12,*) "now:", i
  !      end if
  !      xr1=xg(i)
  !      xr2=yg(i)
  !      xr3=zg(i)
  !      do j=1,Nmax !source
  !         xsa1=xa(j)
  !         xsa2=ya(j)
  !         xsb1=xb(j)
  !         xsb2=yb(j)
  !         xsc1=xc(j)
  !         xsc2=yc(j)
  !         xs3=za(j) !=zb(j)=zc(j)
  !         call Kdyns_all(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
  !         do m=1,Mmax
  !            do dm=1,m-1
  !               dd1=D1(m-dm,j); dd2=D2(m-dm,j); dd3=D3(m-dm,j)
  !               t=Dt*dm
  !               s11(m,i)=s11(m,i)+Kdynf111(t,Dt)*dd1+Kdynf112(t,Dt)*dd2+Kdynf113(t,Dt)*dd3
  !               s12(m,i)=s12(m,i)+Kdynf121(t,Dt)*dd1+Kdynf122(t,Dt)*dd2+Kdynf123(t,Dt)*dd3
  !               s13(m,i)=s13(m,i)+Kdynf131(t,Dt)*dd1+Kdynf132(t,Dt)*dd2+Kdynf133(t,Dt)*dd3
  !               s21(m,i)=s21(m,i)+Kdynf211(t,Dt)*dd1+Kdynf212(t,Dt)*dd2+Kdynf213(t,Dt)*dd3
  !               s22(m,i)=s22(m,i)+Kdynf221(t,Dt)*dd1+Kdynf222(t,Dt)*dd2+Kdynf223(t,Dt)*dd3
  !               s23(m,i)=s23(m,i)+Kdynf231(t,Dt)*dd1+Kdynf232(t,Dt)*dd2+Kdynf233(t,Dt)*dd3
  !               s31(m,i)=s31(m,i)+Kdynf311(t,Dt)*dd1+Kdynf312(t,Dt)*dd2+Kdynf313(t,Dt)*dd3
  !               s32(m,i)=s32(m,i)+Kdynf321(t,Dt)*dd1+Kdynf322(t,Dt)*dd2+Kdynf323(t,Dt)*dd3
  !               s33(m,i)=s33(m,i)+Kdynf331(t,Dt)*dd1+Kdynf332(t,Dt)*dd2+Kdynf333(t,Dt)*dd3
  !            end do
  !         end do
  !      end do
  !   end do
  !   keisu=-mu*0.50_8/beta
  !   do i=1,Nmax
  !      do m=1,Mmax
  !         s11(m,i)=keisu*s11(m,i)
  !         s12(m,i)=keisu*s12(m,i)
  !         s13(m,i)=keisu*s13(m,i)
  !         s21(m,i)=keisu*s21(m,i)
  !         s22(m,i)=keisu*s22(m,i)
  !         s23(m,i)=keisu*s23(m,i)
  !         s31(m,i)=keisu*s31(m,i)
  !         s32(m,i)=keisu*s32(m,i)
  !         s33(m,i)=keisu*s33(m,i)
  !      end do
  !   end do
  !   close(12)
  ! end subroutine exact_solve !old

  subroutine exact_solve(s11,s12,s13,s21,s22,s23,s31,s32,s33) !new
    implicit none
    integer :: i,j,m,dm
    real(8), dimension(Mmax,Nmax), intent(inout) :: s11,s12,s13,s21,s22,s23,s31,s32,s33
    real(8) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t
    real(8) :: dd1,dd2,dd3,keisu
    real(8), dimension(Mmax-1,Nmax,Nmax) :: K111,K112,K113,K121,K122,K123,K131,K132,K133
    real(8), dimension(Mmax-1,Nmax,Nmax) :: K211,K212,K213,K221,K222,K223,K231,K232,K233
    real(8), dimension(Mmax-1,Nmax,Nmax) :: K311,K312,K313,K321,K322,K323,K331,K332,K333
    open(12,file="calc_process_exact.dat",status="replace")
    !===initialization===
    call sigma_ini(s11,s12,s13,s21,s22,s23,s31,s32,s33)
    !===time marching=== honest convolution
    do i=1,Nmax !receiver   
       xr1=xg(i)
       xr2=yg(i)
       xr3=zg(i)
       do j=1,Nmax !source
          xsa1=xa(j)
          xsa2=ya(j)
          xsb1=xb(j)
          xsb2=yb(j)
          xsc1=xc(j)
          xsc2=yc(j)
          xs3=za(j) !=zb(j)=zc(j)
          call Kdyns_all(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
          do dm=1,Mmax-1
             t=Dt*dm
             K111(dm,i,j)=Kdynf111(t,Dt)
             K112(dm,i,j)=Kdynf112(t,Dt)
             K113(dm,i,j)=Kdynf113(t,Dt)
             K121(dm,i,j)=Kdynf121(t,Dt)
             K122(dm,i,j)=Kdynf122(t,Dt)
             K123(dm,i,j)=Kdynf123(t,Dt)
             K131(dm,i,j)=Kdynf131(t,Dt)
             K132(dm,i,j)=Kdynf132(t,Dt)
             K133(dm,i,j)=Kdynf133(t,Dt)
             K211(dm,i,j)=Kdynf211(t,Dt)
             K212(dm,i,j)=Kdynf212(t,Dt)
             K213(dm,i,j)=Kdynf213(t,Dt)
             K221(dm,i,j)=Kdynf221(t,Dt)
             K222(dm,i,j)=Kdynf222(t,Dt)
             K223(dm,i,j)=Kdynf223(t,Dt)
             K231(dm,i,j)=Kdynf231(t,Dt)
             K232(dm,i,j)=Kdynf232(t,Dt)
             K233(dm,i,j)=Kdynf233(t,Dt)
             K311(dm,i,j)=Kdynf311(t,Dt)
             K312(dm,i,j)=Kdynf312(t,Dt)
             K313(dm,i,j)=Kdynf313(t,Dt)
             K321(dm,i,j)=Kdynf321(t,Dt)
             K322(dm,i,j)=Kdynf322(t,Dt)
             K323(dm,i,j)=Kdynf323(t,Dt)
             K331(dm,i,j)=Kdynf331(t,Dt)
             K332(dm,i,j)=Kdynf332(t,Dt)
             K333(dm,i,j)=Kdynf333(t,Dt)
          end do
       end do
    end do
    !
    do i=1,Nmax !receiver
       if (mod(i,50)==0) then
          write(*,*) "now:", i
          write(12,*) "now:", i
       end if
       !xr1=xg(i)
       !xr2=yg(i)
       !xr3=zg(i)
       do j=1,Nmax !source
          !xsa1=xa(j)
          !xsa2=ya(j)
          !xsb1=xb(j)
          !xsb2=yb(j)
          !xsc1=xc(j)
          !xsc2=yc(j)
          !xs3=za(j) !=zb(j)=zc(j)
          do m=1,Mmax
             do dm=1,m-1
                dd1=D1(m-dm,j); dd2=D2(m-dm,j); dd3=D3(m-dm,j)
                !t=Dt*dm
                s11(m,i)=s11(m,i)+K111(dm,i,j)*dd1+K112(dm,i,j)*dd2+K113(dm,i,j)*dd3
                s12(m,i)=s12(m,i)+K121(dm,i,j)*dd1+K122(dm,i,j)*dd2+K123(dm,i,j)*dd3
                s13(m,i)=s13(m,i)+K131(dm,i,j)*dd1+K132(dm,i,j)*dd2+K133(dm,i,j)*dd3
                s21(m,i)=s21(m,i)+K211(dm,i,j)*dd1+K212(dm,i,j)*dd2+K213(dm,i,j)*dd3
                s22(m,i)=s22(m,i)+K221(dm,i,j)*dd1+K222(dm,i,j)*dd2+K223(dm,i,j)*dd3
                s23(m,i)=s23(m,i)+K231(dm,i,j)*dd1+K232(dm,i,j)*dd2+K233(dm,i,j)*dd3
                s31(m,i)=s31(m,i)+K311(dm,i,j)*dd1+K312(dm,i,j)*dd2+K313(dm,i,j)*dd3
                s32(m,i)=s32(m,i)+K321(dm,i,j)*dd1+K322(dm,i,j)*dd2+K323(dm,i,j)*dd3
                s33(m,i)=s33(m,i)+K331(dm,i,j)*dd1+K332(dm,i,j)*dd2+K333(dm,i,j)*dd3
             end do
          end do
       end do
    end do
    keisu=-mu*0.50_8/beta
    do i=1,Nmax
       do m=1,Mmax
          s11(m,i)=keisu*s11(m,i)
          s12(m,i)=keisu*s12(m,i)
          s13(m,i)=keisu*s13(m,i)
          s21(m,i)=keisu*s21(m,i)
          s22(m,i)=keisu*s22(m,i)
          s23(m,i)=keisu*s23(m,i)
          s31(m,i)=keisu*s31(m,i)
          s32(m,i)=keisu*s32(m,i)
          s33(m,i)=keisu*s33(m,i)
       end do
    end do
    close(12)
  end subroutine exact_solve !new

  subroutine fdph_solve(s11,s12,s13,s21,s22,s23,s31,s32,s33) !o(N^2M) new
    implicit none
    integer :: i,j,bl,m,dm,k
    real(8), dimension(Mmax,Nmax), intent(inout) :: s11,s12,s13,s21,s22,s23,s31,s32,s33
    real(8) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,xsg1,xsg2,xsg3,r,t,t1
    real(8) :: krnFp_lis(27),krnFs_lis(27),htmFp_lis(27),htmFs_lis(27)
    real(8) :: krnInsmp_lis(0:Nqmax,27,Nmax,Nmax)
    real(8) :: ks(27)
    real(8) :: dd1,dd2,dd3,keisu
    integer :: istar_lis(Nmax),istar, block_lis(Nmax), jstar !istar_ind
    integer :: ind,nq,mip,mis,mii,mi,mj,msmin,mpmax_now,msmax_now,mm
    integer :: inds
    integer :: msmin_lis(Nmax,Nmax),mip_lis(Nmax,nblock),mis_lis(Nmax,nblock)
    integer :: cal_was(nblock,Nmax), cal_now
    real(8) :: xr1star,xr2star,xr3star,xry,yry,zry,rry,x,y,z,dri,ri,rj,rstar,rtrue
    real(8) :: val
    real(8), dimension(Nmax,Nmax) :: dis
    real(8), dimension(Mmax-1,Nmax,Nmax) :: K111,K112,K113,K121,K122,K123,K131,K132,K133
    real(8), dimension(Mmax-1,Nmax,Nmax) :: K211,K212,K213,K221,K222,K223,K231,K232,K233
    real(8), dimension(Mmax-1,Nmax,Nmax) :: K311,K312,K313,K321,K322,K323,K331,K332,K333
    real(8), dimension(27,Nmax,Nmax) :: Khatp,Khats
    real(8), dimension(0:mpmax,27,nblock,Nmax):: hp
    real(8), dimension(0:msmax-1,27,nblock,Nmax) :: hs
    real(8), dimension(27,Nmax,Nmax) :: ks_lis
    
    !===initialization===
    Ds1(:,:)=0.0_8; Ds2(:,:)=0.0_8; Ds3(:,:)=0.0_8
    call sigma_ini(s11,s12,s13,s21,s22,s23,s31,s32,s33)
    call Dhat_ini()
    call calc_istar(nxb,nyb,nzb,istar_lis,block_lis,block_center) !determine representative points of observation
    krnInsmp_lis(:,:,:,:)=0.0_8
    hp(:,:,:,:)=0.0_8; hs(:,:,:,:)=0.0_8
    cal_was(:,:)=-1
    
    !write(*,*) block_lis
    call check_same(nxb,nyb,nzb,block_lis,block_center)

    !===calculate Kernel and store===
    do j=1,Nmax
       xsa1=xa(j); xsa2=ya(j); xsb1=xb(j); xsb2=yb(j)
       xsc1=xc(j); xsc2=yc(j); xs3=za(j)
       xsg1=xg(j); xsg2=yg(j); xsg3=zg(j)
       do i=1,Nmax
          cal_now=-1
          bl=block_lis(i)
          istar=block_center(bl)
          xr1star=xg(istar); xr2star=yg(istar); xr3star=zg(istar)
          xr1=xg(i); xr2=yg(i); xr3=zg(i)
          call Kdyns_all(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
          r=sqrt((xr1-xsg1)**2+(xr2-xsg2)**2+(xr3-xsg3)**2)
          rj=sqrt((xr1star-xsg1)**2+(xr2star-xsg2)**2+(xr3star-xsg3)**2)
          dis(i,j)=r
          !Khat (Fp)  
          Khatp(:,i,j)=0.0_8
          do dm=0,mpmax-1
             t=r/alpha+dm*Dt-Ds/alpha*0.50_8 !t can be negative, but kernel become 0 when t<0
             Khatp(1,i,j)=Khatp(1,i,j)+Kdynf111(t,Dt)
             Khatp(2,i,j)=Khatp(2,i,j)+Kdynf112(t,Dt)
             Khatp(3,i,j)=Khatp(3,i,j)+Kdynf113(t,Dt)
             Khatp(4,i,j)=Khatp(4,i,j)+Kdynf121(t,Dt)
             Khatp(5,i,j)=Khatp(5,i,j)+Kdynf122(t,Dt)
             Khatp(6,i,j)=Khatp(6,i,j)+Kdynf123(t,Dt)
             Khatp(7,i,j)=Khatp(7,i,j)+Kdynf131(t,Dt)
             Khatp(8,i,j)=Khatp(8,i,j)+Kdynf132(t,Dt)
             Khatp(9,i,j)=Khatp(9,i,j)+Kdynf133(t,Dt)
             Khatp(10,i,j)=Khatp(10,i,j)+Kdynf211(t,Dt)
             Khatp(11,i,j)=Khatp(11,i,j)+Kdynf212(t,Dt)
             Khatp(12,i,j)=Khatp(12,i,j)+Kdynf213(t,Dt)
             Khatp(13,i,j)=Khatp(13,i,j)+Kdynf221(t,Dt)
             Khatp(14,i,j)=Khatp(14,i,j)+Kdynf222(t,Dt)
             Khatp(15,i,j)=Khatp(15,i,j)+Kdynf223(t,Dt)
             Khatp(16,i,j)=Khatp(16,i,j)+Kdynf231(t,Dt)
             Khatp(17,i,j)=Khatp(17,i,j)+Kdynf232(t,Dt)
             Khatp(18,i,j)=Khatp(18,i,j)+Kdynf233(t,Dt)
             Khatp(19,i,j)=Khatp(19,i,j)+Kdynf311(t,Dt)
             Khatp(20,i,j)=Khatp(20,i,j)+Kdynf312(t,Dt)
             Khatp(21,i,j)=Khatp(21,i,j)+Kdynf313(t,Dt)
             Khatp(22,i,j)=Khatp(22,i,j)+Kdynf321(t,Dt)
             Khatp(23,i,j)=Khatp(23,i,j)+Kdynf322(t,Dt)
             Khatp(24,i,j)=Khatp(24,i,j)+Kdynf323(t,Dt)
             Khatp(25,i,j)=Khatp(25,i,j)+Kdynf331(t,Dt)
             Khatp(26,i,j)=Khatp(26,i,j)+Kdynf332(t,Dt)
             Khatp(27,i,j)=Khatp(27,i,j)+Kdynf333(t,Dt)
          end do
          !Khat (Fs)
          msmin=0
          t1=r/alpha+Ds/alpha*0.50_8
          if(t1>r/beta-Ds/beta*0.50_8) then
             msmin=floor((t1-(r/beta-Ds/beta*0.50_8))/Dt)
          end if
          msmin_lis(i,j)=msmin
          !t2=r/beta+msmin*Dt-Ds/beta*0.50_8 !t2>t1
          Khats(:,i,j)=0.0_8
          do dm=msmin,msmax-1,1
             t=r/beta+dm*Dt-Ds/beta*0.50_8
             Khats(1,i,j)=Khats(1,i,j)+Kdynf111(t,Dt)
             Khats(2,i,j)=Khats(2,i,j)+Kdynf112(t,Dt)
             Khats(3,i,j)=Khats(3,i,j)+Kdynf113(t,Dt)
             Khats(4,i,j)=Khats(4,i,j)+Kdynf121(t,Dt)
             Khats(5,i,j)=Khats(5,i,j)+Kdynf122(t,Dt)
             Khats(6,i,j)=Khats(6,i,j)+Kdynf123(t,Dt)
             Khats(7,i,j)=Khats(7,i,j)+Kdynf131(t,Dt)
             Khats(8,i,j)=Khats(8,i,j)+Kdynf132(t,Dt)
             Khats(9,i,j)=Khats(9,i,j)+Kdynf133(t,Dt)
             Khats(10,i,j)=Khats(10,i,j)+Kdynf211(t,Dt)
             Khats(11,i,j)=Khats(11,i,j)+Kdynf212(t,Dt)
             Khats(12,i,j)=Khats(12,i,j)+Kdynf213(t,Dt)
             Khats(13,i,j)=Khats(13,i,j)+Kdynf221(t,Dt)
             Khats(14,i,j)=Khats(14,i,j)+Kdynf222(t,Dt)
             Khats(15,i,j)=Khats(15,i,j)+Kdynf223(t,Dt)
             Khats(16,i,j)=Khats(16,i,j)+Kdynf231(t,Dt)
             Khats(17,i,j)=Khats(17,i,j)+Kdynf232(t,Dt)
             Khats(18,i,j)=Khats(18,i,j)+Kdynf233(t,Dt)
             Khats(19,i,j)=Khats(19,i,j)+Kdynf311(t,Dt)
             Khats(20,i,j)=Khats(20,i,j)+Kdynf312(t,Dt)
             Khats(21,i,j)=Khats(21,i,j)+Kdynf313(t,Dt)
             Khats(22,i,j)=Khats(22,i,j)+Kdynf321(t,Dt)
             Khats(23,i,j)=Khats(23,i,j)+Kdynf322(t,Dt)
             Khats(24,i,j)=Khats(24,i,j)+Kdynf323(t,Dt)
             Khats(25,i,j)=Khats(25,i,j)+Kdynf331(t,Dt)
             Khats(26,i,j)=Khats(26,i,j)+Kdynf332(t,Dt)
             Khats(27,i,j)=Khats(27,i,j)+Kdynf333(t,Dt)            
          end do
          !calculate h
          if (istar==i.and.cal_was(bl,j)==-1) then
             do dm=0,mpmax-1
                t=rj/alpha+dm*Dt-Ds/alpha*0.50_8
                hp(dm,1,bl,j)=cal_h(Kdynf111(t,Dt),Khatp(1,i,j))
                hp(dm,2,bl,j)=cal_h(Kdynf112(t,Dt),Khatp(2,i,j))
                hp(dm,3,bl,j)=cal_h(Kdynf113(t,Dt),Khatp(3,i,j))
                hp(dm,4,bl,j)=cal_h(Kdynf121(t,Dt),Khatp(4,i,j))
                hp(dm,5,bl,j)=cal_h(Kdynf122(t,Dt),Khatp(5,i,j))
                hp(dm,6,bl,j)=cal_h(Kdynf123(t,Dt),Khatp(6,i,j))
                hp(dm,7,bl,j)=cal_h(Kdynf131(t,Dt),Khatp(7,i,j))
                hp(dm,8,bl,j)=cal_h(Kdynf132(t,Dt),Khatp(8,i,j))
                hp(dm,9,bl,j)=cal_h(Kdynf133(t,Dt),Khatp(9,i,j))
                hp(dm,10,bl,j)=cal_h(Kdynf211(t,Dt),Khatp(10,i,j))
                hp(dm,11,bl,j)=cal_h(Kdynf212(t,Dt),Khatp(11,i,j))
                hp(dm,12,bl,j)=cal_h(Kdynf213(t,Dt),Khatp(12,i,j))
                hp(dm,13,bl,j)=cal_h(Kdynf221(t,Dt),Khatp(13,i,j))
                hp(dm,14,bl,j)=cal_h(Kdynf222(t,Dt),Khatp(14,i,j))
                hp(dm,15,bl,j)=cal_h(Kdynf223(t,Dt),Khatp(15,i,j))
                hp(dm,16,bl,j)=cal_h(Kdynf231(t,Dt),Khatp(16,i,j))
                hp(dm,17,bl,j)=cal_h(Kdynf232(t,Dt),Khatp(17,i,j))
                hp(dm,18,bl,j)=cal_h(Kdynf233(t,Dt),Khatp(18,i,j))
                hp(dm,19,bl,j)=cal_h(Kdynf311(t,Dt),Khatp(19,i,j))
                hp(dm,20,bl,j)=cal_h(Kdynf312(t,Dt),Khatp(20,i,j))
                hp(dm,21,bl,j)=cal_h(Kdynf313(t,Dt),Khatp(21,i,j))
                hp(dm,22,bl,j)=cal_h(Kdynf321(t,Dt),Khatp(22,i,j))
                hp(dm,23,bl,j)=cal_h(Kdynf322(t,Dt),Khatp(23,i,j))
                hp(dm,24,bl,j)=cal_h(Kdynf323(t,Dt),Khatp(24,i,j))
                hp(dm,25,bl,j)=cal_h(Kdynf331(t,Dt),Khatp(25,i,j))
                hp(dm,26,bl,j)=cal_h(Kdynf332(t,Dt),Khatp(26,i,j))
                hp(dm,27,bl,j)=cal_h(Kdynf333(t,Dt),Khatp(27,i,j))
             end do
             !msmin=0
             !t1=r/alpha+Ds/alpha*0.50_8
             !if(t1>=r/beta-Ds/beta*0.50_8) then
             !   msmin=floor((t1-(r/beta-Ds/beta*0.50_8))/Dt)
             !end if
             !t2=r/beta+msmin*Dt-Ds/beta*0.50_8 !t2>t1
             do dm=msmin,msmax-1,1
                t=rj/beta+dm*Dt-Ds/beta*0.50_8
                hs(dm,1,bl,j)=cal_h(Kdynf111(t,Dt),Khats(1,i,j))
                hs(dm,2,bl,j)=cal_h(Kdynf112(t,Dt),Khats(2,i,j))
                hs(dm,3,bl,j)=cal_h(Kdynf113(t,Dt),Khats(3,i,j))
                hs(dm,4,bl,j)=cal_h(Kdynf121(t,Dt),Khats(4,i,j))
                hs(dm,5,bl,j)=cal_h(Kdynf122(t,Dt),Khats(5,i,j))
                hs(dm,6,bl,j)=cal_h(Kdynf123(t,Dt),Khats(6,i,j))
                hs(dm,7,bl,j)=cal_h(Kdynf131(t,Dt),Khats(7,i,j))
                hs(dm,8,bl,j)=cal_h(Kdynf132(t,Dt),Khats(8,i,j))
                hs(dm,9,bl,j)=cal_h(Kdynf133(t,Dt),Khats(9,i,j))
                hs(dm,10,bl,j)=cal_h(Kdynf211(t,Dt),Khats(10,i,j))
                hs(dm,11,bl,j)=cal_h(Kdynf212(t,Dt),Khats(11,i,j))
                hs(dm,12,bl,j)=cal_h(Kdynf213(t,Dt),Khats(12,i,j))
                hs(dm,13,bl,j)=cal_h(Kdynf221(t,Dt),Khats(13,i,j))
                hs(dm,14,bl,j)=cal_h(Kdynf222(t,Dt),Khats(14,i,j))
                hs(dm,15,bl,j)=cal_h(Kdynf223(t,Dt),Khats(15,i,j))
                hs(dm,16,bl,j)=cal_h(Kdynf231(t,Dt),Khats(16,i,j))
                hs(dm,17,bl,j)=cal_h(Kdynf232(t,Dt),Khats(17,i,j))
                hs(dm,18,bl,j)=cal_h(Kdynf233(t,Dt),Khats(18,i,j))
                hs(dm,19,bl,j)=cal_h(Kdynf311(t,Dt),Khats(19,i,j))
                hs(dm,20,bl,j)=cal_h(Kdynf312(t,Dt),Khats(20,i,j))
                hs(dm,21,bl,j)=cal_h(Kdynf313(t,Dt),Khats(21,i,j))
                hs(dm,22,bl,j)=cal_h(Kdynf321(t,Dt),Khats(22,i,j))
                hs(dm,23,bl,j)=cal_h(Kdynf322(t,Dt),Khats(23,i,j))
                hs(dm,24,bl,j)=cal_h(Kdynf323(t,Dt),Khats(24,i,j))
                hs(dm,25,bl,j)=cal_h(Kdynf331(t,Dt),Khats(25,i,j))
                hs(dm,26,bl,j)=cal_h(Kdynf332(t,Dt),Khats(26,i,j))
                hs(dm,27,bl,j)=cal_h(Kdynf333(t,Dt),Khats(27,i,j))
             end do
             cal_was(bl,j)=1
          end if
          !calculate kernel of Domain I and store it
          call krnInsmp(i,j,r,krnInsmp_lis)
          !calculate kernel of Domain S and store it
          call calc_ks(i,j,r,ks_lis)
       end do
    end do

    !To calculate mi
    do bl=1,nblock
       jstar=block_center(bl)
       do i=1,Nmax
          istar=istar_lis(i)
          !call calc_mi_31(istar,i,jstar,dri,mip,mis) !Sato & Ando (31)
          !Sato & Ando (42)
          ri=sqrt((xg(i)-xg(jstar))**2+(yg(i)-yg(jstar))**2+(zg(i)-zg(jstar))**2)
          rstar=sqrt((xg(istar)-xg(jstar))**2+(yg(istar)-yg(jstar))**2+(zg(istar)-zg(jstar))**2)
          dri=ri-rstar
          mip=floor(dri/(alpha*Dt))
          mis=floor(dri/(beta*Dt))
          mip_lis(i,bl)=mip
          mis_lis(i,bl)=mis
       end do
    end do
           
    ! !===time marching===
    do m=Mst+1,Mmax
       if (mod(m,10)==0) then
          write(*,*) "now: m=",m
       end if
       !calculate Ds and Dhat
       do j=1,Nmax
          Ds1(m,j)=Ds1(m-1,j)+D1(m,j)!*Dt
          Ds2(m,j)=Ds2(m-1,j)+D2(m,j)!*Dt
          Ds3(m,j)=Ds3(m-1,j)+D3(m,j)!*Dt
          !xsa1=xa(j); xsa2=ya(j); xsb1=xb(j); xsb2=yb(j)
          !xsc1=xc(j); xsc2=yc(j); xs3=za(j)
          !xsg1=xg(j); xsg2=yg(j); xsg3=zg(j)
          do bl=1,nblock
             istar=block_center(bl)
             rj=dis(istar,j)
             !xr1star=xg(istar); xr2star=yg(istar); xr3star=zg(istar)
             !xr1=xg(i); xr2=yg(i); xr3=zg(i)
             !call Kdyns_all(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
             !rtrue=sqrt((xr1-xsg1)**2+(xr2-xsg2)**2+(xr3-xsg3)**2)
             !call krnFp_all(i,j,rtrue,krnFp_lis)
             !call krnFs_all(i,j,rtrue,krnFs_lis
             do dm=0,mpmax-1
                ind=m-dm !Sato & Ando (52)
                if (ind<1 .or. ind>=m) then
                   cycle
                end if
                dd1=D1(ind,j); dd2=D2(ind,j); dd3=D3(ind,j)
                D111hatFp(m,bl,j)=D111hatFp(m,bl,j)+dd1*hp(dm,1,bl,j)
                D112hatFp(m,bl,j)=D112hatFp(m,bl,j)+dd2*hp(dm,2,bl,j)
                D113hatFp(m,bl,j)=D113hatFp(m,bl,j)+dd3*hp(dm,3,bl,j)
                D121hatFp(m,bl,j)=D121hatFp(m,bl,j)+dd1*hp(dm,4,bl,j)
                D122hatFp(m,bl,j)=D122hatFp(m,bl,j)+dd2*hp(dm,5,bl,j)
                D123hatFp(m,bl,j)=D123hatFp(m,bl,j)+dd3*hp(dm,6,bl,j)
                D131hatFp(m,bl,j)=D131hatFp(m,bl,j)+dd1*hp(dm,7,bl,j)
                D132hatFp(m,bl,j)=D132hatFp(m,bl,j)+dd2*hp(dm,8,bl,j)
                D133hatFp(m,bl,j)=D133hatFp(m,bl,j)+dd3*hp(dm,9,bl,j)
                D211hatFp(m,bl,j)=D211hatFp(m,bl,j)+dd1*hp(dm,10,bl,j)
                D212hatFp(m,bl,j)=D212hatFp(m,bl,j)+dd2*hp(dm,11,bl,j)
                D213hatFp(m,bl,j)=D213hatFp(m,bl,j)+dd3*hp(dm,12,bl,j)
                D221hatFp(m,bl,j)=D221hatFp(m,bl,j)+dd1*hp(dm,13,bl,j)
                D222hatFp(m,bl,j)=D222hatFp(m,bl,j)+dd2*hp(dm,14,bl,j)
                D223hatFp(m,bl,j)=D223hatFp(m,bl,j)+dd3*hp(dm,15,bl,j)
                D231hatFp(m,bl,j)=D231hatFp(m,bl,j)+dd1*hp(dm,16,bl,j)
                D232hatFp(m,bl,j)=D232hatFp(m,bl,j)+dd2*hp(dm,17,bl,j)
                D233hatFp(m,bl,j)=D233hatFp(m,bl,j)+dd3*hp(dm,18,bl,j)
                D311hatFp(m,bl,j)=D311hatFp(m,bl,j)+dd1*hp(dm,19,bl,j)
                D312hatFp(m,bl,j)=D312hatFp(m,bl,j)+dd2*hp(dm,20,bl,j)
                D313hatFp(m,bl,j)=D313hatFp(m,bl,j)+dd3*hp(dm,21,bl,j)
                D321hatFp(m,bl,j)=D321hatFp(m,bl,j)+dd1*hp(dm,22,bl,j)
                D322hatFp(m,bl,j)=D322hatFp(m,bl,j)+dd2*hp(dm,23,bl,j)
                D323hatFp(m,bl,j)=D323hatFp(m,bl,j)+dd3*hp(dm,24,bl,j)
                D331hatFp(m,bl,j)=D331hatFp(m,bl,j)+dd1*hp(dm,25,bl,j)
                D332hatFp(m,bl,j)=D332hatFp(m,bl,j)+dd2*hp(dm,26,bl,j)
                D333hatFp(m,bl,j)=D333hatFp(m,bl,j)+dd3*hp(dm,27,bl,j)
             end do

             msmin=0
             t1=rj/alpha+Ds/alpha*0.50_8
             if(t1>rj/beta-Ds/beta*0.50_8) then
                msmin=floor((t1-(rj/beta-Ds/beta*0.50_8))/Dt)
             end if
             do dm=msmin,msmax-1
                ind=m-dm
                if (ind<1 .or. ind>=m) then
                   cycle
                end if
                !call htmFs(istar,j,dm,krnFs_lis,htmFs_lis)
                dd1=D1(ind,j); dd2=D2(ind,j); dd3=D3(ind,j)
                D111hatFs(m,bl,j)=D111hatFs(m,bl,j)+dd1*hs(dm,1,bl,j)
                D112hatFs(m,bl,j)=D112hatFs(m,bl,j)+dd2*hs(dm,2,bl,j)
                D113hatFs(m,bl,j)=D113hatFs(m,bl,j)+dd3*hs(dm,3,bl,j)
                D121hatFs(m,bl,j)=D121hatFs(m,bl,j)+dd1*hs(dm,4,bl,j)
                D122hatFs(m,bl,j)=D122hatFs(m,bl,j)+dd2*hs(dm,5,bl,j)
                D123hatFs(m,bl,j)=D123hatFs(m,bl,j)+dd3*hs(dm,6,bl,j)
                D131hatFs(m,bl,j)=D131hatFs(m,bl,j)+dd1*hs(dm,7,bl,j)
                D132hatFs(m,bl,j)=D132hatFs(m,bl,j)+dd2*hs(dm,8,bl,j)
                D133hatFs(m,bl,j)=D133hatFs(m,bl,j)+dd3*hs(dm,9,bl,j)
                D211hatFs(m,bl,j)=D211hatFs(m,bl,j)+dd1*hs(dm,10,bl,j)
                D212hatFs(m,bl,j)=D212hatFs(m,bl,j)+dd2*hs(dm,11,bl,j)
                D213hatFs(m,bl,j)=D213hatFs(m,bl,j)+dd3*hs(dm,12,bl,j)
                D221hatFs(m,bl,j)=D221hatFs(m,bl,j)+dd1*hs(dm,13,bl,j)
                D222hatFs(m,bl,j)=D222hatFs(m,bl,j)+dd2*hs(dm,14,bl,j)
                D223hatFs(m,bl,j)=D223hatFs(m,bl,j)+dd3*hs(dm,15,bl,j)
                D231hatFs(m,bl,j)=D231hatFs(m,bl,j)+dd1*hs(dm,16,bl,j)
                D232hatFs(m,bl,j)=D232hatFs(m,bl,j)+dd2*hs(dm,17,bl,j)
                D233hatFs(m,bl,j)=D233hatFs(m,bl,j)+dd3*hs(dm,18,bl,j)
                D311hatFs(m,bl,j)=D311hatFs(m,bl,j)+dd1*hs(dm,19,bl,j)
                D312hatFs(m,bl,j)=D312hatFs(m,bl,j)+dd2*hs(dm,20,bl,j)
                D313hatFs(m,bl,j)=D313hatFs(m,bl,j)+dd3*hs(dm,21,bl,j)
                D321hatFs(m,bl,j)=D321hatFs(m,bl,j)+dd1*hs(dm,22,bl,j)
                D322hatFs(m,bl,j)=D322hatFs(m,bl,j)+dd2*hs(dm,23,bl,j)
                D323hatFs(m,bl,j)=D323hatFs(m,bl,j)+dd3*hs(dm,24,bl,j)
                D331hatFs(m,bl,j)=D331hatFs(m,bl,j)+dd1*hs(dm,25,bl,j)
                D332hatFs(m,bl,j)=D332hatFs(m,bl,j)+dd2*hs(dm,26,bl,j)
                D333hatFs(m,bl,j)=D333hatFs(m,bl,j)+dd3*hs(dm,27,bl,j)
             end do
          end do
       end do
       !convolution start !!!
       do i=1,Nmax !receiver
          xr1=xg(i)
          xr2=yg(i)
          xr3=zg(i)
          istar=istar_lis(i)
          bl=block_lis(i)
          xr1star=xg(istar)
          xr2star=yg(istar)
          xr3star=zg(istar)
          do j=1,Nmax !source
             xsa1=xa(j)
             xsa2=ya(j)
             xsb1=xb(j)
             xsb2=yb(j)
             xsc1=xc(j)
             xsc2=yc(j)
             xs3=za(j) !=zb(j)=zc(j)
             xsg1=xg(j); xsg2=yg(j); xsg3=zg(j)
             r=dis(i,j)
             !computation of mi
             rj=dis(istar,j) !sqrt((xr1star-xsg1)**2+(xr2star-xsg2)**2+(xr3star-xsg3)**2)
             mj=max(ceiling((rj-Ds*0.50_8)/(alpha*Dt)),0)
             jstar=istar_lis(j)
             mip=mip_lis(i,block_lis(jstar))
             mis=mis_lis(i,block_lis(jstar))
     
             !Amplitude 
             !rtrue=sqrt((xr1-xsg1)**2+(xr2-xsg2)**2+(xr3-xsg3)**2)
             !call Kdyns_all(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
             !call krnFp_all(i,j,rtrue,krnFp_lis)
             !call krnFs_all(i,j,rtrue,krnFs_lis)
             !domain Fp
             !mj=max(ceiling(rj/(alpha*Dt)-Ds/(alpha*Dt)*0.5),1)
             ind=m+mip
             mm=m-mj
             !if (mj<=m-1.and.1<=ind.and.ind<=Mmax) then
             !write(*,*) m,mip,ind,Mmax
             if (1<=ind.and.ind<=Mmax.and.mm>=1) then
                s11(ind,i)=s11(ind,i)+D111hatFp(mm,bl,j)*Khatp(1,i,j)&
                     & +D112hatFp(mm,bl,j)*Khatp(2,i,j)+D113hatFp(mm,bl,j)*Khatp(3,i,j)
                s12(ind,i)=s12(ind,i)+D121hatFp(mm,bl,j)*Khatp(4,i,j)&
                     & +D122hatFp(mm,bl,j)*Khatp(5,i,j)+D123hatFp(mm,bl,j)*Khatp(6,i,j)
                s13(ind,i)=s13(ind,i)+D131hatFp(mm,bl,j)*Khatp(7,i,j)&
                     & +D132hatFp(mm,bl,j)*Khatp(8,i,j)+D133hatFp(mm,bl,j)*Khatp(9,i,j)
                s21(ind,i)=s21(ind,i)+D211hatFp(mm,bl,j)*Khatp(10,i,j)&
                     & +D212hatFp(mm,bl,j)*Khatp(11,i,j)+D213hatFp(mm,bl,j)*Khatp(12,i,j)
                s22(ind,i)=s22(ind,i)+D221hatFp(mm,bl,j)*Khatp(13,i,j)&
                     & +D222hatFp(mm,bl,j)*Khatp(14,i,j)+D223hatFp(mm,bl,j)*Khatp(15,i,j)
                s23(ind,i)=s23(ind,i)+D231hatFp(mm,bl,j)*Khatp(16,i,j)&
                     & +D232hatFp(mm,bl,j)*Khatp(17,i,j)+D233hatFp(mm,bl,j)*Khatp(18,i,j)
                s31(ind,i)=s31(ind,i)+D311hatFp(mm,bl,j)*Khatp(19,i,j)&
                     & +D312hatFp(mm,bl,j)*Khatp(20,i,j)+D313hatFp(mm,bl,j)*Khatp(21,i,j)
                s32(ind,i)=s32(ind,i)+D321hatFp(mm,bl,j)*Khatp(22,i,j)&
                     & +D322hatFp(mm,bl,j)*Khatp(23,i,j)+D323hatFp(mm,bl,j)*Khatp(24,i,j)
                s33(ind,i)=s33(ind,i)+D331hatFp(mm,bl,j)*Khatp(25,i,j)&
                     & +D332hatFp(mm,bl,j)*Khatp(26,i,j)+D333hatFp(mm,bl,j)*Khatp(27,i,j)
             end if
             !domain Fs
             !mj=max(ceiling(rj/(beta*Dt)-Ds/(beta*Dt)*0.5),1)
             ind=m+mis
             !if (mj<=m-1.and.1<=ind.and.ind<=Mmax) then
             if (1<=ind.and.ind<=Mmax.and.mm>=1) then 
                s11(ind,i)=s11(ind,i)+D111hatFs(mm,bl,j)*Khats(1,i,j)&
                     & +D112hatFs(mm,bl,j)*Khats(2,i,j)+D113hatFs(mm,bl,j)*Khats(3,i,j)
                s12(ind,i)=s12(ind,i)+D121hatFs(mm,bl,j)*Khats(4,i,j)&
                     & +D122hatFs(mm,bl,j)*Khats(5,i,j)+D123hatFs(mm,bl,j)*Khats(6,i,j)
                s13(ind,i)=s13(ind,i)+D131hatFs(mm,bl,j)*Khats(7,i,j)&
                     & +D132hatFs(mm,bl,j)*Khats(8,i,j)+D133hatFs(mm,bl,j)*Khats(9,i,j)
                s21(ind,i)=s21(ind,i)+D211hatFs(mm,bl,j)*Khats(10,i,j)&
                     & +D212hatFs(mm,bl,j)*Khats(11,i,j)+D213hatFs(mm,bl,j)*Khats(12,i,j)
                s22(ind,i)=s22(ind,i)+D221hatFs(mm,bl,j)*Khats(13,i,j)&
                     & +D222hatFs(mm,bl,j)*Khats(14,i,j)+D223hatFs(mm,bl,j)*Khats(15,i,j)
                s23(ind,i)=s23(ind,i)+D231hatFs(mm,bl,j)*Khats(16,i,j)&
                     & +D232hatFs(mm,bl,j)*Khats(17,i,j)+D233hatFs(mm,bl,j)*Khats(18,i,j)    
                s31(ind,i)=s31(ind,i)+D311hatFs(mm,bl,j)*Khats(19,i,j)&
                     & +D312hatFs(mm,bl,j)*Khats(20,i,j)+D313hatFs(mm,bl,j)*Khats(21,i,j)
                s32(ind,i)=s32(ind,i)+D321hatFs(mm,bl,j)*Khats(22,i,j)&
                     & +D322hatFs(mm,bl,j)*Khats(23,i,j)+D323hatFs(mm,bl,j)*Khats(24,i,j)
                s33(ind,i)=s33(ind,i)+D331hatFs(mm,bl,j)*Khats(25,i,j)&
                     & +D332hatFs(mm,bl,j)*Khats(26,i,j)+D333hatFs(mm,bl,j)*Khats(27,i,j)
             end if
             
             !domain I
             ind=m+mip
             !ind=m+mii
             
             if ((r-Ds*0.50_8)/beta>(r+Ds*0.50_8)/alpha) then
                do nq=0,Nqmax
                   dm=ceiling(((rj+Ds*0.50_8)/alpha+((rj-Ds*0.50_8)/beta-(rj+Ds*0.50_8)/alpha)*nq/(Nqmax*1.0_8))/Dt)
                   if (m-dm>=1.and.1<=ind.and.ind<=Mmax) then
                      dd1=Ds1(m-dm,j); dd2=Ds2(m-dm,j); dd3=Ds3(m-dm,j)
                      s11(ind,i)=s11(ind,i)+dd1*krnInsmp_lis(nq,1,i,j)+dd2*krnInsmp_lis(nq,2,i,j)+dd3*krnInsmp_lis(nq,3,i,j)
                      s12(ind,i)=s12(ind,i)+dd1*krnInsmp_lis(nq,4,i,j)+dd2*krnInsmp_lis(nq,5,i,j)+dd3*krnInsmp_lis(nq,6,i,j)
                      s13(ind,i)=s13(ind,i)+dd1*krnInsmp_lis(nq,7,i,j)+dd2*krnInsmp_lis(nq,8,i,j)+dd3*krnInsmp_lis(nq,9,i,j)
                      s21(ind,i)=s21(ind,i)+dd1*krnInsmp_lis(nq,10,i,j)&
                           & +dd2*krnInsmp_lis(nq,11,i,j)+dd3*krnInsmp_lis(nq,12,i,j)
                      s22(ind,i)=s22(ind,i)+dd1*krnInsmp_lis(nq,13,i,j)&
                           & +dd2*krnInsmp_lis(nq,14,i,j)+dd3*krnInsmp_lis(nq,15,i,j)
                      s23(ind,i)=s23(ind,i)+dd1*krnInsmp_lis(nq,16,i,j)&
                           & +dd2*krnInsmp_lis(nq,17,i,j)+dd3*krnInsmp_lis(nq,18,i,j)
                      s31(ind,i)=s31(ind,i)+dd1*krnInsmp_lis(nq,19,i,j)&
                           & +dd2*krnInsmp_lis(nq,20,i,j)+dd3*krnInsmp_lis(nq,21,i,j)
                      s32(ind,i)=s32(ind,i)+dd1*krnInsmp_lis(nq,22,i,j)&
                           & +dd2*krnInsmp_lis(nq,23,i,j)+dd3*krnInsmp_lis(nq,24,i,j)
                      s33(ind,i)=s33(ind,i)+dd1*krnInsmp_lis(nq,25,i,j)&
                           & +dd2*krnInsmp_lis(nq,26,i,j)+dd3*krnInsmp_lis(nq,27,i,j)
                   end if
                end do
             end if
             
             !domain S
             inds=ceiling((rj+Ds*0.50_8)/(beta*Dt))
             ind=m+mis
                
             if (1<=ind.and.ind<=Mmax.and.inds<=m-Mst) then
                dd1=Ds1(m-inds,j)
                dd2=Ds2(m-inds,j)
                dd3=Ds3(m-inds,j)
                s11(ind,i)=s11(ind,i)+dd1*ks_lis(1,i,j)+dd2*ks_lis(2,i,j)+dd3*ks_lis(3,i,j)
                s12(ind,i)=s12(ind,i)+dd1*ks_lis(4,i,j)+dd2*ks_lis(5,i,j)+dd3*ks_lis(6,i,j)
                s13(ind,i)=s13(ind,i)+dd1*ks_lis(7,i,j)+dd2*ks_lis(8,i,j)+dd3*ks_lis(9,i,j)
                s21(ind,i)=s21(ind,i)+dd1*ks_lis(10,i,j)+dd2*ks_lis(11,i,j)+dd3*ks_lis(12,i,j)
                s22(ind,i)=s22(ind,i)+dd1*ks_lis(13,i,j)+dd2*ks_lis(14,i,j)+dd3*ks_lis(15,i,j)
                s23(ind,i)=s23(ind,i)+dd1*ks_lis(16,i,j)+dd2*ks_lis(17,i,j)+dd3*ks_lis(18,i,j)
                s31(ind,i)=s31(ind,i)+dd1*ks_lis(19,i,j)+dd2*ks_lis(20,i,j)+dd3*ks_lis(21,i,j)
                s32(ind,i)=s32(ind,i)+dd1*ks_lis(22,i,j)+dd2*ks_lis(23,i,j)+dd3*ks_lis(24,i,j)
                s33(ind,i)=s33(ind,i)+dd1*ks_lis(25,i,j)+dd2*ks_lis(26,i,j)+dd3*ks_lis(27,i,j)
             end if
           end do
        end do
    end do
    keisu=-mu*0.50_8/beta
    do i=1,Nmax
       do m=1,Mmax
          s11(m,i)=keisu*s11(m,i)
          s12(m,i)=keisu*s12(m,i)
          s13(m,i)=keisu*s13(m,i)
          s21(m,i)=keisu*s21(m,i)
          s22(m,i)=keisu*s22(m,i)
          s23(m,i)=keisu*s23(m,i)
          s31(m,i)=keisu*s31(m,i)
          s32(m,i)=keisu*s32(m,i)
          s33(m,i)=keisu*s33(m,i)
       end do
    end do
          
  end subroutine fdph_solve !new

  ! time dependent part of Kernel of Fp and Fs
  subroutine htmFp(istar,j,dm,krnFp_lis,htmFp_lis) !before you call this, you have to call Kdyns_all and krnFp_all
    implicit none
    integer, intent(in) :: istar,j,dm
    real(8), dimension(27), intent(in) :: krnFp_lis
    real(8), dimension(27), intent(inout) :: htmFp_lis
    real(8) x,y,z,t,r,ri1,ri2,ma
    integer :: k
    real(8) :: val(27)

    x=xg(istar)-xg(j)
    y=yg(istar)-yg(j)
    z=zg(istar)-zg(j)
    r=sqrt(x*x+y*y+z*z)
    t=r/alpha+dm*Dt-Ds/alpha*0.50_8
    !call Kdyns_all(xg(istar),yg(istar),zg(istar),xa(j),ya(j),xb(j),yb(j),xc(j),yc(j),za(j))
      
    val(1)=Kdynf111(t,Dt)
    val(2)=Kdynf112(t,Dt)
    val(3)=Kdynf113(t,Dt)
    val(4)=Kdynf121(t,Dt)
    val(5)=Kdynf122(t,Dt)
    val(6)=Kdynf123(t,Dt)
    val(7)=Kdynf131(t,Dt)
    val(8)=Kdynf132(t,Dt)
    val(9)=Kdynf133(t,Dt)
    val(10)=Kdynf211(t,Dt)
    val(11)=Kdynf212(t,Dt)
    val(12)=Kdynf213(t,Dt)
    val(13)=Kdynf221(t,Dt)
    val(14)=Kdynf222(t,Dt)
    val(15)=Kdynf223(t,Dt)
    val(16)=Kdynf231(t,Dt)
    val(17)=Kdynf232(t,Dt)
    val(18)=Kdynf233(t,Dt)
    val(19)=Kdynf311(t,Dt)
    val(20)=Kdynf312(t,Dt)
    val(21)=Kdynf313(t,Dt)
    val(22)=Kdynf321(t,Dt)
    val(23)=Kdynf322(t,Dt)
    val(24)=Kdynf323(t,Dt)
    val(25)=Kdynf331(t,Dt)
    val(26)=Kdynf332(t,Dt)
    val(27)=Kdynf333(t,Dt)
          
    do k=1,27
       if (abs(krnFp_lis(k))<eps) then
          htmFp_lis(k)=0.0_8
       else
          htmFp_lis(k)=val(k)/krnFp_lis(k)
       end if
    end do 
  end subroutine htmFp

  subroutine htmFs(istar,j,dm,krnFs_lis,htmFs_lis)
    implicit none
    integer, intent(in) :: istar,j,dm
    real(8), dimension(27), intent(in) :: krnFs_lis
    real(8), dimension(27), intent(inout) :: htmFs_lis
    real(8) :: x,y,z,t,r
    integer :: k
    real(8) :: val(27)
    
    x=xg(istar)-xg(j)
    y=yg(istar)-yg(j)
    z=zg(istar)-zg(j)
    r=sqrt(x*x+y*y+z*z)
    t=r/beta+dm*Dt-Ds/beta*0.50_8
 
    val(1)=Kdynf111(t,Dt)
    val(2)=Kdynf112(t,Dt)
    val(3)=Kdynf113(t,Dt)
    val(4)=Kdynf121(t,Dt)
    val(5)=Kdynf122(t,Dt)
    val(6)=Kdynf123(t,Dt)
    val(7)=Kdynf131(t,Dt)
    val(8)=Kdynf132(t,Dt)
    val(9)=Kdynf133(t,Dt)
    val(10)=Kdynf211(t,Dt)
    val(11)=Kdynf212(t,Dt)
    val(12)=Kdynf213(t,Dt)
    val(13)=Kdynf221(t,Dt)
    val(14)=Kdynf222(t,Dt)
    val(15)=Kdynf223(t,Dt)
    val(16)=Kdynf231(t,Dt)
    val(17)=Kdynf232(t,Dt)
    val(18)=Kdynf233(t,Dt)
    val(19)=Kdynf311(t,Dt)
    val(20)=Kdynf312(t,Dt)
    val(21)=Kdynf313(t,Dt)
    val(22)=Kdynf321(t,Dt)
    val(23)=Kdynf322(t,Dt)
    val(24)=Kdynf323(t,Dt)
    val(25)=Kdynf331(t,Dt)
    val(26)=Kdynf332(t,Dt)
    val(27)=Kdynf333(t,Dt)
    do k=1,27
       if (abs(krnFs_lis(k))==0.0_8) then !<eps) then !eps=10**(-5.0_8)
          htmFs_lis(k)=0.0_8
       else
          htmFs_lis(k)=val(k)/krnFs_lis(k)
       end if
    end do  
  end subroutine htmFs

  function cal_h(K,Khat)
    implicit none
    real(8) :: K,Khat,cal_h
    if (abs(Khat)<=0.0_8) then !eps) then
       cal_h=0.0_8
    else
       cal_h=K/Khat
    end if
  end function cal_h

  subroutine krnInsmp(i,j,r,krnInsmp_lis) !calculate kernel of Domain I and store in KrnIn_lis
    implicit none
    integer, intent(in) :: i,j
    real(8), intent(in) :: r
    real(8), intent(inout) :: krnInsmp_lis(0:Nqmax,27,Nmax,Nmax)
    integer :: k,nq
    real(8) :: ts,te,dtofi,t
    real(8), dimension(27,Nqmax) :: ki
    
    ts=(r+Ds*0.50_8)/alpha
    te=(r-Ds*0.50_8)/beta
    if (ts<=te) then
       dtofi=te-ts
       do k=1,Nqmax
          t=ts+dtofi*(2*k-1)*0.50_8/(Nqmax)
          ki(1,k)=Kdynf111(t,Dt); ki(2,k)=Kdynf112(t,Dt); ki(3,k)=Kdynf113(t,Dt)
          ki(4,k)=Kdynf121(t,Dt); ki(5,k)=Kdynf122(t,Dt); ki(6,k)=Kdynf123(t,Dt)
          ki(7,k)=Kdynf131(t,Dt); ki(8,k)=Kdynf132(t,Dt); ki(9,k)=Kdynf133(t,Dt)
          ki(10,k)=Kdynf211(t,Dt); ki(11,k)=Kdynf212(t,Dt); ki(12,k)=Kdynf213(t,Dt)
          ki(13,k)=Kdynf221(t,Dt); ki(14,k)=Kdynf222(t,Dt); ki(15,k)=Kdynf223(t,Dt)
          ki(16,k)=Kdynf231(t,Dt); ki(17,k)=Kdynf232(t,Dt); ki(18,k)=Kdynf233(t,Dt)
          ki(19,k)=Kdynf311(t,Dt); ki(20,k)=Kdynf312(t,Dt); ki(21,k)=Kdynf313(t,Dt)
          ki(22,k)=Kdynf321(t,Dt); ki(23,k)=Kdynf322(t,Dt); ki(24,k)=Kdynf323(t,Dt)
          ki(25,k)=Kdynf331(t,Dt); ki(26,k)=Kdynf332(t,Dt); ki(27,k)=Kdynf333(t,Dt)
       end do
       !nq=0 & nq=Nqmax
       do k=1,27
          krnInsmp_lis(0,k,i,j)=-ki(k,1)
          krnInsmp_lis(Nqmax,k,i,j)=ki(k,Nqmax)
       end do
       !nq=1 to Nqmax-1
       do k=1,27
          do nq=1,Nqmax-1
             krnInsmp_lis(nq,k,i,j)=ki(k,nq)-ki(k,nq+1)
          end do
       end do
    end if
  end subroutine krnInsmp

  subroutine calc_ks(i,j,r,ks_lis)
    implicit none
    integer, intent(in) :: i,j
    real(8), intent(in) :: r
    real(8), intent(inout) :: ks_lis(:,:,:)
    integer :: inds
    real(8) :: t

    inds=min(ceiling((r+Ds*0.50_8)/(beta*Dt)+1),Mmax)
    t=inds*Dt
    ks_lis(1,i,j)=Kdynf111(t,Dt)
    ks_lis(2,i,j)=Kdynf112(t,Dt)
    ks_lis(3,i,j)=Kdynf113(t,Dt)
    ks_lis(4,i,j)=Kdynf121(t,Dt)
    ks_lis(5,i,j)=Kdynf122(t,Dt)
    ks_lis(6,i,j)=Kdynf123(t,Dt)
    ks_lis(7,i,j)=Kdynf131(t,Dt)
    ks_lis(8,i,j)=Kdynf132(t,Dt)
    ks_lis(9,i,j)=Kdynf133(t,Dt)
    ks_lis(10,i,j)=Kdynf211(t,Dt)
    ks_lis(11,i,j)=Kdynf212(t,Dt)
    ks_lis(12,i,j)=Kdynf213(t,Dt)
    ks_lis(13,i,j)=Kdynf221(t,Dt)
    ks_lis(14,i,j)=Kdynf222(t,Dt)
    ks_lis(15,i,j)=Kdynf223(t,Dt)
    ks_lis(16,i,j)=Kdynf231(t,Dt)
    ks_lis(17,i,j)=Kdynf232(t,Dt)
    ks_lis(18,i,j)=Kdynf233(t,Dt)
    ks_lis(19,i,j)=Kdynf311(t,Dt)
    ks_lis(20,i,j)=Kdynf312(t,Dt)
    ks_lis(21,i,j)=Kdynf313(t,Dt)
    ks_lis(22,i,j)=Kdynf321(t,Dt)
    ks_lis(23,i,j)=Kdynf322(t,Dt)
    ks_lis(24,i,j)=Kdynf323(t,Dt)
    ks_lis(25,i,j)=Kdynf331(t,Dt)
    ks_lis(26,i,j)=Kdynf332(t,Dt)
    ks_lis(27,i,j)=Kdynf333(t,Dt)
  end subroutine calc_ks

  subroutine error()
    implicit none
    integer :: i,j,bl,m,dm,k
    real(8), dimension(Nmax,Nmax) :: dis
    real(8) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,r,t1
    real(8) :: xsg1,xsg2,xsg3,xs1star,xs2star,xs3star,xr1star,xr2star,xr3star
    !real(8) :: val_ex,val_app,gosa2 !val=Kernel(t)*D(m-t)
    character (len=2) :: dom
    real(8) :: tt !travel-time
    real(8) :: krnInsmp_lis(0:Nqmax,27,Nmax,Nmax)
    real(8) :: ks(27)
    real(8) :: dd1,dd2,dd3,keisu
    integer :: istar_lis(Nmax),istar, block_lis(Nmax), jstar !istar_ind
    integer :: ind,nq,mip,mis,mii,mi,mj,msmin,mpmax_now,msmax_now,mm
    integer :: inds
    integer :: msmin_lis(Nmax,Nmax),mip_lis(Nmax,nblock),mis_lis(Nmax,nblock)
    integer :: cal_was(nblock,Nmax), cal_now
    real(8) :: xry,yry,zry,rry,x,y,z,ri,rj,dri,rtrue,rstar
    real(8) :: val
    real(8), dimension(Mmax-1,Nmax,Nmax) :: K111,K112,K113,K121,K122,K123,K131,K132,K133
    real(8), dimension(Mmax-1,Nmax,Nmax) :: K211,K212,K213,K221,K222,K223,K231,K232,K233
    real(8), dimension(Mmax-1,Nmax,Nmax) :: K311,K312,K313,K321,K322,K323,K331,K332,K333
    real(8), dimension(27,Nmax,Nmax) :: Khatp,Khats
    real(8), dimension(0:mpmax,27,nblock,Nmax):: hp
    real(8), dimension(0:msmax-1,27,nblock,Nmax) :: hs
    real(8), dimension(27,Nmax,Nmax) :: ks_lis
    real(8), dimension(Mmax,Nmax,3) :: err11,ex11,app11,err31,ex31,app31
    real(8) :: val1,val2,val3,val4,val5,val6
 
    !===initialization===
    Ds1(:,:)=0.0_8; Ds2(:,:)=0.0_8; Ds3(:,:)=0.0_8
    call Dhat_ini()
    call calc_istar(nxb,nyb,nzb,istar_lis,block_lis,block_center) !determine representative points of observation
    krnInsmp_lis(:,:,:,:)=0.0_8
    hp(:,:,:,:)=0.0_8; hs(:,:,:,:)=0.0_8
    cal_was(:,:)=-1
    err11(:,:,:)=0.0_8; ex11(:,:,:)=0.0_8; app11(:,:,:)=0.0_8
    err31(:,:,:)=0.0_8; ex31(:,:,:)=0.0_8; app31(:,:,:)=0.0_8
    
    !write(*,*) block_lis
    call check_same(nxb,nyb,nzb,block_lis,block_center)

    !===calculate Kernel and store===
    do j=1,Nmax
       xsa1=xa(j); xsa2=ya(j); xsb1=xb(j); xsb2=yb(j)
       xsc1=xc(j); xsc2=yc(j); xs3=za(j)
       xsg1=xg(j); xsg2=yg(j); xsg3=zg(j)
       do i=1,Nmax
          cal_now=-1
          bl=block_lis(i)
          istar=block_center(bl)
          xr1star=xg(istar); xr2star=yg(istar); xr3star=zg(istar)
          xr1=xg(i); xr2=yg(i); xr3=zg(i)
          call Kdyns_all(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
          r=sqrt((xr1-xsg1)**2+(xr2-xsg2)**2+(xr3-xsg3)**2)
          rj=sqrt((xr1star-xsg1)**2+(xr2star-xsg2)**2+(xr3star-xsg3)**2)
          dis(i,j)=r
          !Khat (Fp)  
          Khatp(:,i,j)=0.0_8
          do dm=0,mpmax-1
             t=r/alpha+dm*Dt-Ds/alpha*0.50_8 !t can be negative, but kernel become 0 when t<0
             Khatp(1,i,j)=Khatp(1,i,j)+Kdynf111(t,Dt)
             Khatp(2,i,j)=Khatp(2,i,j)+Kdynf112(t,Dt)
             Khatp(3,i,j)=Khatp(3,i,j)+Kdynf113(t,Dt)
             Khatp(4,i,j)=Khatp(4,i,j)+Kdynf121(t,Dt)
             Khatp(5,i,j)=Khatp(5,i,j)+Kdynf122(t,Dt)
             Khatp(6,i,j)=Khatp(6,i,j)+Kdynf123(t,Dt)
             Khatp(7,i,j)=Khatp(7,i,j)+Kdynf131(t,Dt)
             Khatp(8,i,j)=Khatp(8,i,j)+Kdynf132(t,Dt)
             Khatp(9,i,j)=Khatp(9,i,j)+Kdynf133(t,Dt)
             Khatp(10,i,j)=Khatp(10,i,j)+Kdynf211(t,Dt)
             Khatp(11,i,j)=Khatp(11,i,j)+Kdynf212(t,Dt)
             Khatp(12,i,j)=Khatp(12,i,j)+Kdynf213(t,Dt)
             Khatp(13,i,j)=Khatp(13,i,j)+Kdynf221(t,Dt)
             Khatp(14,i,j)=Khatp(14,i,j)+Kdynf222(t,Dt)
             Khatp(15,i,j)=Khatp(15,i,j)+Kdynf223(t,Dt)
             Khatp(16,i,j)=Khatp(16,i,j)+Kdynf231(t,Dt)
             Khatp(17,i,j)=Khatp(17,i,j)+Kdynf232(t,Dt)
             Khatp(18,i,j)=Khatp(18,i,j)+Kdynf233(t,Dt)
             Khatp(19,i,j)=Khatp(19,i,j)+Kdynf311(t,Dt)
             Khatp(20,i,j)=Khatp(20,i,j)+Kdynf312(t,Dt)
             Khatp(21,i,j)=Khatp(21,i,j)+Kdynf313(t,Dt)
             Khatp(22,i,j)=Khatp(22,i,j)+Kdynf321(t,Dt)
             Khatp(23,i,j)=Khatp(23,i,j)+Kdynf322(t,Dt)
             Khatp(24,i,j)=Khatp(24,i,j)+Kdynf323(t,Dt)
             Khatp(25,i,j)=Khatp(25,i,j)+Kdynf331(t,Dt)
             Khatp(26,i,j)=Khatp(26,i,j)+Kdynf332(t,Dt)
             Khatp(27,i,j)=Khatp(27,i,j)+Kdynf333(t,Dt)
          end do
          !Khat (Fs)
          msmin=0
          t1=r/alpha+Ds/alpha*0.50_8
          if(t1>r/beta-Ds/beta*0.50_8) then
             msmin=floor((t1-(r/beta-Ds/beta*0.50_8))/Dt)
          end if
          msmin_lis(i,j)=msmin
          !t2=r/beta+msmin*Dt-Ds/beta*0.50_8 !t2>t1
          Khats(:,i,j)=0.0_8
          do dm=msmin,msmax-1,1
             t=r/beta+dm*Dt-Ds/beta*0.50_8
             Khats(1,i,j)=Khats(1,i,j)+Kdynf111(t,Dt)
             Khats(2,i,j)=Khats(2,i,j)+Kdynf112(t,Dt)
             Khats(3,i,j)=Khats(3,i,j)+Kdynf113(t,Dt)
             Khats(4,i,j)=Khats(4,i,j)+Kdynf121(t,Dt)
             Khats(5,i,j)=Khats(5,i,j)+Kdynf122(t,Dt)
             Khats(6,i,j)=Khats(6,i,j)+Kdynf123(t,Dt)
             Khats(7,i,j)=Khats(7,i,j)+Kdynf131(t,Dt)
             Khats(8,i,j)=Khats(8,i,j)+Kdynf132(t,Dt)
             Khats(9,i,j)=Khats(9,i,j)+Kdynf133(t,Dt)
             Khats(10,i,j)=Khats(10,i,j)+Kdynf211(t,Dt)
             Khats(11,i,j)=Khats(11,i,j)+Kdynf212(t,Dt)
             Khats(12,i,j)=Khats(12,i,j)+Kdynf213(t,Dt)
             Khats(13,i,j)=Khats(13,i,j)+Kdynf221(t,Dt)
             Khats(14,i,j)=Khats(14,i,j)+Kdynf222(t,Dt)
             Khats(15,i,j)=Khats(15,i,j)+Kdynf223(t,Dt)
             Khats(16,i,j)=Khats(16,i,j)+Kdynf231(t,Dt)
             Khats(17,i,j)=Khats(17,i,j)+Kdynf232(t,Dt)
             Khats(18,i,j)=Khats(18,i,j)+Kdynf233(t,Dt)
             Khats(19,i,j)=Khats(19,i,j)+Kdynf311(t,Dt)
             Khats(20,i,j)=Khats(20,i,j)+Kdynf312(t,Dt)
             Khats(21,i,j)=Khats(21,i,j)+Kdynf313(t,Dt)
             Khats(22,i,j)=Khats(22,i,j)+Kdynf321(t,Dt)
             Khats(23,i,j)=Khats(23,i,j)+Kdynf322(t,Dt)
             Khats(24,i,j)=Khats(24,i,j)+Kdynf323(t,Dt)
             Khats(25,i,j)=Khats(25,i,j)+Kdynf331(t,Dt)
             Khats(26,i,j)=Khats(26,i,j)+Kdynf332(t,Dt)
             Khats(27,i,j)=Khats(27,i,j)+Kdynf333(t,Dt)            
          end do
          !calculate h
          if (istar==i.and.cal_was(bl,j)==-1) then
             do dm=0,mpmax-1
                t=rj/alpha+dm*Dt-Ds/alpha*0.50_8
                hp(dm,1,bl,j)=cal_h(Kdynf111(t,Dt),Khatp(1,i,j))
                hp(dm,2,bl,j)=cal_h(Kdynf112(t,Dt),Khatp(2,i,j))
                hp(dm,3,bl,j)=cal_h(Kdynf113(t,Dt),Khatp(3,i,j))
                hp(dm,4,bl,j)=cal_h(Kdynf121(t,Dt),Khatp(4,i,j))
                hp(dm,5,bl,j)=cal_h(Kdynf122(t,Dt),Khatp(5,i,j))
                hp(dm,6,bl,j)=cal_h(Kdynf123(t,Dt),Khatp(6,i,j))
                hp(dm,7,bl,j)=cal_h(Kdynf131(t,Dt),Khatp(7,i,j))
                hp(dm,8,bl,j)=cal_h(Kdynf132(t,Dt),Khatp(8,i,j))
                hp(dm,9,bl,j)=cal_h(Kdynf133(t,Dt),Khatp(9,i,j))
                hp(dm,10,bl,j)=cal_h(Kdynf211(t,Dt),Khatp(10,i,j))
                hp(dm,11,bl,j)=cal_h(Kdynf212(t,Dt),Khatp(11,i,j))
                hp(dm,12,bl,j)=cal_h(Kdynf213(t,Dt),Khatp(12,i,j))
                hp(dm,13,bl,j)=cal_h(Kdynf221(t,Dt),Khatp(13,i,j))
                hp(dm,14,bl,j)=cal_h(Kdynf222(t,Dt),Khatp(14,i,j))
                hp(dm,15,bl,j)=cal_h(Kdynf223(t,Dt),Khatp(15,i,j))
                hp(dm,16,bl,j)=cal_h(Kdynf231(t,Dt),Khatp(16,i,j))
                hp(dm,17,bl,j)=cal_h(Kdynf232(t,Dt),Khatp(17,i,j))
                hp(dm,18,bl,j)=cal_h(Kdynf233(t,Dt),Khatp(18,i,j))
                hp(dm,19,bl,j)=cal_h(Kdynf311(t,Dt),Khatp(19,i,j))
                hp(dm,20,bl,j)=cal_h(Kdynf312(t,Dt),Khatp(20,i,j))
                hp(dm,21,bl,j)=cal_h(Kdynf313(t,Dt),Khatp(21,i,j))
                hp(dm,22,bl,j)=cal_h(Kdynf321(t,Dt),Khatp(22,i,j))
                hp(dm,23,bl,j)=cal_h(Kdynf322(t,Dt),Khatp(23,i,j))
                hp(dm,24,bl,j)=cal_h(Kdynf323(t,Dt),Khatp(24,i,j))
                hp(dm,25,bl,j)=cal_h(Kdynf331(t,Dt),Khatp(25,i,j))
                hp(dm,26,bl,j)=cal_h(Kdynf332(t,Dt),Khatp(26,i,j))
                hp(dm,27,bl,j)=cal_h(Kdynf333(t,Dt),Khatp(27,i,j))
             end do
             !msmin=0
             !t1=r/alpha+Ds/alpha*0.50_8
             !if(t1>=r/beta-Ds/beta*0.50_8) then
             !   msmin=floor((t1-(r/beta-Ds/beta*0.50_8))/Dt)
             !end if
             !t2=r/beta+msmin*Dt-Ds/beta*0.50_8 !t2>t1
             do dm=msmin,msmax-1,1
                t=rj/beta+dm*Dt-Ds/beta*0.50_8
                hs(dm,1,bl,j)=cal_h(Kdynf111(t,Dt),Khats(1,i,j))
                hs(dm,2,bl,j)=cal_h(Kdynf112(t,Dt),Khats(2,i,j))
                hs(dm,3,bl,j)=cal_h(Kdynf113(t,Dt),Khats(3,i,j))
                hs(dm,4,bl,j)=cal_h(Kdynf121(t,Dt),Khats(4,i,j))
                hs(dm,5,bl,j)=cal_h(Kdynf122(t,Dt),Khats(5,i,j))
                hs(dm,6,bl,j)=cal_h(Kdynf123(t,Dt),Khats(6,i,j))
                hs(dm,7,bl,j)=cal_h(Kdynf131(t,Dt),Khats(7,i,j))
                hs(dm,8,bl,j)=cal_h(Kdynf132(t,Dt),Khats(8,i,j))
                hs(dm,9,bl,j)=cal_h(Kdynf133(t,Dt),Khats(9,i,j))
                hs(dm,10,bl,j)=cal_h(Kdynf211(t,Dt),Khats(10,i,j))
                hs(dm,11,bl,j)=cal_h(Kdynf212(t,Dt),Khats(11,i,j))
                hs(dm,12,bl,j)=cal_h(Kdynf213(t,Dt),Khats(12,i,j))
                hs(dm,13,bl,j)=cal_h(Kdynf221(t,Dt),Khats(13,i,j))
                hs(dm,14,bl,j)=cal_h(Kdynf222(t,Dt),Khats(14,i,j))
                hs(dm,15,bl,j)=cal_h(Kdynf223(t,Dt),Khats(15,i,j))
                hs(dm,16,bl,j)=cal_h(Kdynf231(t,Dt),Khats(16,i,j))
                hs(dm,17,bl,j)=cal_h(Kdynf232(t,Dt),Khats(17,i,j))
                hs(dm,18,bl,j)=cal_h(Kdynf233(t,Dt),Khats(18,i,j))
                hs(dm,19,bl,j)=cal_h(Kdynf311(t,Dt),Khats(19,i,j))
                hs(dm,20,bl,j)=cal_h(Kdynf312(t,Dt),Khats(20,i,j))
                hs(dm,21,bl,j)=cal_h(Kdynf313(t,Dt),Khats(21,i,j))
                hs(dm,22,bl,j)=cal_h(Kdynf321(t,Dt),Khats(22,i,j))
                hs(dm,23,bl,j)=cal_h(Kdynf322(t,Dt),Khats(23,i,j))
                hs(dm,24,bl,j)=cal_h(Kdynf323(t,Dt),Khats(24,i,j))
                hs(dm,25,bl,j)=cal_h(Kdynf331(t,Dt),Khats(25,i,j))
                hs(dm,26,bl,j)=cal_h(Kdynf332(t,Dt),Khats(26,i,j))
                hs(dm,27,bl,j)=cal_h(Kdynf333(t,Dt),Khats(27,i,j))
             end do
             cal_was(bl,j)=1
          end if
          !calculate kernel of Domain I and store it
          call krnInsmp(i,j,r,krnInsmp_lis)
          !calculate kernel of Domain S and store it
          call calc_ks(i,j,r,ks_lis)
       end do
    end do

    !To calculate mi
    do bl=1,nblock
       jstar=block_center(bl)
       xs1star=xg(jstar); xs2star=yg(jstar); xs3star=zg(jstar)
       do i=1,Nmax
          istar=istar_lis(i)
          !call calc_ri_31(istar,i,jstar,ri,mip,mis)
          xr1=xg(i); xr2=yg(i); xr3=zg(i)
          xr1star=xg(istar); xr2star=yg(istar); xr3star=zg(istar)
          ri=sqrt((xr1-xs1star)**2+(xr2-xs2star)**2+(xr3-xs3star)**2)
          rstar=sqrt((xr1star-xs1star)**2+(xr2star-xs2star)**2+(xr3star-xs3star)**2)
          dri=ri-rstar
          mip=floor(dri/(alpha*Dt))
          mis=floor(dri/(beta*Dt))
          mip_lis(i,bl)=mip
          mis_lis(i,bl)=mis
       end do
    end do
           
    ! !===time marching===
    do m=Mst+1,Mmax
       if (mod(m,10)==0) then
          write(*,*) "now: m=",m
       end if
       !calculate Ds and Dhat
       do j=1,Nmax
          Ds1(m,j)=Ds1(m-1,j)+D1(m,j)
          Ds2(m,j)=Ds2(m-1,j)+D2(m,j)
          Ds3(m,j)=Ds3(m-1,j)+D3(m,j)
          !xsa1=xa(j); xsa2=ya(j); xsb1=xb(j); xsb2=yb(j)
          !xsc1=xc(j); xsc2=yc(j); xs3=za(j)
          !xsg1=xg(j); xsg2=yg(j); xsg3=zg(j)
          do bl=1,nblock
             istar=block_center(bl)
             rj=dis(istar,j)
             !xr1star=xg(istar); xr2star=yg(istar); xr3star=zg(istar)
             !xr1=xg(i); xr2=yg(i); xr3=zg(i)
             !call Kdyns_all(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
             !rtrue=sqrt((xr1-xsg1)**2+(xr2-xsg2)**2+(xr3-xsg3)**2)
             !call krnFp_all(i,j,rtrue,krnFp_lis)
             !call krnFs_all(i,j,rtrue,krnFs_lis
             do dm=0,mpmax-1
                ind=m-dm !Sato & Ando (52)
                if (ind<1 .or. ind>=m) then
                   cycle
                end if
                dd1=D1(ind,j); dd2=D2(ind,j); dd3=D3(ind,j)
                D111hatFp(m,bl,j)=D111hatFp(m,bl,j)+dd1*hp(dm,1,bl,j)
                D112hatFp(m,bl,j)=D112hatFp(m,bl,j)+dd2*hp(dm,2,bl,j)
                D113hatFp(m,bl,j)=D113hatFp(m,bl,j)+dd3*hp(dm,3,bl,j)
                D121hatFp(m,bl,j)=D121hatFp(m,bl,j)+dd1*hp(dm,4,bl,j)
                D122hatFp(m,bl,j)=D122hatFp(m,bl,j)+dd2*hp(dm,5,bl,j)
                D123hatFp(m,bl,j)=D123hatFp(m,bl,j)+dd3*hp(dm,6,bl,j)
                D131hatFp(m,bl,j)=D131hatFp(m,bl,j)+dd1*hp(dm,7,bl,j)
                D132hatFp(m,bl,j)=D132hatFp(m,bl,j)+dd2*hp(dm,8,bl,j)
                D133hatFp(m,bl,j)=D133hatFp(m,bl,j)+dd3*hp(dm,9,bl,j)
                D211hatFp(m,bl,j)=D211hatFp(m,bl,j)+dd1*hp(dm,10,bl,j)
                D212hatFp(m,bl,j)=D212hatFp(m,bl,j)+dd2*hp(dm,11,bl,j)
                D213hatFp(m,bl,j)=D213hatFp(m,bl,j)+dd3*hp(dm,12,bl,j)
                D221hatFp(m,bl,j)=D221hatFp(m,bl,j)+dd1*hp(dm,13,bl,j)
                D222hatFp(m,bl,j)=D222hatFp(m,bl,j)+dd2*hp(dm,14,bl,j)
                D223hatFp(m,bl,j)=D223hatFp(m,bl,j)+dd3*hp(dm,15,bl,j)
                D231hatFp(m,bl,j)=D231hatFp(m,bl,j)+dd1*hp(dm,16,bl,j)
                D232hatFp(m,bl,j)=D232hatFp(m,bl,j)+dd2*hp(dm,17,bl,j)
                D233hatFp(m,bl,j)=D233hatFp(m,bl,j)+dd3*hp(dm,18,bl,j)
                D311hatFp(m,bl,j)=D311hatFp(m,bl,j)+dd1*hp(dm,19,bl,j)
                D312hatFp(m,bl,j)=D312hatFp(m,bl,j)+dd2*hp(dm,20,bl,j)
                D313hatFp(m,bl,j)=D313hatFp(m,bl,j)+dd3*hp(dm,21,bl,j)
                D321hatFp(m,bl,j)=D321hatFp(m,bl,j)+dd1*hp(dm,22,bl,j)
                D322hatFp(m,bl,j)=D322hatFp(m,bl,j)+dd2*hp(dm,23,bl,j)
                D323hatFp(m,bl,j)=D323hatFp(m,bl,j)+dd3*hp(dm,24,bl,j)
                D331hatFp(m,bl,j)=D331hatFp(m,bl,j)+dd1*hp(dm,25,bl,j)
                D332hatFp(m,bl,j)=D332hatFp(m,bl,j)+dd2*hp(dm,26,bl,j)
                D333hatFp(m,bl,j)=D333hatFp(m,bl,j)+dd3*hp(dm,27,bl,j)
             end do

             msmin=0
             t1=rj/alpha+Ds/alpha*0.50_8
             if(t1>rj/beta-Ds/beta*0.50_8) then
                msmin=floor((t1-(rj/beta-Ds/beta*0.50_8))/Dt)
             end if
             do dm=msmin,msmax-1
                ind=m-dm
                if (ind<1 .or. ind>=m) then
                   cycle
                end if
                !call htmFs(istar,j,dm,krnFs_lis,htmFs_lis)
                dd1=D1(ind,j); dd2=D2(ind,j); dd3=D3(ind,j)
                D111hatFs(m,bl,j)=D111hatFs(m,bl,j)+dd1*hs(dm,1,bl,j)
                D112hatFs(m,bl,j)=D112hatFs(m,bl,j)+dd2*hs(dm,2,bl,j)
                D113hatFs(m,bl,j)=D113hatFs(m,bl,j)+dd3*hs(dm,3,bl,j)
                D121hatFs(m,bl,j)=D121hatFs(m,bl,j)+dd1*hs(dm,4,bl,j)
                D122hatFs(m,bl,j)=D122hatFs(m,bl,j)+dd2*hs(dm,5,bl,j)
                D123hatFs(m,bl,j)=D123hatFs(m,bl,j)+dd3*hs(dm,6,bl,j)
                D131hatFs(m,bl,j)=D131hatFs(m,bl,j)+dd1*hs(dm,7,bl,j)
                D132hatFs(m,bl,j)=D132hatFs(m,bl,j)+dd2*hs(dm,8,bl,j)
                D133hatFs(m,bl,j)=D133hatFs(m,bl,j)+dd3*hs(dm,9,bl,j)
                D211hatFs(m,bl,j)=D211hatFs(m,bl,j)+dd1*hs(dm,10,bl,j)
                D212hatFs(m,bl,j)=D212hatFs(m,bl,j)+dd2*hs(dm,11,bl,j)
                D213hatFs(m,bl,j)=D213hatFs(m,bl,j)+dd3*hs(dm,12,bl,j)
                D221hatFs(m,bl,j)=D221hatFs(m,bl,j)+dd1*hs(dm,13,bl,j)
                D222hatFs(m,bl,j)=D222hatFs(m,bl,j)+dd2*hs(dm,14,bl,j)
                D223hatFs(m,bl,j)=D223hatFs(m,bl,j)+dd3*hs(dm,15,bl,j)
                D231hatFs(m,bl,j)=D231hatFs(m,bl,j)+dd1*hs(dm,16,bl,j)
                D232hatFs(m,bl,j)=D232hatFs(m,bl,j)+dd2*hs(dm,17,bl,j)
                D233hatFs(m,bl,j)=D233hatFs(m,bl,j)+dd3*hs(dm,18,bl,j)
                D311hatFs(m,bl,j)=D311hatFs(m,bl,j)+dd1*hs(dm,19,bl,j)
                D312hatFs(m,bl,j)=D312hatFs(m,bl,j)+dd2*hs(dm,20,bl,j)
                D313hatFs(m,bl,j)=D313hatFs(m,bl,j)+dd3*hs(dm,21,bl,j)
                D321hatFs(m,bl,j)=D321hatFs(m,bl,j)+dd1*hs(dm,22,bl,j)
                D322hatFs(m,bl,j)=D322hatFs(m,bl,j)+dd2*hs(dm,23,bl,j)
                D323hatFs(m,bl,j)=D323hatFs(m,bl,j)+dd3*hs(dm,24,bl,j)
                D331hatFs(m,bl,j)=D331hatFs(m,bl,j)+dd1*hs(dm,25,bl,j)
                D332hatFs(m,bl,j)=D332hatFs(m,bl,j)+dd2*hs(dm,26,bl,j)
                D333hatFs(m,bl,j)=D333hatFs(m,bl,j)+dd3*hs(dm,27,bl,j)
             end do
          end do
       end do
       !convolution start !!!
       do i=1,Nmax !receiver
          xr1=xg(i)
          xr2=yg(i)
          xr3=zg(i)
          istar=istar_lis(i)
          bl=block_lis(i)
          xr1star=xg(istar)
          xr2star=yg(istar)
          xr3star=zg(istar)
          do j=1,Nmax !source
             xsa1=xa(j)
             xsa2=ya(j)
             xsb1=xb(j)
             xsb2=yb(j)
             xsc1=xc(j)
             xsc2=yc(j)
             xs3=za(j) !=zb(j)=zc(j)
             xsg1=xg(j); xsg2=yg(j); xsg3=zg(j)
             r=dis(i,j)
             !computation of mi
             rj=dis(istar,j) !sqrt((xr1star-xsg1)**2+(xr2star-xsg2)**2+(xr3star-xsg3)**2)
             mj=max(ceiling((rj-Ds*0.50_8)/(alpha*Dt)),0)
             jstar=istar_lis(j)
             mip=mip_lis(i,block_lis(jstar))
             mis=mis_lis(i,block_lis(jstar))
     
             !domain Fp
             !mj=max(ceiling(rj/(alpha*Dt)-Ds/(alpha*Dt)*0.5),1)
             ind=m+mip
             mm=m-mj
             
             if (1<=ind.and.ind<=Mmax.and.mm>=1) then
                app11(m,i,1)=app11(m,i,1)+D111hatFp(mm,bl,j)*Khatp(1,i,j)&
                     & +D112hatFp(mm,bl,j)*Khatp(2,i,j)+D113hatFp(mm,bl,j)*Khatp(3,i,j)
                app31(m,i,1)=app31(m,i,1)+D311hatFp(mm,bl,j)*Khatp(19,i,j)&
                     & +D312hatFp(mm,bl,j)*Khatp(20,i,j)+D313hatFp(mm,bl,j)*Khatp(21,i,j)
             end if
             !domain Fs
             !mj=max(ceiling(rj/(beta*Dt)-Ds/(beta*Dt)*0.5),1)
             ind=m+mis
             if (1<=ind.and.ind<=Mmax.and.mm>=1) then
                app11(m,i,1)=app11(m,i,1)+D111hatFs(mm,bl,j)*Khats(1,i,j)&
                     & +D112hatFs(mm,bl,j)*Khats(2,i,j)+D113hatFs(mm,bl,j)*Khats(3,i,j)
                app31(m,i,1)=app31(m,i,1)+D311hatFs(mm,bl,j)*Khats(19,i,j)&
                     & +D312hatFs(mm,bl,j)*Khats(20,i,j)+D313hatFs(mm,bl,j)*Khats(21,i,j)
             end if
             
             !domain I
             ind=m+mip
             !ind=m+mii
             
             if ((r-Ds*0.50_8)/beta>(r+Ds*0.50_8)/alpha) then
                do nq=0,Nqmax
                   dm=ceiling(((rj+Ds*0.50_8)/alpha+((rj-Ds*0.50_8)/beta-(rj+Ds*0.50_8)/alpha)*nq/(Nqmax*1.0_8))/Dt)
                   if (m-dm>=1.and.1<=ind.and.ind<=Mmax) then
                      dd1=Ds1(m-dm,j); dd2=Ds2(m-dm,j); dd3=Ds3(m-dm,j)
                      app11(m,i,2)=app11(m,i,2)+dd1*krnInsmp_lis(nq,1,i,j)+dd2*krnInsmp_lis(nq,2,i,j)+dd3*krnInsmp_lis(nq,3,i,j)
                      app31(m,i,2)=app31(m,i,2)+dd1*krnInsmp_lis(nq,19,i,j) &
                           & +dd2*krnInsmp_lis(nq,20,i,j)+dd3*krnInsmp_lis(nq,21,i,j)
                   end if
                end do
             end if
             
             !domain S
             inds=ceiling((rj+Ds*0.50_8)/(beta*Dt))
             ind=m+mis
                
             if (1<=ind.and.ind<=Mmax.and.inds<=m-Mst) then
                dd1=Ds1(m-inds,j); dd2=Ds2(m-inds,j); dd3=Ds3(m-inds,j)
                app11(m,i,3)=app11(m,i,3)+dd1*ks_lis(1,i,j)+dd2*ks_lis(2,i,j)+dd3*ks_lis(3,i,j)
                app31(m,i,3)=app31(m,i,3)+dd1*ks_lis(19,i,j)+dd2*ks_lis(20,i,j)+dd3*ks_lis(21,i,j)
             end if
           end do
        end do
    end do
    keisu=-mu*0.50_8/beta
    app11=app11*keisu
    app31=app31*keisu

    !ex
      do i=1,Nmax !receiver   
       xr1=xg(i)
       xr2=yg(i)
       xr3=zg(i)
       do j=1,Nmax !source
          xsa1=xa(j)
          xsa2=ya(j)
          xsb1=xb(j)
          xsb2=yb(j)
          xsc1=xc(j)
          xsc2=yc(j)
          xs3=za(j) !=zb(j)=zc(j)
          call Kdyns_all(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
          do dm=1,Mmax-1
             t=Dt*dm
             K111(dm,i,j)=Kdynf111(t,Dt)
             K112(dm,i,j)=Kdynf112(t,Dt)
             K113(dm,i,j)=Kdynf113(t,Dt)
             K311(dm,i,j)=Kdynf311(t,Dt)
             K312(dm,i,j)=Kdynf312(t,Dt)
             K313(dm,i,j)=Kdynf313(t,Dt)
          end do
       end do
    end do
    do i=1,Nmax
       do j=1,Nmax
          r=dis(i,j)
          do m=1,Mmax
             do dm=1,m-1
                t=dm*Dt
                call domain(r,t,dom)
                if (dom=="Fp".or.dom=="Fs") then
                   ex11(m,i,1)=ex11(m,i,1)+K111(dm,i,j)*D1(m-dm,j)+K112(dm,i,j)*D2(m-dm,j)+K113(dm,i,j)*D3(m-dm,j)
                   ex31(m,i,1)=ex31(m,i,1)+K311(dm,i,j)*D1(m-dm,j)+K312(dm,i,j)*D2(m-dm,j)+K313(dm,i,j)*D3(m-dm,j)
                else if (dom=="I ") then
                   ex11(m,i,2)=ex11(m,i,2)+K111(dm,i,j)*D1(m-dm,j)+K112(dm,i,j)*D2(m-dm,j)+K113(dm,i,j)*D3(m-dm,j)
                   ex31(m,i,2)=ex31(m,i,2)+K311(dm,i,j)*D1(m-dm,j)+K312(dm,i,j)*D2(m-dm,j)+K313(dm,i,j)*D3(m-dm,j)
                else if (dom=="S ") then
                   ex11(m,i,3)=ex11(m,i,3)+K111(dm,i,j)*D1(m-dm,j)+K112(dm,i,j)*D2(m-dm,j)+K113(dm,i,j)*D3(m-dm,j)
                   ex31(m,i,3)=ex31(m,i,3)+K311(dm,i,j)*D1(m-dm,j)+K312(dm,i,j)*D2(m-dm,j)+K313(dm,i,j)*D3(m-dm,j)
                end if
             end do
          end do
       end do
    end do
    
    keisu=-mu*0.50_8/beta
    ! f11ex=f11ex*keisu
    ! f31ex=f31ex*keisu
    ! i11ex=i11ex*keisu
    ! i31ex=i31ex*keisu
    ! s11ex=s11ex*keisu
    ! s31ex=s31ex*keisu
    ex11=ex11*keisu
    ex31=ex31*keisu
    
    !open(18,file="f_err.dat",status="replace")
    !open(19,file="i_err.dat",status="replace")
    !open(20,file="s_err.dat",status="replace")
    open(21,file="err11.dat",status="replace")
    open(22,file="err31.dat",status="replace")

    do m=1,Mmax
       val1=0.0_8; val2=0.0_8; val3=0.0_8
       val4=0.0_8; val5=0.0_8; val6=0.0_8
       do i=1,1,Nmax
          do k=1,3
             err11(m,i,k)=err11(m,i,k)+(ex11(m,i,k)-app11(m,i,k))**2
             err31(m,i,k)=err31(m,i,k)+(ex31(m,i,k)-app31(m,i,k))**2
             val1=val1+err11(m,i,1)
             val2=val2+err11(m,i,2)
             val3=val3+err11(m,i,3)
             val4=val4+err31(m,i,1)
             val5=val5+err31(m,i,2)
             val6=val6+err31(m,i,3)
          end do
       end do
       write(21,*) m,val1,val2,val3
       write(22,*) m,val4,val5,val6
    end do
    !write(21,*) f11ex,f11app,i11ex,i11app,s11ex,s11app
    !write(22,*) f31ex,f31app,i31ex,i31app,s31ex,s31app           
    
    !close(18)
    !close(19)
    !close(20)
    close(21)
    close(22)
  end subroutine error

  subroutine domain(r,t,dom)
    implicit none
    real(8), intent(in) :: r,t
    character(len=2), intent(out) :: dom
    real(8) :: tp1,tp2,ts1,ts2
    tp1=(r-Ds*0.50_8)/alpha
    tp2=(r+Ds*0.50_8)/alpha
    ts1=(r-Ds*0.50_8)/beta
    ts2=(r+Ds*0.50_8)/beta
    if (tp1<=t.and.t<=tp2) then
       dom="Fp"
       return
    end if
    if(tp2<=ts1) then
       if (ts1<=t.and.t<=ts2) then
          dom="Fs"
          return
       end if
       if (tp2<=t.and.t<=ts1) then
          dom="I "
          return
       end if
    else
       if (tp2<=t.and.t<=ts2) then
          dom="Fs"
          return
       end if
    end if
    if (ts2<=t) then
       dom="S "
       return
    else
       if (t>=tp1) then
          write(*,*) "Which domain??"
       end if
    end if
  end subroutine domain
end module solve_module
