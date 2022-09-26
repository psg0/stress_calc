module solve_module
  use setting_module
  use kernel_module
  use func_module
  use test_module
  implicit none
  private
  public exact_solve, fdph_solve
contains
  subroutine exact_solve(s11,s12,s13,s21,s22,s23,s31,s32,s33)
    implicit none
    integer :: i,j,m,dm
    real(8), dimension(Mmax,Nmax), intent(inout) :: s11,s12,s13,s21,s22,s23,s31,s32,s33
    real(8) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t
    real(8) :: dd1,dd2,dd3,keisu
    open(12,file="calc_process_exact.dat",status="replace")
    !===initialization===
    call sigma_ini(s11,s12,s13,s21,s22,s23,s31,s32,s33)
    !===time marching=== honest convolution
    do i=1,Nmax !receiver   
       if (mod(i,50)==0) then
          write(*,*) "now:", i
          write(12,*) "now:", i
       end if
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
          do m=1,Mmax
             do dm=1,m-1
                dd1=D1(m-dm,j); dd2=D2(m-dm,j); dd3=D3(m-dm,j)
                t=Dt*dm
                s11(m,i)=s11(m,i)+Kdynf111(t,Dt)*dd1+Kdynf112(t,Dt)*dd2+Kdynf113(t,Dt)*dd3
                s12(m,i)=s12(m,i)+Kdynf121(t,Dt)*dd1+Kdynf122(t,Dt)*dd2+Kdynf123(t,Dt)*dd3
                s13(m,i)=s13(m,i)+Kdynf131(t,Dt)*dd1+Kdynf132(t,Dt)*dd2+Kdynf133(t,Dt)*dd3
                s21(m,i)=s21(m,i)+Kdynf211(t,Dt)*dd1+Kdynf212(t,Dt)*dd2+Kdynf213(t,Dt)*dd3
                s22(m,i)=s22(m,i)+Kdynf221(t,Dt)*dd1+Kdynf222(t,Dt)*dd2+Kdynf223(t,Dt)*dd3
                s23(m,i)=s23(m,i)+Kdynf231(t,Dt)*dd1+Kdynf232(t,Dt)*dd2+Kdynf233(t,Dt)*dd3
                s31(m,i)=s31(m,i)+Kdynf311(t,Dt)*dd1+Kdynf312(t,Dt)*dd2+Kdynf313(t,Dt)*dd3
                s32(m,i)=s32(m,i)+Kdynf321(t,Dt)*dd1+Kdynf322(t,Dt)*dd2+Kdynf323(t,Dt)*dd3
                s33(m,i)=s33(m,i)+Kdynf331(t,Dt)*dd1+Kdynf332(t,Dt)*dd2+Kdynf333(t,Dt)*dd3
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
  end subroutine exact_solve

 subroutine fdph_solve(s11,s12,s13,s21,s22,s23,s31,s32,s33) !o(N^2M)
    implicit none
    integer :: i,j,bl,m,dm
    real(8), dimension(Mmax,Nmax), intent(inout) :: s11,s12,s13,s21,s22,s23,s31,s32,s33
    real(8) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,xsg1,xsg2,xsg3,r,t,t1
    real(8) :: krnFp_lis(27),krnFs_lis(27),htmFp_lis(27),htmFs_lis(27)
    real(8) :: krnInsmp_lis(27)
    real(8) :: ks(27)
    real(8) :: dd1,dd2,dd3,keisu
    integer :: istar_lis(Nmax),istar, block_lis(Nmax) !istar_ind
    integer :: ind,nq,mip,mis,mii,mi,mj,msmin
    integer :: inds
    real(8) :: xr1star,xr2star,xr3star,xry,yry,zry,rry,x,y,z,ri,rj,rtrue
    integer :: cal_was
    
    open(13,file="calc_process_fdph.dat",status="replace")
    ! open(70,file="111.dat",status="replace")
    ! open(71,file="112.dat",status="replace")
    ! open(72,file="113.dat",status="replace")
    ! open(73,file="121.dat",status="replace")
    ! open(74,file="122.dat",status="replace")
    ! open(75,file="123.dat",status="replace")
    ! open(76,file="131.dat",status="replace")
    ! open(77,file="132.dat",status="replace")
    ! open(78,file="133.dat",status="replace")
    ! open(79,file="211.dat",status="replace")
    ! open(80,file="212.dat",status="replace")
    ! open(81,file="213.dat",status="replace")
    ! open(82,file="221.dat",status="replace")
    ! open(83,file="222.dat",status="replace")
    ! open(84,file="223.dat",status="replace")
    ! open(85,file="231.dat",status="replace")
    ! open(86,file="232.dat",status="replace")
    ! open(87,file="233.dat",status="replace")
    ! open(88,file="311.dat",status="replace")
    ! open(89,file="312.dat",status="replace")
    ! open(90,file="313.dat",status="replace")
    ! open(91,file="321.dat",status="replace")
    ! open(92,file="322.dat",status="replace")
    ! open(93,file="323.dat",status="replace")
    ! open(94,file="331.dat",status="replace")
    ! open(95,file="332.dat",status="replace")
    ! open(96,file="333.dat",status="replace")
    open(14,file="k111_s.dat",status="replace")
    !open(70,file="i2.dat",status="replace")
    !===initialization===
    Ds1(:,:)=0.0_8; Ds2(:,:)=0.0_8; Ds3(:,:)=0.0_8
    call sigma_ini(s11,s12,s13,s21,s22,s23,s31,s32,s33)
    call Dhat_ini()
    call calc_istar(nxb,nyb,nzb,istar_lis,block_lis,block_center) !determine representative points of observation
    cal_was=-1
    !write(*,*) block_lis
    call check_same(nxb,nyb,nzb,block_lis,block_center)
    
    !===time marching===
    do m=Mst+1,Mmax
       if (mod(m,10)==0) then
          write(*,*) "now: m=",m
       end if
       !calculate Ds and Dhat
       do j=1,Nmax
          Ds1(m,j)=Ds1(m-1,j)+D1(m-1,j)*Dt
          Ds2(m,j)=Ds2(m-1,j)+D2(m-1,j)*Dt
          Ds3(m,j)=Ds3(m-1,j)+D3(m-1,j)*Dt
          xsa1=xa(j); xsa2=ya(j); xsb1=xb(j); xsb2=yb(j)
          xsc1=xc(j); xsc2=yc(j); xs3=za(j)
          xsg1=xg(j); xsg2=yg(j); xsg3=zg(j)
          do i=1,Nmax
             bl=block_lis(i)
             istar=block_center(bl)
             xr1star=xg(istar); xr2star=yg(istar); xr3star=zg(istar)
             xr1=xg(i); xr2=yg(i); xr3=zg(i)
             !call Kdyns_all(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
             !rtrue=sqrt((xr1-xsg1)**2+(xr2-xsg2)**2+(xr3-xsg3)**2)
             !call krnFp_all(i,j,rtrue,krnFp_lis)
             !call krnFs_all(i,j,rtrue,krnFs_lis)

             ! To calculate h (call htmFp & htmFs)
             if (cal_was/=bl) then
                cal_was=bl
                call Kdyns_all(xr1star,xr2star,xr3star,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
                r=sqrt((xr1star-xsg1)**2+(xr2star-xsg2)**2+(xr3star-xsg3)**2)
                call krnFp_all(istar,j,r,krnFp_lis)
                call krnFs_all(istar,j,r,krnFs_lis)
             end if
             
             do dm=0,mpmax
                ind=m-ceiling(((r-Ds*0.50_8)/alpha+dm*Dt)/Dt)
                if (ind<1 .or. ind>Mmax) then
                   cycle
                end if
                call htmFp(istar,j,dm,krnFp_lis,htmFp_lis)
                dd1=D1(ind,j); dd2=D2(ind,j); dd3=D3(ind,j)
                D111hatFp(m,bl,j)=D111hatFp(m,bl,j)+dd1*htmFp_lis(1)
                D112hatFp(m,bl,j)=D112hatFp(m,bl,j)+dd2*htmFp_lis(2)
                D113hatFp(m,bl,j)=D113hatFp(m,bl,j)+dd3*htmFp_lis(3)
                D121hatFp(m,bl,j)=D121hatFp(m,bl,j)+dd1*htmFp_lis(4)
                D122hatFp(m,bl,j)=D122hatFp(m,bl,j)+dd2*htmFp_lis(5)
                D123hatFp(m,bl,j)=D123hatFp(m,bl,j)+dd3*htmFp_lis(6)
                D131hatFp(m,bl,j)=D131hatFp(m,bl,j)+dd1*htmFp_lis(7)
                D132hatFp(m,bl,j)=D132hatFp(m,bl,j)+dd2*htmFp_lis(8)
                D133hatFp(m,bl,j)=D133hatFp(m,bl,j)+dd3*htmFp_lis(9)
                D211hatFp(m,bl,j)=D211hatFp(m,bl,j)+dd1*htmFp_lis(10)
                D212hatFp(m,bl,j)=D212hatFp(m,bl,j)+dd2*htmFp_lis(11)
                D213hatFp(m,bl,j)=D213hatFp(m,bl,j)+dd3*htmFp_lis(12)
                D221hatFp(m,bl,j)=D221hatFp(m,bl,j)+dd1*htmFp_lis(13)
                D222hatFp(m,bl,j)=D222hatFp(m,bl,j)+dd2*htmFp_lis(14)
                D223hatFp(m,bl,j)=D223hatFp(m,bl,j)+dd3*htmFp_lis(15)
                D231hatFp(m,bl,j)=D231hatFp(m,bl,j)+dd1*htmFp_lis(16)
                D232hatFp(m,bl,j)=D232hatFp(m,bl,j)+dd2*htmFp_lis(17)
                D233hatFp(m,bl,j)=D233hatFp(m,bl,j)+dd3*htmFp_lis(18)
                D311hatFp(m,bl,j)=D311hatFp(m,bl,j)+dd1*htmFp_lis(19)
                D312hatFp(m,bl,j)=D312hatFp(m,bl,j)+dd2*htmFp_lis(20)
                D313hatFp(m,bl,j)=D313hatFp(m,bl,j)+dd3*htmFp_lis(21)
                D321hatFp(m,bl,j)=D321hatFp(m,bl,j)+dd1*htmFp_lis(22)
                D322hatFp(m,bl,j)=D322hatFp(m,bl,j)+dd2*htmFp_lis(23)
                D323hatFp(m,bl,j)=D323hatFp(m,bl,j)+dd3*htmFp_lis(24)
                D331hatFp(m,bl,j)=D331hatFp(m,bl,j)+dd1*htmFp_lis(25)
                D332hatFp(m,bl,j)=D332hatFp(m,bl,j)+dd2*htmFp_lis(26)
                D333hatFp(m,bl,j)=D333hatFp(m,bl,j)+dd3*htmFp_lis(27)
             end do

             msmin=0
             t1=r/alpha+mpmax*Dt-Ds/alpha*0.50_8
             if(t1>=r/beta-Ds/beta*0.50_8) then
                msmin=floor((t1-(r/beta-Ds/beta*0.50_8))/Dt)
                msmin=msmin+1 !There should not be any overlap between domain Fp and Fs
             end if
             do dm=msmin,msmax-1
                ind=m-ceiling(((r-Ds*0.50_8)/beta+dm*Dt)/Dt)
                if (ind<1 .or. ind>Mmax) then
                   cycle
                end if
                call htmFs(istar,j,dm,krnFs_lis,htmFs_lis)
                dd1=D1(ind,j); dd2=D2(ind,j); dd3=D3(ind,j)
                D111hatFs(m,bl,j)=D111hatFs(m,bl,j)+dd1*htmFs_lis(1)
                D112hatFs(m,bl,j)=D112hatFs(m,bl,j)+dd2*htmFs_lis(2)
                D113hatFs(m,bl,j)=D113hatFs(m,bl,j)+dd3*htmFs_lis(3)
                D121hatFs(m,bl,j)=D121hatFs(m,bl,j)+dd1*htmFs_lis(4)
                D122hatFs(m,bl,j)=D122hatFs(m,bl,j)+dd2*htmFs_lis(5)
                D123hatFs(m,bl,j)=D123hatFs(m,bl,j)+dd3*htmFs_lis(6)
                D131hatFs(m,bl,j)=D131hatFs(m,bl,j)+dd1*htmFs_lis(7)
                D132hatFs(m,bl,j)=D132hatFs(m,bl,j)+dd2*htmFs_lis(8)
                D133hatFs(m,bl,j)=D133hatFs(m,bl,j)+dd3*htmFs_lis(9)
                D211hatFs(m,bl,j)=D211hatFs(m,bl,j)+dd1*htmFs_lis(10)
                D212hatFs(m,bl,j)=D212hatFs(m,bl,j)+dd2*htmFs_lis(11)
                D213hatFs(m,bl,j)=D213hatFs(m,bl,j)+dd3*htmFs_lis(12)
                D221hatFs(m,bl,j)=D221hatFs(m,bl,j)+dd1*htmFs_lis(13)
                D222hatFs(m,bl,j)=D222hatFs(m,bl,j)+dd2*htmFs_lis(14)
                D223hatFs(m,bl,j)=D223hatFs(m,bl,j)+dd3*htmFs_lis(15)
                D231hatFs(m,bl,j)=D231hatFs(m,bl,j)+dd1*htmFs_lis(16)
                D232hatFs(m,bl,j)=D232hatFs(m,bl,j)+dd2*htmFs_lis(17)
                D233hatFs(m,bl,j)=D233hatFs(m,bl,j)+dd3*htmFs_lis(18)
                D311hatFs(m,bl,j)=D311hatFs(m,bl,j)+dd1*htmFs_lis(19)
                D312hatFs(m,bl,j)=D312hatFs(m,bl,j)+dd2*htmFs_lis(20)
                D313hatFs(m,bl,j)=D313hatFs(m,bl,j)+dd3*htmFs_lis(21)
                D321hatFs(m,bl,j)=D321hatFs(m,bl,j)+dd1*htmFs_lis(22)
                D322hatFs(m,bl,j)=D322hatFs(m,bl,j)+dd2*htmFs_lis(23)
                D323hatFs(m,bl,j)=D323hatFs(m,bl,j)+dd3*htmFs_lis(24)
                D331hatFs(m,bl,j)=D331hatFs(m,bl,j)+dd1*htmFs_lis(25)
                D332hatFs(m,bl,j)=D332hatFs(m,bl,j)+dd2*htmFs_lis(26)
                D333hatFs(m,bl,j)=D333hatFs(m,bl,j)+dd3*htmFs_lis(27)
             end do
          end do
       end do

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
             !computation of mi
             call calc_mi("2",istar,i,j,mip,mis,mii,ri,rj)
     
             !Amplitude 
             rtrue=sqrt((xr1-xsg1)**2+(xr2-xsg2)**2+(xr3-xsg3)**2)
             call Kdyns_all(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
             call krnFp_all(i,j,rtrue,krnFp_lis)
             call krnFs_all(i,j,rtrue,krnFs_lis)
             !domain Fp
             !mj=max(ceiling(rj/(alpha*Dt)-Ds/(alpha*Dt)*0.5),1)
             ind=m+mip
             !if (mj<=m-1.and.1<=ind.and.ind<=Mmax) then
             !write(*,*) m,mip,ind,Mmax
             if (1<=ind.and.ind<=Mmax) then
                s11(ind,i)=s11(ind,i)+D111hatFp(m,bl,j)*krnFp_lis(1)&
                     & +D112hatFp(m,bl,j)*krnFp_lis(2)+D113hatFp(m,bl,j)*krnFp_lis(3)
                s12(ind,i)=s12(ind,i)+D121hatFp(m,bl,j)*krnFp_lis(4)&
                     & +D122hatFp(m,bl,j)*krnFp_lis(5)+D123hatFp(m,bl,j)*krnFp_lis(6)
                s13(ind,i)=s13(ind,i)+D131hatFp(m,bl,j)*krnFp_lis(7)&
                     & +D132hatFp(m,bl,j)*krnFp_lis(8)+D133hatFp(m,bl,j)*krnFp_lis(9)
                s21(ind,i)=s21(ind,i)+D211hatFp(m,bl,j)*krnFp_lis(10)&
                     & +D212hatFp(m,bl,j)*krnFp_lis(11)+D213hatFp(m,bl,j)*krnFp_lis(12)
                s22(ind,i)=s22(ind,i)+D221hatFp(m,bl,j)*krnFp_lis(13)&
                     & +D222hatFp(m,bl,j)*krnFp_lis(14)+D223hatFp(m,bl,j)*krnFp_lis(15)
                s23(ind,i)=s23(ind,i)+D231hatFp(m,bl,j)*krnFp_lis(16)&
                     & +D232hatFp(m,bl,j)*krnFp_lis(17)+D233hatFp(m,bl,j)*krnFp_lis(18)
                s31(ind,i)=s31(ind,i)+D311hatFp(m,bl,j)*krnFp_lis(19)&
                     & +D312hatFp(m,bl,j)*krnFp_lis(20)+D313hatFp(m,bl,j)*krnFp_lis(21)
                s32(ind,i)=s32(ind,i)+D321hatFp(m,bl,j)*krnFp_lis(22)&
                     & +D322hatFp(m,bl,j)*krnFp_lis(23)+D323hatFp(m,bl,j)*krnFp_lis(24)
                s33(ind,i)=s33(ind,i)+D331hatFp(m,bl,j)*krnFp_lis(25)&
                     & +D332hatFp(m,bl,j)*krnFp_lis(26)+D333hatFp(m,bl,j)*krnFp_lis(27)
             end if
             !domain Fs
             !mj=max(ceiling(rj/(beta*Dt)-Ds/(beta*Dt)*0.5),1)
             ind=m+mis
             !if (mj<=m-1.and.1<=ind.and.ind<=Mmax) then
             if (1<=ind.and.ind<=Mmax) then
                s11(ind,i)=s11(ind,i)+D111hatFs(m,bl,j)*krnFs_lis(1)&
                     & +D112hatFs(m,bl,j)*krnFs_lis(2)+D113hatFs(m,bl,j)*krnFs_lis(3)
                s12(ind,i)=s12(ind,i)+D121hatFs(m,bl,j)*krnFs_lis(4)&
                     & +D122hatFs(m,bl,j)*krnFs_lis(5)+D123hatFs(m,bl,j)*krnFs_lis(6)
                s13(ind,i)=s13(ind,i)+D131hatFs(m,bl,j)*krnFs_lis(7)&
                     & +D132hatFs(m,bl,j)*krnFs_lis(8)+D133hatFs(m,bl,j)*krnFs_lis(9)
                s21(ind,i)=s21(ind,i)+D211hatFs(m,bl,j)*krnFs_lis(10)&
                     & +D212hatFs(m,bl,j)*krnFs_lis(11)+D213hatFs(m,bl,j)*krnFs_lis(12)
                s22(ind,i)=s22(ind,i)+D221hatFs(m,bl,j)*krnFs_lis(13)&
                     & +D222hatFs(m,bl,j)*krnFs_lis(14)+D223hatFs(m,bl,j)*krnFs_lis(15)
                s23(ind,i)=s23(ind,i)+D231hatFs(m,bl,j)*krnFs_lis(16)&
                     & +D232hatFs(m,bl,j)*krnFs_lis(17)+D233hatFs(m,bl,j)*krnFs_lis(18)    
                s31(ind,i)=s31(ind,i)+D311hatFs(m,bl,j)*krnFs_lis(19)&
                     & +D312hatFs(m,bl,j)*krnFs_lis(20)+D313hatFs(m,bl,j)*krnFs_lis(21)
                s32(ind,i)=s32(ind,i)+D321hatFs(m,bl,j)*krnFs_lis(22)&
                     & +D322hatFs(m,bl,j)*krnFs_lis(23)+D323hatFs(m,bl,j)*krnFs_lis(24)
                s33(ind,i)=s33(ind,i)+D331hatFs(m,bl,j)*krnFs_lis(25)&
                     & +D332hatFs(m,bl,j)*krnFs_lis(26)+D333hatFs(m,bl,j)*krnFs_lis(27)
             end if
             
             !domain I
             ind=m+mip
             !ind=m+mii
             if ((rj-Ds*0.50_8)/beta>(rj+Ds*0.50_8)/alpha) then
                do nq=0,Nqmax
                   dm=ceiling(((rj+Ds*0.50_8)/alpha+((rj-Ds*0.50_8)/beta-(rj+Ds*0.50_8)/alpha)*nq/(Nqmax*1.0_8))/Dt)
                   call krnInsmp(rj,nq,krnInsmp_lis)
                   if (m-dm>=1.and.m-dm<=m-1.and.1<=ind.and.ind<=Mmax) then
                      dd1=Ds1(m-dm,j); dd2=Ds2(m-dm,j); dd3=Ds3(m-dm,j)
                      sigma11(ind,i)=sigma11(ind,i)+dd1*krnInsmp_lis(1)+dd2*krnInsmp_lis(2)+dd3*krnInsmp_lis(3)
                      sigma12(ind,i)=sigma12(ind,i)+dd1*krnInsmp_lis(4)+dd2*krnInsmp_lis(5)+dd3*krnInsmp_lis(6)
                      sigma13(ind,i)=sigma13(ind,i)+dd1*krnInsmp_lis(7)+dd2*krnInsmp_lis(8)+dd3*krnInsmp_lis(9)
                      sigma21(ind,i)=sigma21(ind,i)+dd1*krnInsmp_lis(10)+dd2*krnInsmp_lis(11)+dd3*krnInsmp_lis(12)
                      sigma22(ind,i)=sigma22(ind,i)+dd1*krnInsmp_lis(13)+dd2*krnInsmp_lis(14)+dd3*krnInsmp_lis(15)
                      sigma23(ind,i)=sigma23(ind,i)+dd1*krnInsmp_lis(16)+dd2*krnInsmp_lis(17)+dd3*krnInsmp_lis(18)
                      sigma31(ind,i)=sigma31(ind,i)+dd1*krnInsmp_lis(19)+dd2*krnInsmp_lis(20)+dd3*krnInsmp_lis(21)
                      sigma32(ind,i)=sigma32(ind,i)+dd1*krnInsmp_lis(22)+dd2*krnInsmp_lis(23)+dd3*krnInsmp_lis(24)
6                     sigma33(ind,i)=sigma33(ind,i)+dd1*krnInsmp_lis(25)+dd2*krnInsmp_lis(26)+dd3*krnInsmp_lis(27)
                   end if
                end do
             end if
             
             !domain S
             inds=ceiling((rj+Ds*0.50_8)/(beta*Dt))
             ind=m+mis
             if (1<=ind.and.ind<=Mmax) then
                call calc_ks(inds*Dt,ks) !Kdyns_all has already calculated!
                do dm=inds,m-Mst
                   t=dm*Dt
                   dd1=D1(m-dm,j); dd2=D2(m-dm,j); dd3=D3(m-dm,j)
                   
                   sigma11(ind,i)=sigma11(ind,i)+dd1*ks(1)+dd2*ks(2)+dd3*ks(3)
                   sigma12(ind,i)=sigma12(ind,i)+dd1*ks(4)+dd2*ks(5)+dd3*ks(6)
                   sigma13(ind,i)=sigma13(ind,i)+dd1*ks(7)+dd2*ks(8)+dd3*ks(9)
                   sigma21(ind,i)=sigma21(ind,i)+dd1*ks(10)+dd2*ks(11)+dd3*ks(12)
                   sigma22(ind,i)=sigma22(ind,i)+dd1*ks(13)+dd2*ks(14)+dd3*ks(15)
                   sigma23(ind,i)=sigma23(ind,i)+dd1*ks(16)+dd2*ks(17)+dd3*ks(18)
                   sigma31(ind,i)=sigma31(ind,i)+dd1*ks(19)+dd2*ks(20)+dd3*ks(21)
                   sigma32(ind,i)=sigma32(ind,i)+dd1*ks(22)+dd2*ks(23)+dd3*ks(24)
                   sigma33(ind,i)=sigma33(ind,i)+dd1*ks(25)+dd2*ks(26)+dd3*ks(27)
                end do
             end if
             !write(14,*)
           end do
           !write(14,*)
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
    close(13)
    close(14)
    !close(70)
  end subroutine fdph_solve

  ! domain Fp
  subroutine krnFp_all(i,j,r,krnFp_lis) !Indeed, i&j is not needed
    implicit none
    integer, intent(in) :: i,j
    real(8), intent(in) :: r
    real(8), dimension(27), intent(inout) :: krnFp_lis
    integer :: l,dm
    real(8) :: t

    krnFp_lis(:)=0.0_8
    !write(*,*) "mpmax:",mpmax
    do dm=0,mpmax
       t=r/alpha+dm*Dt-Ds/alpha*0.50_8 !t can be negative, but kernel become 0 when t<0
       krnFp_lis(1)=krnFp_lis(1)+Kdynf111(t,Dt)*Dt
       krnFp_lis(2)=krnFp_lis(2)+Kdynf112(t,Dt)*Dt
       krnFp_lis(3)=krnFp_lis(3)+Kdynf113(t,Dt)*Dt
       krnFp_lis(4)=krnFp_lis(4)+Kdynf121(t,Dt)*Dt
       krnFp_lis(5)=krnFp_lis(5)+Kdynf122(t,Dt)*Dt
       krnFp_lis(6)=krnFp_lis(6)+Kdynf123(t,Dt)*Dt
       krnFp_lis(7)=krnFp_lis(7)+Kdynf131(t,Dt)*Dt
       krnFp_lis(8)=krnFp_lis(8)+Kdynf132(t,Dt)*Dt
       krnFp_lis(9)=krnFp_lis(9)+Kdynf133(t,Dt)*Dt
       krnFp_lis(10)=krnFp_lis(10)+Kdynf211(t,Dt)*Dt
       krnFp_lis(11)=krnFp_lis(11)+Kdynf212(t,Dt)*Dt
       krnFp_lis(12)=krnFp_lis(12)+Kdynf213(t,Dt)*Dt
       krnFp_lis(13)=krnFp_lis(13)+Kdynf221(t,Dt)*Dt
       krnFp_lis(14)=krnFp_lis(14)+Kdynf222(t,Dt)*Dt
       krnFp_lis(15)=krnFp_lis(15)+Kdynf223(t,Dt)*Dt
       krnFp_lis(16)=krnFp_lis(16)+Kdynf231(t,Dt)*Dt
       krnFp_lis(17)=krnFp_lis(17)+Kdynf232(t,Dt)*Dt
       krnFp_lis(18)=krnFp_lis(18)+Kdynf233(t,Dt)*Dt
       krnFp_lis(19)=krnFp_lis(19)+Kdynf311(t,Dt)*Dt
       krnFp_lis(20)=krnFp_lis(20)+Kdynf312(t,Dt)*Dt
       krnFp_lis(21)=krnFp_lis(21)+Kdynf313(t,Dt)*Dt
       krnFp_lis(22)=krnFp_lis(22)+Kdynf321(t,Dt)*Dt
       krnFp_lis(23)=krnFp_lis(23)+Kdynf322(t,Dt)*Dt
       krnFp_lis(24)=krnFp_lis(24)+Kdynf323(t,Dt)*Dt
       krnFp_lis(25)=krnFp_lis(25)+Kdynf331(t,Dt)*Dt
       krnFp_lis(26)=krnFp_lis(26)+Kdynf332(t,Dt)*Dt
       krnFp_lis(27)=krnFp_lis(27)+Kdynf333(t,Dt)*Dt
    end do
  end subroutine krnFp_all

  ! domain Fs
  subroutine krnFs_all(i,j,r,krnFs_lis)
    implicit none
    integer, intent(in) :: i,j
    real(8), intent(in) :: r
    real(8), dimension(27), intent(inout) :: krnFs_lis
    integer :: l,dm, msmin
    real(8) :: t,t1,t2
    
    msmin=0
    t1=r/alpha+mpmax*Dt-Ds/alpha*0.50_8
    if(t1>=r/beta-Ds/beta*0.50_8) then
       msmin=floor((t1-(r/beta-Ds/beta*0.50_8))/Dt)
       msmin=msmin+1 !There should not be any overlap between domain Fp and Fs
    end if
    !t2=r/beta+msmin*Dt-Ds/beta*0.50_8 !t2>t1
    krnFs_lis(:)=0.0_8
    do dm=msmin,msmax-1,1
       t=r/beta+dm*Dt-Ds/beta*0.50_8
       krnFs_lis(1)=krnFs_lis(1)+Kdynf111(t,Dt)*Dt
       krnFs_lis(2)=krnFs_lis(2)+Kdynf112(t,Dt)*Dt
       krnFs_lis(3)=krnFs_lis(3)+Kdynf113(t,Dt)*Dt
       krnFs_lis(4)=krnFs_lis(4)+Kdynf121(t,Dt)*Dt
       krnFs_lis(5)=krnFs_lis(5)+Kdynf122(t,Dt)*Dt
       krnFs_lis(6)=krnFs_lis(6)+Kdynf123(t,Dt)*Dt
       krnFs_lis(7)=krnFs_lis(7)+Kdynf131(t,Dt)*Dt
       krnFs_lis(8)=krnFs_lis(8)+Kdynf132(t,Dt)*Dt
       krnFs_lis(9)=krnFs_lis(9)+Kdynf133(t,Dt)*Dt
       krnFs_lis(10)=krnFs_lis(10)+Kdynf211(t,Dt)*Dt
       krnFs_lis(11)=krnFs_lis(11)+Kdynf212(t,Dt)*Dt
       krnFs_lis(12)=krnFs_lis(12)+Kdynf213(t,Dt)*Dt
       krnFs_lis(13)=krnFs_lis(13)+Kdynf221(t,Dt)*Dt
       krnFs_lis(14)=krnFs_lis(14)+Kdynf222(t,Dt)*Dt
       krnFs_lis(15)=krnFs_lis(15)+Kdynf223(t,Dt)*Dt
       krnFs_lis(16)=krnFs_lis(16)+Kdynf231(t,Dt)*Dt
       krnFs_lis(17)=krnFs_lis(17)+Kdynf232(t,Dt)*Dt
       krnFs_lis(18)=krnFs_lis(18)+Kdynf233(t,Dt)*Dt
       krnFs_lis(19)=krnFs_lis(19)+Kdynf311(t,Dt)*Dt
       krnFs_lis(20)=krnFs_lis(20)+Kdynf312(t,Dt)*Dt
       krnFs_lis(21)=krnFs_lis(21)+Kdynf313(t,Dt)*Dt
       krnFs_lis(22)=krnFs_lis(22)+Kdynf321(t,Dt)*Dt
       krnFs_lis(23)=krnFs_lis(23)+Kdynf322(t,Dt)*Dt
       krnFs_lis(24)=krnFs_lis(24)+Kdynf323(t,Dt)*Dt
       krnFs_lis(25)=krnFs_lis(25)+Kdynf331(t,Dt)*Dt
       krnFs_lis(26)=krnFs_lis(26)+Kdynf332(t,Dt)*Dt
       krnFs_lis(27)=krnFs_lis(27)+Kdynf333(t,Dt)*Dt
    end do
  end subroutine krnFs_all

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
       if (abs(krnFs_lis(k))<eps) then
          htmFs_lis(k)=0.0_8
       else
          htmFs_lis(k)=val(k)/krnFs_lis(k)
       end if
    end do  
  end subroutine htmFs

  subroutine krnInsmp(rj,nq,krnInsmp_lis)
    implicit none
    integer, intent(in) :: nq
    real(8), intent(in) :: rj
    real(8), dimension(27), intent(inout) :: krnInsmp_lis
    real(8) :: dtofi,tnq0,tnq1,tnq2,tcntq0,tcntq1,kkq0(27),kkq1(27)
    real(8) :: ts,te
    integer :: k

    ts=(rj+Ds*0.50_8)/alpha
    te=(rj-Ds*0.50_8)/beta
    dtofi=te-ts !time range of domain I
    
    !nq from 0 to Nqmax
    tnq0=ts+dtofi*(nq-1)/(Nqmax*1.0_8)
    tnq1=ts+dtofi*(nq)/(Nqmax*1.0_8)
    tnq2=ts+dtofi*(nq+1)/(Nqmax*1.0_8)
    tcntq0=(tnq0+tnq1)*0.50_8
    tcntq1=(tnq1+tnq2)*0.50_8
    kkq0(:)=0.0_8
    kkq1(:)=0.0_8
    if(ts<tcntq0.and.tcntq0<te) then
       kkq0(1)=Kdynf111(tcntq0,Dt)
       kkq0(2)=Kdynf112(tcntq0,Dt)
       kkq0(3)=Kdynf113(tcntq0,Dt)
       kkq0(4)=Kdynf121(tcntq0,Dt)
       kkq0(5)=Kdynf122(tcntq0,Dt)
       kkq0(6)=Kdynf123(tcntq0,Dt)
       kkq0(7)=Kdynf131(tcntq0,Dt)
       kkq0(8)=Kdynf132(tcntq0,Dt)
       kkq0(9)=Kdynf133(tcntq0,Dt)      
       kkq0(10)=Kdynf211(tcntq0,Dt)
       kkq0(11)=Kdynf212(tcntq0,Dt)
       kkq0(12)=Kdynf213(tcntq0,Dt)
       kkq0(13)=Kdynf221(tcntq0,Dt)
       kkq0(14)=Kdynf222(tcntq0,Dt)
       kkq0(15)=Kdynf223(tcntq0,Dt)
       kkq0(16)=Kdynf231(tcntq0,Dt)
       kkq0(17)=Kdynf232(tcntq0,Dt)
       kkq0(18)=Kdynf233(tcntq0,Dt)
       kkq0(19)=Kdynf311(tcntq0,Dt)
       kkq0(20)=Kdynf312(tcntq0,Dt)
       kkq0(21)=Kdynf313(tcntq0,Dt)
       kkq0(22)=Kdynf321(tcntq0,Dt)
       kkq0(23)=Kdynf322(tcntq0,Dt)
       kkq0(24)=Kdynf323(tcntq0,Dt)
       kkq0(25)=Kdynf331(tcntq0,Dt)
       kkq0(26)=Kdynf332(tcntq0,Dt)
       kkq0(27)=Kdynf333(tcntq0,Dt)    
    end if
    if(ts<tcntq1.and.tcntq1<te) then
       kkq1(1)=Kdynf111(tcntq1,Dt)
       kkq1(2)=Kdynf112(tcntq1,Dt)
       kkq1(3)=Kdynf113(tcntq1,Dt)
       kkq1(4)=Kdynf121(tcntq1,Dt)
       kkq1(5)=Kdynf122(tcntq1,Dt)
       kkq1(6)=Kdynf123(tcntq1,Dt)
       kkq1(7)=Kdynf131(tcntq1,Dt)
       kkq1(8)=Kdynf132(tcntq1,Dt)
       kkq1(9)=Kdynf133(tcntq1,Dt)      
       kkq1(10)=Kdynf211(tcntq1,Dt)
       kkq1(11)=Kdynf212(tcntq1,Dt)
       kkq1(12)=Kdynf213(tcntq1,Dt)
       kkq1(13)=Kdynf221(tcntq1,Dt)
       kkq1(14)=Kdynf222(tcntq1,Dt)
       kkq1(15)=Kdynf223(tcntq1,Dt)
       kkq1(16)=Kdynf231(tcntq1,Dt)
       kkq1(17)=Kdynf232(tcntq1,Dt)
       kkq1(18)=Kdynf233(tcntq1,Dt)
       kkq1(19)=Kdynf311(tcntq1,Dt)
       kkq1(20)=Kdynf312(tcntq1,Dt)
       kkq1(21)=Kdynf313(tcntq1,Dt)
       kkq1(22)=Kdynf321(tcntq1,Dt)
       kkq1(23)=Kdynf322(tcntq1,Dt)
       kkq1(24)=Kdynf323(tcntq1,Dt)
       kkq1(25)=Kdynf331(tcntq1,Dt)
       kkq1(26)=Kdynf332(tcntq1,Dt)
       kkq1(27)=Kdynf333(tcntq1,Dt)    
    end if
    do k=1,27
       krnInsmp_lis(k)=kkq0(k)-kkq1(k)
    end do
  end subroutine krnInsmp  

end module solve_module
