module test_module
  use setting_module
  use kernel_module
  implicit none
  private
  public fig6_7,output,dattovtk
contains
  !(subroutine fig6_7: To check that kernels are correct!)
  !please confirm that alpha=sqrt(3.0_8); beta=1.0_8
  subroutine fig6_7()
    implicit none
    integer(8) :: i,j
    integer(8), parameter :: n=2000
    real(8) :: xr(3)=(/0.6,0.3,0.2/)
    real(8) :: xsa(2),xsb(2),xsc(2),xs3,t,dtime
    real(8) :: L6(1:10,0:n,2), L7(1:10,0:n,2) !L111(1)~L333(10)
    real(8), parameter :: t1=0.50_8
    open(10,file='fig6/L111.dat', status='replace')!No1
    open(11,file='fig6/L121.dat', status='replace')!No2
    open(12,file='fig6/L221.dat', status='replace')!No3
    open(13,file='fig6/L311.dat', status='replace')!No4
    open(14,file='fig6/L321.dat', status='replace')!No5
    open(15,file='fig6/L331.dat', status='replace')!No6
    open(16,file='fig6/L113.dat', status='replace')!No7
    open(17,file='fig6/L123.dat', status='replace')!No8
    open(18,file='fig6/L313.dat', status='replace')!No9
    open(19,file='fig6/L333.dat', status='replace')!No10
    open(20,file='fig7/L111.dat', status='replace')!No1
    open(21,file='fig7/L121.dat', status='replace')!No2
    open(22,file='fig7/L221.dat', status='replace')!No3
    open(23,file='fig7/L311.dat', status='replace')!No4
    open(24,file='fig7/L321.dat', status='replace')!No5
    open(25,file='fig7/L331.dat', status='replace')!No6
    open(26,file='fig7/L113.dat', status='replace')!No7
    open(27,file='fig7/L123.dat', status='replace')!No8
    open(28,file='fig7/L313.dat', status='replace')!No9
    open(29,file='fig7/L333.dat', status='replace')!No10
    write(*,*) "test(fig6&fig7)"
    write(*,*) "alpha:",alpha
    write(*,*) "beta:",beta
    dtime=1.0_8/n
    !===ABC===
    xsa(1)=0.0_8; xsa(2)=0.0_8; xsb(1)=1.0_8; xsb(2)=0.0_8
    xsc(1)=1.0_8; xsc(2)=1.0_8; xs3=0.0_8
    call Kdyns_all(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    !No1 (a)L111 (ABC)
    !call Kdyns111(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(1,i,1)=Kdynf111(0.0_8,t)     !(t,t)
       L7(1,i,1)=Kdynf111(t-t1,t1) !(t,t1)
    end do
    !No2 (b)L121 (ABC)
    !call Kdyns121(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(2,i,1)=Kdynf121(0.0_8,t)
       L7(2,i,1)=Kdynf121(t-t1,t1)
    end do
    !No3 (c)L221 (ABC)
    !call Kdyns221(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(3,i,1)=Kdynf221(0.0_8,t)
       L7(3,i,1)=Kdynf221(t-t1,t1)
    end do
    !No4 (d)L311 (ABC)
    !call Kdyns311(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(4,i,1)=Kdynf311(0.0_8,t)
       L7(4,i,1)=Kdynf311(t-t1,t1)
    end do
    !No5 (e)L321 (ABC)
    !call Kdyns321(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(5,i,1)=Kdynf321(0.0_8,t)
       L7(5,i,1)=Kdynf321(t-t1,t1)
    end do
    !No6 (f)L331 (ABC)
    !call Kdyns331(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(6,i,1)=Kdynf331(0.0_8,t)
       L7(6,i,1)=Kdynf331(t-t1,t1)
    end do
    !No7 (g)L113 (ABC)
    !call Kdyns113(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(7,i,1)=Kdynf113(0.0_8,t)
       L7(7,i,1)=Kdynf113(t-t1,t1)
    end do
    !No8 (h)L123 (ABC)
    !call Kdyns123(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(8,i,1)=Kdynf123(0.0_8,t)
       L7(8,i,1)=Kdynf123(t-t1,t1)
    end do
    !No9 (i)L313 (ABC)
    !call Kdyns313(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(9,i,1)=Kdynf313(0.0_8,t)
       L7(9,i,1)=Kdynf313(t-t1,t1)
    end do
    !No10 (j)L333 (ABC)
    !call Kdyns333(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(10,i,1)=Kdynf333(0.0_8,t)
       L7(10,i,1)=Kdynf333(t-t1,t1)
    end do

    !===CDA===
    xsa(1)=1.0_8; xsa(2)=1.0_8; xsb(1)=0.0_8; xsb(2)=1.0_8
    xsc(1)=0.0_8; xsc(2)=0.0_8; xs3=0.0_8
    call Kdyns_all(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    !No1 (a)L111 (CDA)
    !call Kdyns111(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(1,i,2)=Kdynf111(0.0_8,t)
       L7(1,i,2)=Kdynf111(t-t1,t1)
    end do
    !No2 (b)L121 (CDA)
    !call Kdyns121(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(2,i,2)=Kdynf121(0.0_8,t)
       L7(2,i,2)=Kdynf121(t-t1,t1)
    end do
    !No3 (c)L221 (CDA)
    !call Kdyns221(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(3,i,2)=Kdynf221(0.0_8,t)
       L7(3,i,2)=Kdynf221(t-t1,t1)
    end do
    !No4 (d)L311 (CDA)
    !call Kdyns311(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(4,i,2)=Kdynf311(0.0_8,t)
       L7(4,i,2)=Kdynf311(t-t1,t1)
    end do
    !No5 (e)L321 (CDA)
    !call Kdyns321(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(5,i,2)=Kdynf321(0.0_8,t)
       L7(5,i,2)=Kdynf321(t-t1,t1)
    end do
    !No6 (f)L331 (CDA)
    !call Kdyns331(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(6,i,2)=Kdynf331(0.0_8,t)
       L7(6,i,2)=Kdynf331(t-t1,t1)
    end do
    !No7 (g)L113 (CDA)
    !call Kdyns113(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(7,i,2)=Kdynf113(0.0_8,t)
       L7(7,i,2)=Kdynf113(t-t1,t1)
    end do
    !No8 (c)L123 (CDA)
    !call Kdyns123(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(8,i,2)=Kdynf123(0.0_8,t)
       L7(8,i,2)=Kdynf123(t-t1,t1)
    end do
    !No9 (d)L313 (CDA)
    !call Kdyns313(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(9,i,2)=Kdynf313(0.0_8,t)
       L7(9,i,2)=Kdynf313(t-t1,t1)
    end do
    !No10 (e)L333 (CDA)
    !call Kdyns333(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(10,i,2)=Kdynf333(0.0_8,t)
       L7(10,i,2)=Kdynf333(t-t1,t1)
    end do
  
    !===to data file===
    do i=0,n
       t=dtime*i
       write(10,*) t,L6(1,i,1),L6(1,i,2),L6(1,i,1)+L6(1,i,2)
       write(11,*) t,L6(2,i,1),L6(2,i,2),L6(2,i,1)+L6(2,i,2)
       write(12,*) t,L6(3,i,1),L6(3,i,2),L6(3,i,1)+L6(3,i,2)
       write(13,*) t,L6(4,i,1),L6(4,i,2),L6(4,i,1)+L6(4,i,2)
       write(14,*) t,L6(5,i,1),L6(5,i,2),L6(5,i,1)+L6(5,i,2)
       write(15,*) t,L6(6,i,1),L6(6,i,2),L6(6,i,1)+L6(6,i,2)
       write(16,*) t,L6(7,i,1),L6(7,i,2),L6(7,i,1)+L6(7,i,2)
       write(17,*) t,L6(8,i,1),L6(8,i,2),L6(8,i,1)+L6(8,i,2)
       write(18,*) t,L6(9,i,1),L6(9,i,2),L6(9,i,1)+L6(9,i,2)
       write(19,*) t,L6(10,i,1),L6(10,i,2),L6(10,i,1)+L6(10,i,2)
       write(20,*) t,L7(1,i,1),L7(1,i,2),L7(1,i,1)+L7(1,i,2)
       write(21,*) t,L7(2,i,1),L7(2,i,2),L7(2,i,1)+L7(2,i,2)
       write(22,*) t,L7(3,i,1),L7(3,i,2),L7(3,i,1)+L7(3,i,2)
       write(23,*) t,L7(4,i,1),L7(4,i,2),L7(4,i,1)+L7(4,i,2)
       write(24,*) t,L7(5,i,1),L7(5,i,2),L7(5,i,1)+L7(5,i,2)
       write(25,*) t,L7(6,i,1),L7(6,i,2),L7(6,i,1)+L7(6,i,2)
       write(26,*) t,L7(7,i,1),L7(7,i,2),L7(7,i,1)+L7(7,i,2)
       write(27,*) t,L7(8,i,1),L7(8,i,2),L7(8,i,1)+L7(8,i,2)
       write(28,*) t,L7(9,i,1),L7(9,i,2),L7(9,i,1)+L7(9,i,2)
       write(29,*) t,L7(10,i,1),L7(10,i,2),L7(10,i,1)+L7(10,i,2)
    end do
    close(10)
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    close(16)
    close(17)
    close(18)
    close(19)
    close(20)
    close(21)
    close(22)
    close(23)
    close(24)
    close(25)
    close(26)
    close(27)
    close(28)
    close(29)
  end subroutine fig6_7

  !(subroutine output: To output calculated stress data!)
  subroutine output(s11,s12,s13,s21,s22,s23,s31,s32,s33,ex_num,dat_num,vtk_num)
    implicit none
    integer :: i,j
    integer, intent(in) :: ex_num,dat_num,vtk_num
    real(8), dimension(Mmax,Nmax), intent(inout) :: s11,s12,s13,s21,s22,s23,s31,s32,s33
    real(8) :: t,xx,yy,zz
    if (dat_num==1 .and. ex_num/=1) then
       open(10,file='./data/sigma11.dat', status='replace')
       open(11,file='./data/sigma12.dat', status='replace')
       open(12,file='./data/sigma13.dat', status='replace')
       open(13,file='./data/sigma21.dat', status='replace')
       open(14,file='./data/sigma22.dat', status='replace')
       open(15,file='./data/sigma23.dat', status='replace')
       open(16,file='./data/sigma31.dat', status='replace')
       open(17,file='./data/sigma32.dat', status='replace')
       open(18,file='./data/sigma33.dat', status='replace')
    else if (dat_num==1) then
       open(10,file='./data/sigma11ex.dat', status='replace')
       open(11,file='./data/sigma12ex.dat', status='replace')
       open(12,file='./data/sigma13ex.dat', status='replace')
       open(13,file='./data/sigma21ex.dat', status='replace')
       open(14,file='./data/sigma22ex.dat', status='replace')
       open(15,file='./data/sigma23ex.dat', status='replace')
       open(16,file='./data/sigma31ex.dat', status='replace')
       open(17,file='./data/sigma32ex.dat', status='replace')
       open(18,file='./data/sigma33ex.dat', status='replace')
    end if
    if (dat_num==1) then
       do i=1,Nmax
          do j=1,MMax
             t=Dt*(j-1)
             xx=xg(i); yy=yg(i); zz=zg(i)
             write(10,*) t,xx,yy,zz,s11(j,i)
             write(11,*) t,xx,yy,zz,s12(j,i)
             write(12,*) t,xx,yy,zz,s13(j,i)
             write(13,*) t,xx,yy,zz,s21(j,i)
             write(14,*) t,xx,yy,zz,s22(j,i)
             write(15,*) t,xx,yy,zz,s23(j,i)
             write(16,*) t,xx,yy,zz,s31(j,i)
             write(17,*) t,xx,yy,zz,s32(j,i)
             write(18,*) t,xx,yy,zz,s33(j,i)
          end do
          write(10,*)
          write(11,*)
          write(12,*)
          write(13,*)
          write(14,*)
          write(15,*)
          write(16,*)
          write(17,*)
          write(18,*)
       end do
    end if
    close(10)
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    close(16)
    close(17)
    close(18)
  end subroutine output

  subroutine dattovtk() !(s11,s12,s13,s21,s22,s23,s31,s32,s33)
    implicit none
    integer :: i,j,k,l,gn,ind,num
    character(len=40) :: file11,file12,flie13,file21,file22,file23,flie31,flie32,file33
    real(8), dimension(Mmax,Nmax) :: s11!,s12,s13,s21,s22,s23,s31,s32,s33
    real(8) :: t,xx,yy,zz
    !open(10,file="data/sigma11ex.dat")
    open(10,file="data/sigma11.dat")
    do i=1,Nmax
       do j=1,Mmax
          read(10,*) t,xx,yy,zz,s11(j,i)
       end do
       read(10,*)
    end do
    gn=(Nx+2)*Ny*Nz
    do i=1,Mmax
       !write(file11,"(a,i3.3,a)") "vtk/sigma11ex_",i,".vtk"
       write(file11,"(a,i3.3,a)") "vtk/sigma11_",i,".vtk"
       open(11,file=file11)
       write(11,"('# vtk DataFile Version 4.0')")
       write(11,"('sigma11')")
       write(11,"('ASCII')")
       write(11,"('DATASET UNSTRUCTURED_GRID')")
       write(11,"('POINTS',i8,'float')") gn
       do k=1,Nz !coordinates of the grids
          do j=1,Ny
             ind=(k-1)*(Nx*Ny)+(j-1)*Nx+1
             write(11,"(3(f9.5,1x))") xa(ind),ya(ind),za(ind)
             write(11,"(3(f9.5,1x))") xb(ind),yb(ind),zb(ind)
             write(11,"(3(f9.5,1x))") xc(ind),yc(ind),zc(ind)
             do l=2,Nx
                if (mod(l,2)==0) then
                   write(11,"(3(f9.5,1x))") xb(ind-1+l),yb(ind-1+l),zb(ind-1+l)
                else
                   write(11,"(3(f9.5,1x))") xc(ind-1+l),yc(ind-1+l),zc(ind-1+l)
                end if
             end do
          end do
       end do
       write(11,"(a,i8,i8)") "CELLS",Nmax,Nmax*4
       do k=1,Nz !calculate grid number for each cell (0~gn)
          do j=1,Ny
             do l=1,Nx
                if (l==1) then
                   num=(k-1)*(Nx+2)*Ny+(j-1)*(Nx+2)
                   write(11,"(4i8)") 3,num,num+1,num+2
                else if (l==2) then
                   num=(k-1)*(Nx+2)*Ny+(j-1)*(Nx+2)+2
                   write(11,"(4i8)") 3,num,num+1,num-2
                else if (l==3) then
                   num=(k-1)*(Nx+2)*Ny+(j-1)*(Nx+2)+3
                   write(11,"(4i8)") 3,num,num-1,num+1
                else
                   if (mod(l,2)==0) then
                      num=(k-1)*(Nx+2)*Ny+(j-1)*(Nx+2)+l
                      write(11,"(4i8)") 3,num,num+1,num-1
                   else
                      num=(k-1)*(Nx+2)*Ny+(j-1)*(Nx+2)+l
                      write(11,"(4i8)") 3,num,num-1,num+1
                   end if
                end if
             end do
          end do
       end do
       write(11,"(a,i8)") "CELL_TYPES", Nmax
       do j=1,Nmax
          write(11,*) 5 !triangle
       end do
       write(11,"(a,i8)") "CELL_DATA", Nmax
       write(11,"('SCALARS sigma11 float')") 
       write(11,"('LOOKUP_TABLE default')")
       do j=1,Nmax
          write(11,*) s11(i,j)
       end do
    end do
    close(10)
    close(11)
  end subroutine dattovtk
  
end module test_module
