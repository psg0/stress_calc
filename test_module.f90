module test_module
  use setting_module
  use kernel_module
  implicit none
  private
  public fig6_7
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
    !No1 (a)L111 (ABC)
    call Kdyns111(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(1,i,1)=Kdynf111(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t)
       L7(1,i,1)=Kdynf111(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t1)
    end do
    !No2 (b)L121 (ABC)
    call Kdyns121(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(2,i,1)=Kdynf121(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t)
       L7(2,i,1)=Kdynf121(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t1)
    end do
    !No3 (c)L221 (ABC)
    call Kdyns221(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(3,i,1)=Kdynf221(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t)
       L7(3,i,1)=Kdynf221(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t1)
    end do
    !No4 (d)L311 (ABC)
    call Kdyns311(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(4,i,1)=Kdynf311(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t)
       L7(4,i,1)=Kdynf311(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t1)
    end do
    !No5 (e)L321 (ABC)
    call Kdyns321(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(5,i,1)=Kdynf321(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t)
       L7(5,i,1)=Kdynf321(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t1)
    end do
    !No6 (f)L331 (ABC)
    call Kdyns331(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(6,i,1)=Kdynf331(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t)
       L7(6,i,1)=Kdynf331(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t1)
    end do
    !No7 (g)L113 (ABC)
    call Kdyns113(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(7,i,1)=Kdynf113(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t)
       L7(7,i,1)=Kdynf113(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t1)
    end do
    !No8 (h)L123 (ABC)
    call Kdyns123(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(8,i,1)=Kdynf123(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t)
       L7(8,i,1)=Kdynf123(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t1)
    end do
    !No9 (i)L313 (ABC)
    call Kdyns313(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(9,i,1)=Kdynf313(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t)
       L7(9,i,1)=Kdynf313(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t1)
    end do
    !No10 (j)L333 (ABC)
    call Kdyns333(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(10,i,1)=Kdynf333(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t)
       L7(10,i,1)=Kdynf333(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t1)
    end do

    !===CDA===
    xsa(1)=1.0_8; xsa(2)=1.0_8; xsb(1)=0.0_8; xsb(2)=1.0_8
    xsc(1)=0.0_8; xsc(2)=0.0_8; xs3=0.0_8
    !No1 (a)L111 (CDA)
    call Kdyns111(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(1,i,2)=Kdynf111(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t)
       L7(1,i,2)=Kdynf111(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t1)
    end do
    !No2 (b)L121 (CDA)
    call Kdyns121(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(2,i,2)=Kdynf121(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t)
       L7(2,i,2)=Kdynf121(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t1)
    end do
    !No3 (c)L221 (CDA)
    call Kdyns221(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(3,i,2)=Kdynf221(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t)
       L7(3,i,2)=Kdynf221(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t1)
    end do
    !No4 (d)L311 (CDA)
    call Kdyns311(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(4,i,2)=Kdynf311(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t)
       L7(4,i,2)=Kdynf311(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t1)
    end do
    !No5 (e)L321 (CDA)
    call Kdyns321(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(5,i,2)=Kdynf321(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t)
       L7(5,i,2)=Kdynf321(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t1)
    end do
    !No6 (f)L331 (CDA)
    call Kdyns331(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(6,i,2)=Kdynf331(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t)
       L7(6,i,2)=Kdynf331(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t1)
    end do
    !No7 (g)L113 (CDA)
    call Kdyns113(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(7,i,2)=Kdynf113(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t)
       L7(7,i,2)=Kdynf113(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t1)
    end do
    !No8 (c)L123 (CDA)
    call Kdyns123(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(8,i,2)=Kdynf123(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t)
       L7(8,i,2)=Kdynf123(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t1)
    end do
    !No9 (d)L313 (CDA)
    call Kdyns313(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(9,i,2)=Kdynf313(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t)
       L7(9,i,2)=Kdynf313(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t1)
    end do
    !No10 (e)L333 (CDA)
    call Kdyns333(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3)
    do i=0,n
       t=dtime*i
       L6(10,i,2)=Kdynf333(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t)
       L7(10,i,2)=Kdynf333(xr(1),xr(2),xr(3),xsa(1),xsa(2),xsb(1),xsb(2),xsc(1),xsc(2),xs3,t,t1)
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
  end subroutine fig6_7
end module test_module
