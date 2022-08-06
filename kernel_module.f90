module kernel_module !Tada 2006
  use setting_module, only: pi,pi2,pi3,hpi,Cl=>alpha,Ct=>beta,Cl2=>alpha2,Ct2=>beta2,Cl3=>alpha3,Ct3=>beta3 
  implicit none
  private
  public Kdyns111, Kdyns112, Kdyns113, Kdyns121, Kdyns122, Kdyns123
  public Kdyns131, Kdyns132, Kdyns133, Kdyns211, Kdyns212, Kdyns213
  public Kdyns221, Kdyns222, Kdyns223, Kdyns231, Kdyns232, Kdyns233
  public Kdyns311, Kdyns312, Kdyns313, Kdyns321, Kdyns322, Kdyns323
  public Kdyns331, Kdyns332, Kdyns333
  public Kdynf111, Kdynf112, Kdynf113, Kdynf121, Kdynf122, Kdynf123
  public Kdynf131, Kdynf132, Kdynf133, Kdynf211, Kdynf212, Kdynf213
  public Kdynf221, Kdynf222, Kdynf223, Kdynf231, Kdynf232, Kdynf233
  public Kdynf311, Kdynf312, Kdynf313, Kdynf321, Kdynf322, Kdynf323
  public Kdynf331, Kdynf332, Kdynf333
  real(8),parameter :: p=Ct/Cl, p2=p*p, Caa3=(1.0_8-2.0_8*p2)/p
  real(8) :: lds111(0:4-1,0:6-1)
  real(8) :: lds112(0:4-1,0:6-1)
  real(8) :: lds113(0:7-1,0:6-1)
  real(8) :: lds121(0:4-1,0:6-1)
  real(8) :: lds122(0:4-1,0:6-1)
  real(8) :: lds123(0:6-1,0:6-1)
  real(8) :: lds131(0:7-1,0:6-1)
  real(8) :: lds132(0:6-1,0:6-1)
  real(8) :: lds133(0:4-1,0:6-1)
  real(8) :: lds211(0:4-1,0:6-1)
  real(8) :: lds212(0:4-1,0:6-1)
  real(8) :: lds213(0:6-1,0:6-1)
  real(8) :: lds221(0:4-1,0:6-1)
  real(8) :: lds222(0:4-1,0:6-1)
  real(8) :: lds223(0:7-1,0:6-1)
  real(8) :: lds231(0:6-1,0:6-1)
  real(8) :: lds232(0:7-1,0:6-1)
  real(8) :: lds233(0:4-1,0:6-1)
  real(8) :: lds311(0:7-1,0:6-1)
  real(8) :: lds312(0:6-1,0:6-1)
  real(8) :: lds313(0:4-1,0:6-1)
  real(8) :: lds321(0:6-1,0:6-1)
  real(8) :: lds322(0:7-1,0:6-1)
  real(8) :: lds323(0:4-1,0:6-1)
  real(8) :: lds331(0:4-1,0:6-1)
  real(8) :: lds332(0:4-1,0:6-1)
  real(8) :: lds333(0:5-1,0:6-1)
  
 
  real(8) :: zii1s(0:3-1)
  real(8) :: zi1s(0:3-1)
  real(8) :: z2s(0:3-1)
  real(8) :: z3s(0:3-1)
contains

  ! (subroutine Rotcoordnt: calculate local coordinate values)
  subroutine Rotcoordnt(xr1,xr2,xr3,xsb1,xsb2,xsc1,xsc2,xs3,costa,sinta,zb1a,zc1a,z2a,z3a)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsb1,xsb2,xsc1,xsc2,xs3
    real(8),intent(out) :: costa,sinta,zb1a,zc1a,z2a,z3a
    real(8) :: xr1a,xr2a,xsb1a,xsb2a,xsc1a,xsc2a,bc1,bc2,bcr

    !xr : coordinates of receiver,  xs: coordinates of source   
    bc1=xsc1-xsb1
    bc2=xsc2-xsb2
    bcr=sqrt(bc1*bc1+bc2*bc2)
    costa=bc1/bcr
    sinta=bc2/bcr
    !rotation
    xr1a=xr1*(costa)+xr2*(sinta)
    xr2a=-xr1*(sinta)+xr2*(costa)
    xsb1a=xsb1*(costa)+xsb2*(sinta)
    xsb2a=-xsb1*(sinta)+xsb2*(costa)
    xsc1a=xsc1*(costa)+xsc2*(sinta)
    xsc2a=-xsc1*(sinta)+xsc2*(costa)
    !relative location in the local coordinate system
    zb1a=xr1a-xsb1a
    zc1a=xr1a-xsc1a
    z2a=xr2a-xsb2a ! or  xr2a-xsc2a
    z3a=xr3-xs3
  end subroutine Rotcoordnt

  real(8) function Heviw(x)
    implicit none  
    real(8),intent(in) :: x
    real(8) :: ret

    if(x>0.0_8) then
       ret=1.0_8
    else
       ret=0.0_8
    end if
    Heviw=ret
  end function Heviw

  real(8) function Hevix(x)
    implicit none  
    real(8),intent(in) :: x
    real(8) :: ret

    if(x>0.0_8) then
       ret=1.0_8
    else if (x<0.0_8) then
       ret=0.0_8
    else
       ret=0.50_8
    end if
    Hevix=ret
  end function Hevix

  real(8) function Sgn(x)
    implicit none
    real(8),intent(in) :: x
    real(8) :: ret

    if (x>0.0_8) then
       ret=1.0_8
    else if (x<0.0_8) then
       ret=-1.0_8
    else
       ret=0.0_8
    end if
    Sgn=ret
  end function Sgn

  ! (subroutine Calg: calculate the time dependent part of L0)
  subroutine Calg142(z1a,z2a,t,rbc,ra,g1l,g2l,g4l,g1t,g2t,g4t)
    implicit none
    real(8),intent(in) :: z1a,z2a,t,rbc,ra
    real(8),intent(out) :: g1l,g2l,g4l,g1t,g2t,g4t
    real(8) :: t2,t3,rbc2,ra2,ra3
    real(8) :: hvrat,hvral,hvrbct,hvrbcl,hvza,sgza,hvsgt,hvsgl
    real(8) :: trbct,trbcl,sqtrbct,sqtrbcl
    
    t2=t*t
    t3=t2*t
    rbc2=rbc*rbc
    ra2=ra*ra
    ra3=ra2*ra

    hvral=Heviw(t-ra/Cl)
    hvrat=Heviw(t-ra/Ct)
    hvrbcl=Heviw(t-rbc/Cl)
    hvrbct=Heviw(t-rbc/Ct)
    hvza=Hevix(z1a)
    sgza=Sgn(z1a)
    
    trbcl=t2-rbc2/Cl2
    trbct=t2-rbc2/Ct2
    if(trbcl>0.0_8) then
       sqtrbcl=sqrt(trbcl)
    else
       sqtrbcl=0.0_8
    end if
    if(trbct>0.0_8) then
       sqtrbct=sqrt(trbct)
    else
       sqtrbct=0.0_8
    end if
    
    hvsgl=2.0_8*hvza*hvrbcl-sgza*hvral
    hvsgt=2.0_8*hvza*hvrbct-sgza*hvrat

    g1l=trbcl*sqtrbcl*hvsgl+ &
         (z1a*(2.0_8*ra2+rbc2)/(2.0_8*ra3)*t3-3.0_8*z1a*rbc2/(2.0_8*Cl2*ra)*t)*hvral
    g1t=trbct*sqtrbct*hvsgt+ &
         (z1a*(2.0_8*ra2+rbc2)/(2.0_8*ra3)*t3-3.0_8*z1a*rbc2/(2.0d0*Ct2*ra)*t)*hvrat
    g2l=sqtrbcl*hvsgl+z1a/ra*t*hvral
    g2t=sqtrbct*hvsgt+z1a/ra*t*hvrat
    g4l=(t3-ra2/Cl2*t)*hvral
    g4t=(t3-ra2/Ct2*t)*hvrat
  end subroutine Calg142

  subroutine Calg3(z1a,z2a,t,rbc,ra,g3l,g3t)
    implicit none
    real(8),intent(in) :: z1a,z2a,t,rbc,ra
    real(8),intent(out) :: g3l,g3t
    real(8) :: t2,rbc2,hvral,hvrat,hvrbcl,hvrbct,hvza,sgza
    real(8) :: trbcl,trbct,sqtrbcl,sqtrbct,hvsgl,hvsgt,atanl,atant,atanz

    t2=t*t
    rbc2=rbc*rbc

    hvral=Heviw(t-ra/Cl)
    hvrat=Heviw(t-ra/Ct)
    hvrbcl=Heviw(t-rbc/Cl)
    hvrbct=Heviw(t-rbc/Ct)
    hvza=Hevix(z1a)
    sgza=Sgn( z1a )
    
    trbcl=t2-rbc2/Cl2
    trbct=t2-rbc2/Ct2
    if(trbcl>0.0_8) then
       sqtrbcl=sqrt(trbcl)
    else
       sqtrbcl=0.0_8
    end if
    if(trbct>0.0_8) then
       sqtrbct=sqrt(trbct)
    else
       sqtrbct=0.0_8
    end if
    hvsgl=(2.0_8*hvza*hvrbcl-sgza*hvral)
    hvsgt=(2.0_8*hvza*hvrbct-sgza*hvrat)
    if(hvsgl>0.0_8) then
       if(z2a/=0.0_8) then
          atanl=atan(sqtrbcl/(z2a/Cl))
          atant=atan(sqtrbct/(z2a/Ct))
       else
          atanl=hpi
          atant=hpi
       end if
    else
       atanl=0.0_8
       atant=0.0_8
    end if
    if(hvral>0.0_8) then
       if(z2a/=0.0_8) then
          atanz=atan(z1a/z2a)
       else
          atanz=hpi
       end if
    else
       atanz=0.0_8
    end if
    
    g3l=atanl*hvsgl+atanz*hvral
    g3t=atant*hvsgt+atanz*hvrat
  end subroutine Calg3

  subroutine Calg5(z1a,z2a,t,rbc,ra,g5l,g5t)
    implicit none
    real(8),intent(in) :: z1a,z2a,t,rbc,ra
    real(8),intent(out) :: g5l,g5t
    real(8) :: ral,rat,hvrat,hvral

    ral=t-ra/Cl
    rat=t-ra/Ct
    hvral=Heviw(ral)
    hvrat=Heviw(rat)

    g5l=ral*hvral
    g5t=rat* hvrat
  end subroutine Calg5

  ! (subroutine Ldyns: calculate the time independent part of L0 only once and assign the value to lds)
  subroutine Ldyns111(z1a,z2a,z3a,costa,sinta,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,costa,sinta
    integer,intent(in) :: i
    real(8) :: rbc,ra 
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2,sinta3
    real(8) :: tmp

    z2a2=z2a*z2a
    z3a2=z3a*z3a
    rbc=sqrt(z2a2+z3a2)
    rbc2=rbc*rbc
    rbc6=rbc2*rbc2*rbc2
    ra=sqrt(z1a*z1a+rbc2)
    ra2=ra*ra
    ra5=ra2*ra2*ra
    sinta2=sinta*sinta
    sinta3=sinta2*sinta
    costa2=costa*costa

    lds111(0,i)=-2.0_8*Ct3*z3a*(3.0_8*z2a2-z3a2)/(3.0_8*pi*rbc6)*sinta3
    lds111(1,i)=-Ct3*z3a/(pi*ra5)*sinta*(z1a*z2a2/rbc2*sinta2+2.0_8*z2a*sinta*costa-z1a*costa2)
    tmp=-Ct*z3a/(pi*rbc2)*sinta
    lds111(2,i)=tmp*(1.0_8-2.0_8*p2+2.0_8*p2*z2a2/rbc2*sinta2)
    lds111(3,i)=tmp*(1.0_8-2.0_8*z2a2/rbc2*sinta2)
  end subroutine Ldyns111

  subroutine Ldyns112(z1a,z2a,z3a,costa,sinta,i) ! derive from Ldyns221
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,costa,sinta
    integer,intent(in) :: i
    real(8) :: rbc,ra
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2
    real(8) :: tmp

    z2a2=z2a*z2a
    z3a2=z3a*z3a
    rbc=sqrt(z2a2+z3a2)
    rbc2=rbc*rbc
    rbc6=rbc2*rbc2*rbc2
    ra=sqrt(z1a*z1a+rbc2)
    ra2=ra*ra
    ra5=ra2*ra2*ra
    sinta2=sinta*sinta
    costa2=costa*costa

    lds112(0,i)=2.0_8*Ct3*z3a*(3.0_8*z2a2-z3a2)/(3.0_8*pi*rbc6)*costa*sinta2
    lds112(1,i)=Ct3*z3a/(pi*ra5)*costa*(z1a*z2a2/rbc2*sinta2+2.0_8*z2a*sinta*costa-z1a*costa2)
    tmp=Ct*z3a/(pi*rbc2)*costa
    lds112(2,i)=tmp*(1.0_8-2.0_8*p2+2.0_8*p2*z2a2/rbc2*sinta2)
    lds112(3,i)=tmp*(-2.0_8*z2a2/rbc2*sinta2)
  end subroutine Ldyns112

  subroutine Ldyns113(z1a,z2a,z3a,costa,sinta,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,costa,sinta
    integer,intent(in) :: i
    real(8) :: rbc,ra 
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2
    real(8) :: tmp1,tmp2

    z2a2=z2a*z2a
    z3a2=z3a*z3a
    rbc=sqrt(z2a2+z3a2)
    rbc2=rbc*rbc
    rbc6=rbc2*rbc2*rbc2
    ra=sqrt(z1a*z1a+rbc2)
    ra2=ra*ra
    ra5=ra2*ra2*ra
    sinta2=sinta*sinta
    costa2=costa*costa

    lds113(0,i)=2.0_8*Ct3*z2a*(3.0_8*z3a2-z2a2)/(3.0_8*pi*rbc6)*sinta2
    lds113(1,i)=Ct3/(pi*ra5)*sinta*(z1a*z2a*z3a2/rbc2*sinta+(z3a2-ra2/3.0_8)*costa)
    tmp1=-Ct*z2a/(pi*rbc2)
    lds113(2,i)=tmp1*((1.0_8-2.0_8*p2)*costa2-2.0_8*p2*z3a2/rbc2*sinta2)
    lds113(3,i)=tmp1*(2.0_8*z3a2/rbc2*sinta2)
    tmp2=Ct/(pi*ra)*sinta*costa
    lds113(4,i)=tmp2*(1.0_8-4.0_8/3.0_8*p2)
    lds113(5,i)=tmp2*(-2.0_8/3.0_8)
    lds113(6,i)=(1.0_8-2.0_8*p2)/(2.0_8*pi*p)
  end subroutine Ldyns113

  subroutine Ldyns121(z1a,z2a,z3a,costa,sinta,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,costa,sinta
    integer,intent(in) :: i
    real(8) :: rbc,ra
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2
    real(8) :: tmp

    z2a2=z2a*z2a
    z3a2=z3a*z3a
    rbc=sqrt(z2a2+z3a2)
    rbc2=rbc*rbc
    rbc6=rbc2*rbc2*rbc2
    ra=sqrt(z1a*z1a+rbc2)
    ra2=ra*ra
    ra5=ra2*ra2*ra
    sinta2=sinta*sinta
    costa2=costa*costa

    lds121(0,i)=2.0_8*Ct3*z3a*(3.0_8*z2a2-z3a2)/(3.0_8*pi*rbc6)*sinta2*costa
    lds121(1,i)=Ct3*z3a/(pi*ra5)*costa*(z1a*z2a2/rbc2*sinta2+2.0_8*z2a*sinta*costa-z1a*costa2)
    tmp=Ct*z3a/(pi*rbc2)*costa
    lds121(2,i)=tmp*(2.0_8*p2*z2a2/rbc2*sinta2)
    lds121(3,i)=tmp*(0.5_8-2.0_8*z2a2/rbc2*sinta2)
  end subroutine Ldyns121

  subroutine Ldyns122(z1a,z2a,z3a,costa,sinta,i) !derive from Ldyns211=121
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,costa,sinta
    integer,intent(in) :: i
    real(8) :: rbc,ra
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2
    real(8) :: tmp

    z2a2=z2a*z2a
    z3a2=z3a*z3a
    rbc=sqrt(z2a2+z3a2)
    rbc2=rbc*rbc
    rbc6=rbc2*rbc2*rbc2
    ra=sqrt(z1a*z1a+rbc2)
    ra2=ra*ra
    ra5=ra2*ra2*ra
    sinta2=sinta*sinta
    costa2=costa*costa

    lds122(0,i)=-2.0_8*Ct3*z3a*(3.0_8*z2a2-z3a2)/(3.0_8*pi*rbc6)*costa2*sinta
    lds122(1,i)=-Ct3*z3a/(pi*ra5)*sinta*(z1a*z2a2/rbc2*costa2-2.0_8*z2a*sinta*costa-z1a*sinta2)
    tmp=-Ct*z3a/(pi*rbc2)*sinta
    lds122(2,i)=tmp*(2.0_8*p2*z2a2/rbc2*costa2)
    lds122(3,i)=tmp*(0.5_8-2.0_8*z2a2/rbc2*costa2)
  end subroutine Ldyns122

  subroutine Ldyns123(z1a,z2a,z3a,costa,sinta,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,costa,sinta
    integer,intent(in) :: i
    real(8) :: rbc,ra
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2
    real(8) :: tmp1,tmp2
    z2a2=z2a*z2a
    z3a2=z3a*z3a
    rbc=sqrt(z2a2+z3a2)
    rbc2=rbc*rbc
    rbc6=rbc2*rbc2*rbc2
    ra=sqrt(z1a*z1a+rbc2)
    ra2=ra*ra
    ra5=ra2*ra2*ra
    sinta2=sinta*sinta
    costa2=costa*costa

    lds123(0,i)=-2.0_8*Ct3*z2a*(3.0_8*z3a2-z2a2)/(3.0_8*pi*rbc6)*sinta*costa
    lds123(1,i)=-Ct3/(pi*ra5)*costa*(z1a*z2a*z3a2/rbc2*sinta+(z3a2-ra2/3.0_8)*costa)
    tmp1=-Ct*z2a/(pi*rbc2)*sinta*costa
    lds123(2,i)=tmp1*(1.0_8-2.0_8*p2+2.0_8*p2*z3a2/rbc2)
    lds123(3,i)=tmp1*(-2.0_8*z3a2/rbc2)
    tmp2=-Ct/(pi*ra)*costa2
    lds123(4,i)=tmp2*(1.0_8-4.0_8/3.0_8*p2)
    lds123(5,i)=tmp2*(-2.0_8/3.0_8)
  end subroutine Ldyns123

  subroutine Ldyns131(z1a,z2a,z3a,costa,sinta,i) !derive from Ldyns311
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,costa,sinta
    integer,intent(in) :: i
    real(8) :: rbc,ra 
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2
    real(8) :: tmp1,tmp2
    
    z2a2=z2a*z2a
    z3a2=z3a*z3a
    rbc=sqrt(z2a2+z3a2)
    rbc2=rbc*rbc
    rbc6=rbc2*rbc2*rbc2
    ra=sqrt(z1a*z1a+rbc2)
    ra2=ra*ra
    ra5=ra2*ra2*ra
    sinta2=sinta*sinta
    costa2=costa*costa

    lds131(0,i)=2.0_8*Ct3*z2a*(3.0_8*z3a2-z2a2)/(3.0_8*pi*rbc6)*sinta2
    lds131(1,i)=Ct3/(pi*ra5)*sinta*(z1a*z2a*z3a2/rbc2*sinta+(z3a2-ra2/3.0_8 )*costa)
    tmp1=Ct*z2a/(pi*rbc2)
    lds131(2,i)=tmp1*(2.0_8*p2*z3a2/rbc2*sinta2)
    lds131(3,i)=tmp1*(-(0.5_8*costa2+2.0_8*z3a2/rbc2*sinta2))
    tmp2=Ct/(pi*ra)*sinta*costa
    lds131(4,i)=tmp2*(2.0_8*p2/3.0_8)
    lds131(5,i)=tmp2*(-1.0_8/6.0_8)
    lds131(6,i)=0.5_8/pi
  end subroutine Ldyns131

  subroutine Ldyns132(z1a,z2a,z3a,costa,sinta,i) !derive from Ldyns312
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,costa,sinta
    integer,intent(in) :: i
    real(8) :: rbc,ra
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2
    real(8) :: tmp1,tmp2

    z2a2=z2a*z2a
    z3a2=z3a*z3a
    rbc=sqrt(z2a2+z3a2)
    rbc2=rbc*rbc
    rbc6=rbc2*rbc2*rbc2
    ra=sqrt(z1a*z1a+rbc2)
    ra2=ra*ra
    ra5=ra2*ra2*ra
    sinta2=sinta*sinta
    costa2=costa*costa

    lds132(0,i)=-2.0_8*Ct3*z2a*(3.0_8*z3a2-z2a2)/(3.0_8*pi*rbc6)*sinta*costa
    lds132(1,i)=-Ct3/(pi*ra5)*costa*(z1a*z2a*z3a2/rbc2*sinta+(z3a2-ra2/3.0_8)*costa)
    tmp1=-Ct*z2a/(pi*rbc2)*sinta*costa
    lds132(2,i)=tmp1*(2.0_8*p2*z3a2/rbc2)
    lds132(3,i)=tmp1*(0.5_8-2.0_8*z3a2/rbc2)
    tmp2=-Ct/(pi*ra)*costa2
    lds132(4,i)=tmp2*(2.0_8*p2/3.0_8)
    lds132(5,i)=tmp2*(-1.0_8/6.0_8)
  end subroutine Ldyns132

  subroutine Ldyns133(z1a,z2a,z3a,costa,sinta,i) !derive from Ldyns313
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,costa,sinta
    integer,intent(in) :: i
    real(8) :: rbc,ra 
    real(8) :: z2a2,z3a2,z3a3,rbc2,rbc6,ra2,ra5
    real(8) :: tmp

    z2a2=z2a*z2a
    z3a2=z3a*z3a
    z3a3=z3a*z3a2
    rbc=sqrt(z2a2+z3a2)
    rbc2=rbc*rbc
    rbc6=rbc2*rbc2*rbc2
    ra=sqrt(z1a*z1a+rbc2)
    ra2=ra*ra
    ra5=ra2*ra2*ra

    lds133(0,i)=2.0_8*Ct3*z3a*(3.0_8*z2a2-z3a2)/(3.0_8*pi*rbc6)*sinta
    lds133(1,i)=-Ct3*z1a*z3a3/(pi*ra5*rbc2)*sinta
    tmp=-Ct*z3a/(pi*rbc2)*sinta
    lds133(2,i)=tmp*(1.0_8-2.0_8*p2+2.0_8*p2*z3a2/rbc2)
    lds133(3,i)=tmp*(1.0_8-2.0_8*z3a2/rbc2)
  end subroutine Ldyns133

  subroutine Ldyns211(z1a,z2a,z3a,costa,sinta,i) !derive from Ldyns121
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,costa,sinta
    integer,intent(in) :: i
    real(8) :: rbc,ra
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2
    real(8) :: tmp

    z2a2=z2a*z2a
    z3a2=z3a*z3a
    rbc=sqrt(z2a2+z3a2)
    rbc2=rbc*rbc
    rbc6=rbc2*rbc2*rbc2
    ra=sqrt(z1a*z1a+rbc2)
    ra2=ra*ra
    ra5=ra2*ra2*ra
    sinta2=sinta*sinta
    costa2=costa*costa

    lds211(0,i)=2.0_8*Ct3*z3a*(3.0_8*z2a2-z3a2)/(3.0_8*pi*rbc6)*sinta2*costa
    lds211(1,i)=Ct3*z3a/(pi*ra5)*costa*(z1a*z2a2/rbc2*sinta2+2.0_8*z2a*sinta*costa-z1a*costa2)
    tmp=Ct*z3a/(pi*rbc2)*costa
    lds211(2,i)=tmp*(2.0_8*p2*z2a2/rbc2*sinta2)
    lds211(3,i)=tmp*(0.5_8-2.0_8*z2a2/rbc2*sinta2)
  end subroutine Ldyns211

  subroutine Ldyns212(z1a,z2a,z3a,costa,sinta,i) !derive from Ldyns122
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,costa,sinta
    integer,intent(in) :: i
    real(8) :: rbc,ra
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2
    real(8) :: tmp

    z2a2=z2a*z2a
    z3a2=z3a*z3a
    rbc=sqrt(z2a2+z3a2)
    rbc2=rbc*rbc
    rbc6=rbc2*rbc2*rbc2
    ra=sqrt(z1a*z1a+rbc2)
    ra2=ra*ra
    ra5=ra2*ra2*ra
    sinta2=sinta*sinta
    costa2=costa*costa

    lds212(0,i)=-2.0_8*Ct3*z3a*(3.0_8*z2a2-z3a2)/(3.0_8*pi*rbc6)*costa2*sinta
    lds212(1,i)=-Ct3*z3a/(pi*ra5)*sinta*(z1a*z2a2/rbc2*costa2-2.0_8*z2a*sinta*costa-z1a*sinta2)
    tmp=-Ct*z3a/(pi*rbc2)*sinta
    lds212(2,i)=tmp*(2.0_8*p2*z2a2/rbc2*costa2)
    lds212(3,i)=tmp*(0.5_8-2.0_8*z2a2/rbc2*costa2)
  end subroutine Ldyns212

  subroutine Ldyns213(z1a,z2a,z3a,costa,sinta,i) !derive from Ldyns123
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,costa,sinta
    integer,intent(in) :: i
    real(8) :: rbc,ra
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2
    real(8) :: tmp1,tmp2
    z2a2=z2a*z2a
    z3a2=z3a*z3a
    rbc=sqrt(z2a2+z3a2)
    rbc2=rbc*rbc
    rbc6=rbc2*rbc2*rbc2
    ra=sqrt(z1a*z1a+rbc2)
    ra2=ra*ra
    ra5=ra2*ra2*ra
    sinta2=sinta*sinta
    costa2=costa*costa

    lds213(0,i)=-2.0_8*Ct3*z2a*(3.0_8*z3a2-z2a2)/(3.0_8*pi*rbc6)*sinta*costa
    lds213(1,i)=-Ct3/(pi*ra5)*costa*(z1a*z2a*z3a2/rbc2*sinta+(z3a2-ra2/3.0_8)*costa)
    tmp1=-Ct*z2a/(pi*rbc2)*sinta*costa
    lds213(2,i)=tmp1*(1.0_8-2.0_8*p2+2.0_8*p2*z3a2/rbc2)
    lds213(3,i)=tmp1*(-2.0_8*z3a2/rbc2)
    tmp2=-Ct/(pi*ra)*costa2
    lds213(4,i)=tmp2*(1.0_8-4.0_8/3.0_8*p2)
    lds213(5,i)=tmp2*(-2.0_8/3.0_8)
  end subroutine Ldyns213
  
  subroutine Ldyns221(z1a,z2a,z3a,costa,sinta,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,costa,sinta
    integer,intent(in) :: i
    real(8) :: rbc,ra
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2
    real(8) :: tmp

    z2a2=z2a*z2a
    z3a2=z3a*z3a
    rbc=sqrt(z2a2+z3a2)
    rbc2=rbc*rbc
    rbc6=rbc2*rbc2*rbc2
    ra=sqrt(z1a*z1a+rbc2)
    ra2=ra*ra
    ra5=ra2*ra2*ra
    sinta2=sinta*sinta
    costa2=costa*costa

    lds221(0,i)=-2.0_8*Ct3*z3a*(3.0_8*z2a2-z3a2)/(3.0_8*pi*rbc6)*sinta*costa2
    lds221(1,i)=-Ct3*z3a/(pi*ra5)*sinta*(z1a*z2a2/rbc2*costa2-2.0_8*z2a*sinta*costa-z1a*sinta2)
    tmp=-Ct*z3a/(pi*rbc2)*sinta
    lds221(2,i)=tmp*(1.0_8-2.0_8*p2+2.0_8*p2*z2a2/rbc2*costa2)
    lds221(3,i)=tmp*(-2.0_8*z2a2/rbc2*costa2)
  end subroutine Ldyns221

  subroutine Ldyns222(z1a,z2a,z3a,costa,sinta,i) !derive from Ldyns111
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,costa,sinta
    integer,intent(in) :: i
    real(8) :: rbc,ra 
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2,sinta3,costa3
    real(8) :: tmp

    z2a2=z2a*z2a
    z3a2=z3a*z3a
    rbc=sqrt(z2a2+z3a2)
    rbc2=rbc*rbc
    rbc6=rbc2*rbc2*rbc2
    ra=sqrt(z1a*z1a+rbc2)
    ra2=ra*ra
    ra5=ra2*ra2*ra
    sinta2=sinta*sinta
    sinta3=sinta2*sinta
    costa2=costa*costa
    costa3=costa2*costa

    lds222(0,i)=2.0_8*Ct3*z3a*(3.0_8*z2a2-z3a2)/(3.0_8*pi*rbc6)*costa3
    lds222(1,i)=Ct3*z3a/(pi*ra5)*costa*(z1a*z2a2/rbc2*costa2-2.0_8*z2a*sinta*costa-z1a*sinta2)
    tmp=Ct*z3a/(pi*rbc2)*costa
    lds222(2,i)=tmp*(1.0_8-2.0_8*p2+2.0_8*p2*z2a2/rbc2*costa2)
    lds222(3,i)=tmp*(1.0_8-2.0_8*z2a2/rbc2*costa2)
  end subroutine Ldyns222

  subroutine Ldyns223(z1a,z2a,z3a,costa,sinta,i) !derive from Ldyns113
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,costa,sinta
    integer,intent(in) :: i
    real(8) :: rbc,ra 
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2
    real(8) :: tmp1,tmp2

    z2a2=z2a*z2a
    z3a2=z3a*z3a
    rbc=sqrt(z2a2+z3a2)
    rbc2=rbc*rbc
    rbc6=rbc2*rbc2*rbc2
    ra=sqrt(z1a*z1a+rbc2)
    ra2=ra*ra
    ra5=ra2*ra2*ra
    sinta2=sinta*sinta
    costa2=costa*costa

    lds223(0,i)=2.0_8*Ct3*z2a*(3.0_8*z3a2-z2a2)/(3.0_8*pi*rbc6)*costa2
    lds223(1,i)=-Ct3/(pi*ra5)*costa*(-z1a*z2a*z3a2/rbc2*costa+(z3a2-ra2/3.0_8)*sinta)
    tmp1=-Ct*z2a/(pi*rbc2)
    lds223(2,i)=tmp1*((1.0_8-2.0_8*p2)*sinta2-2.0_8*p2*z3a2/rbc2*costa2)
    lds223(3,i)=tmp1*(2.0_8*z3a2/rbc2*costa2)
    tmp2=-Ct/(pi*ra)*sinta*costa
    lds223(4,i)=tmp2*(1.0_8-4.0_8/3.0_8*p2)
    lds223(5,i)=tmp2*(-2.0_8/3.0_8)
    lds223(6,i)=(1.0_8-2.0_8*p2)/(2.0_8*pi*p)
  end subroutine Ldyns223

  subroutine Ldyns231(z1a,z2a,z3a,costa,sinta,i) !derive from Ldyns321
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,costa,sinta
    integer,intent(in) :: i
    real(8) :: rbc,ra
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2
    real(8) :: tmp1,tmp2

    z2a2=z2a*z2a
    z3a2=z3a*z3a
    rbc=sqrt(z2a2+z3a2)
    rbc2=rbc*rbc
    rbc6=rbc2*rbc2*rbc2
    ra=sqrt(z1a*z1a+rbc2)
    ra2=ra*ra
    ra5=ra2*ra2*ra
    sinta2=sinta*sinta
    costa2=costa*costa

    lds231(0,i)=-2.0_8*Ct3*z2a*(3.0_8*z3a2-z2a2)/(3.0_8*pi*rbc6)*sinta*costa
    lds231(1,i)=-Ct3/(pi*ra5)*sinta*(z1a*z2a*z3a2/rbc2*costa-(z3a2-ra2/3.0_8)*sinta)
    tmp1=-Ct*z2a/(pi*rbc2)*sinta*costa
    lds231(2,i)=tmp1*(2.0_8*p2*z3a2/rbc2)
    lds231(3,i)=tmp1*(0.5_8-2.0_8*z3a2/rbc2)
    tmp2=Ct/(pi*ra)*sinta2
    lds231(4,i)=tmp2*(2.0_8*p2/3.0_8)
    lds231(5,i)=tmp2*(-1.0_8/6.0_8)
  end subroutine Ldyns231

  subroutine Ldyns232(z1a,z2a,z3a,costa,sinta,i) ! derive from Ldyns322
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,costa,sinta
    integer,intent(in) :: i
    real(8) :: rbc,ra
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2
    real(8) :: tmp1,tmp2
    
    z2a2=z2a*z2a
    z3a2=z3a*z3a
    rbc=sqrt(z2a2+z3a2)
    rbc2=rbc*rbc
    rbc6=rbc2*rbc2*rbc2
    ra=sqrt(z1a*z1a+rbc2)
    ra2=ra*ra
    ra5=ra2*ra2*ra
    sinta2=sinta*sinta
    costa2=costa*costa

    lds232(0,i)=2.0_8*Ct3*z2a*(3.0_8*z3a2-z2a2)/(3.0_8*pi*rbc6)*costa2
    lds232(1,i)=-Ct3/(pi*ra5)*costa*(-z1a*z2a*z3a2/rbc2*costa+(z3a2-ra2/3.0_8)*sinta)
    tmp1=Ct*z2a/(pi*rbc2)
    lds232(2,i)=tmp1*(2.0_8*p2*z3a2/rbc2*costa2)
    lds232(3,i)=tmp1*(-(0.5_8*sinta2+2.0_8*z3a2/rbc2*costa2))
    tmp2=-Ct/(pi*ra)*sinta*costa
    lds232(4,i)=tmp2*(2.0_8*p2/3.0_8)
    lds232(5,i)=tmp2*(-1.0_8/6.0_8)
    lds232(6,i)=0.5_8/pi
  end subroutine Ldyns232

  subroutine Ldyns233(z1a,z2a,z3a,costa,sinta,i) !derive from Ldyns323
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,costa,sinta
    integer,intent(in) :: i
    real(8) :: rbc,ra 
    real(8) :: z2a2,z3a2,z3a3,rbc2,rbc6,ra2,ra5
    real(8) :: tmp

    z2a2=z2a*z2a
    z3a2=z3a*z3a
    z3a3=z3a*z3a2
    rbc=sqrt(z2a2+z3a2)
    rbc2=rbc*rbc
    rbc6=rbc2*rbc2*rbc2
    ra=sqrt(z1a*z1a+rbc2)
    ra2=ra*ra
    ra5=ra2*ra2*ra

    lds233(0,i)=-2.0_8*Ct3*z3a*(3.0_8*z2a2-z3a2)/(3.0_8*pi*rbc6)*costa
    lds233(1,i)=Ct3*z1a*z3a3/(pi*ra5*rbc2)*costa
    tmp=Ct*z3a/(pi*rbc2)*costa
    lds233(2,i)=tmp*(1.0_8-2.0_8*p2+2.0_8*p2*z3a2/rbc2)
    lds233(3,i)=tmp*(1.0_8-2.0_8*z3a2/rbc2)
  end subroutine Ldyns233

  subroutine Ldyns311(z1a,z2a,z3a,costa,sinta,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,costa,sinta
    integer,intent(in) :: i
    real(8) :: rbc,ra 
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2
    real(8) :: tmp1,tmp2
    
    z2a2=z2a*z2a
    z3a2=z3a*z3a
    rbc=sqrt(z2a2+z3a2)
    rbc2=rbc*rbc
    rbc6=rbc2*rbc2*rbc2
    ra=sqrt(z1a*z1a+rbc2)
    ra2=ra*ra
    ra5=ra2*ra2*ra
    sinta2=sinta*sinta
    costa2=costa*costa

    lds311(0,i)=2.0_8*Ct3*z2a*(3.0_8*z3a2-z2a2)/(3.0_8*pi*rbc6)*sinta2
    lds311(1,i)=Ct3/(pi*ra5)*sinta*(z1a*z2a*z3a2/rbc2*sinta+(z3a2-ra2/3.0_8 )*costa)
    tmp1=Ct*z2a/(pi*rbc2)
    lds311(2,i)=tmp1*(2.0_8*p2*z3a2/rbc2*sinta2)
    lds311(3,i)=tmp1*(-(0.5_8*costa2+2.0_8*z3a2/rbc2*sinta2))
    tmp2=Ct/(pi*ra)*sinta*costa
    lds311(4,i)=tmp2*(2.0_8*p2/3.0_8)
    lds311(5,i)=tmp2*(-1.0_8/6.0_8)
    lds311(6,i)=0.5_8/pi
  end subroutine Ldyns311

  subroutine Ldyns312(z1a,z2a,z3a,costa,sinta,i) !derive from Ldyns321
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,costa,sinta
    integer,intent(in) :: i
    real(8) :: rbc,ra
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2
    real(8) :: tmp1,tmp2

    z2a2=z2a*z2a
    z3a2=z3a*z3a
    rbc=sqrt(z2a2+z3a2)
    rbc2=rbc*rbc
    rbc6=rbc2*rbc2*rbc2
    ra=sqrt(z1a*z1a+rbc2)
    ra2=ra*ra
    ra5=ra2*ra2*ra
    sinta2=sinta*sinta
    costa2=costa*costa

    lds312(0,i)=-2.0_8*Ct3*z2a*(3.0_8*z3a2-z2a2)/(3.0_8*pi*rbc6)*sinta*costa
    lds312(1,i)=-Ct3/(pi*ra5)*costa*(z1a*z2a*z3a2/rbc2*sinta+(z3a2-ra2/3.0_8)*costa)
    tmp1=-Ct*z2a/(pi*rbc2)*sinta*costa
    lds312(2,i)=tmp1*(2.0_8*p2*z3a2/rbc2)
    lds312(3,i)=tmp1*(0.5_8-2.0_8*z3a2/rbc2)
    tmp2=-Ct/(pi*ra)*costa2
    lds312(4,i)=tmp2*(2.0_8*p2/3.0_8)
    lds312(5,i)=tmp2*(-1.0_8/6.0_8)
  end subroutine Ldyns312

  subroutine Ldyns313(z1a,z2a,z3a,costa,sinta,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,costa,sinta
    integer,intent(in) :: i
    real(8) :: rbc,ra 
    real(8) :: z2a2,z3a2,z3a3,rbc2,rbc6,ra2,ra5
    real(8) :: tmp

    z2a2=z2a*z2a
    z3a2=z3a*z3a
    z3a3=z3a*z3a2
    rbc=sqrt(z2a2+z3a2)
    rbc2=rbc*rbc
    rbc6=rbc2*rbc2*rbc2
    ra=sqrt(z1a*z1a+rbc2)
    ra2=ra*ra
    ra5=ra2*ra2*ra

    lds313(0,i)=2.0_8*Ct3*z3a*(3.0_8*z2a2-z3a2)/(3.0_8*pi*rbc6)*sinta
    lds313(1,i)=-Ct3*z1a*z3a3/(pi*ra5*rbc2)*sinta
    tmp=-Ct*z3a/(pi*rbc2)*sinta
    lds313(2,i)=tmp*(1.0_8-2.0_8*p2+2.0_8*p2*z3a2/rbc2)
    lds313(3,i)=tmp*(1.0_8-2.0_8*z3a2/rbc2)
  end subroutine Ldyns313

  subroutine Ldyns321(z1a,z2a,z3a,costa,sinta,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,costa,sinta
    integer,intent(in) :: i
    real(8) :: rbc,ra
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2
    real(8) :: tmp1,tmp2

    z2a2=z2a*z2a
    z3a2=z3a*z3a
    rbc=sqrt(z2a2+z3a2)
    rbc2=rbc*rbc
    rbc6=rbc2*rbc2*rbc2
    ra=sqrt(z1a*z1a+rbc2)
    ra2=ra*ra
    ra5=ra2*ra2*ra
    sinta2=sinta*sinta
    costa2=costa*costa

    lds321(0,i)=-2.0_8*Ct3*z2a*(3.0_8*z3a2-z2a2)/(3.0_8*pi*rbc6)*sinta*costa
    lds321(1,i)=-Ct3/(pi*ra5)*sinta*(z1a*z2a*z3a2/rbc2*costa-(z3a2-ra2/3.0_8)*sinta)
    tmp1=-Ct*z2a/(pi*rbc2)*sinta*costa
    lds321(2,i)=tmp1*(2.0_8*p2*z3a2/rbc2)
    lds321(3,i)=tmp1*(0.5_8-2.0_8*z3a2/rbc2)
    tmp2=Ct/(pi*ra)*sinta2
    lds321(4,i)=tmp2*(2.0_8*p2/3.0_8)
    lds321(5,i)=tmp2*(-1.0_8/6.0_8)
  end subroutine Ldyns321
  
  subroutine Ldyns322(z1a,z2a,z3a,costa,sinta,i) ! derive from Ldyns311
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,costa,sinta
    integer,intent(in) :: i
    real(8) :: rbc,ra
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2
    real(8) :: tmp1,tmp2
    
    z2a2=z2a*z2a
    z3a2=z3a*z3a
    rbc=sqrt(z2a2+z3a2)
    rbc2=rbc*rbc
    rbc6=rbc2*rbc2*rbc2
    ra=sqrt(z1a*z1a+rbc2)
    ra2=ra*ra
    ra5=ra2*ra2*ra
    sinta2=sinta*sinta
    costa2=costa*costa

    lds322(0,i)=2.0_8*Ct3*z2a*(3.0_8*z3a2-z2a2)/(3.0_8*pi*rbc6)*costa2
    lds322(1,i)=-Ct3/(pi*ra5)*costa*(-z1a*z2a*z3a2/rbc2*costa+(z3a2-ra2/3.0_8)*sinta)
    tmp1=Ct*z2a/(pi*rbc2)
    lds322(2,i)=tmp1*(2.0_8*p2*z3a2/rbc2*costa2)
    lds322(3,i)=tmp1*(-(0.5_8*sinta2+2.0_8*z3a2/rbc2*costa2))
    tmp2=-Ct/(pi*ra)*sinta*costa
    lds322(4,i)=tmp2*(2.0_8*p2/3.0_8)
    lds322(5,i)=tmp2*(-1.0_8/6.0_8)
    lds322(6,i)=0.5_8/pi
  end subroutine Ldyns322

  subroutine Ldyns323(z1a,z2a,z3a,costa,sinta,i) !derive from Ldyns313
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,costa,sinta
    integer,intent(in) :: i
    real(8) :: rbc,ra 
    real(8) :: z2a2,z3a2,z3a3,rbc2,rbc6,ra2,ra5
    real(8) :: tmp

    z2a2=z2a*z2a
    z3a2=z3a*z3a
    z3a3=z3a*z3a2
    rbc=sqrt(z2a2+z3a2)
    rbc2=rbc*rbc
    rbc6=rbc2*rbc2*rbc2
    ra=sqrt(z1a*z1a+rbc2)
    ra2=ra*ra
    ra5=ra2*ra2*ra

    lds323(0,i)=-2.0_8*Ct3*z3a*(3.0_8*z2a2-z3a2)/(3.0_8*pi*rbc6)*costa
    lds323(1,i)=Ct3*z1a*z3a3/(pi*ra5*rbc2)*costa
    tmp=Ct*z3a/(pi*rbc2)*costa
    lds323(2,i)=tmp*(1.0_8-2.0_8*p2+2.0_8*p2*z3a2/rbc2)
    lds323(3,i)=tmp*(1.0_8-2.0_8*z3a2/rbc2)
  end subroutine Ldyns323

  subroutine Ldyns331(z1a,z2a,z3a,costa,sinta,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,costa,sinta
    integer,intent(in) :: i
    real(8) :: rbc,ra 
    real(8) :: z2a2,z3a2,z3a3,rbc2,rbc6,ra2,ra5
    real(8) :: tmp

    z2a2=z2a*z2a
    z3a2=z3a*z3a
    z3a3=z3a*z3a2
    rbc=sqrt(z2a2+z3a2)
    rbc2=rbc*rbc
    rbc6=rbc2*rbc2*rbc2
    ra=sqrt(z1a*z1a+rbc2)
    ra2=ra*ra
    ra5=ra2*ra2*ra

    lds331(0,i)=2.0_8*Ct3*z3a*(3.0_8*z2a2-z3a2)/(3.0_8*pi*rbc6)*sinta
    lds331(1,i)=-Ct3*z1a*z3a3/(pi*ra5*rbc2)*sinta
    tmp=-Ct*z3a/(pi*rbc2)*sinta
    lds331(2,i)=tmp*(1.0_8-2.0_8*p2+2.0_8*p2*z3a2/rbc2)
    lds331(3,i)=tmp*(1.0_8-2.0_8*z3a2/rbc2)
  end subroutine Ldyns331

  subroutine Ldyns332(z1a,z2a,z3a,costa,sinta,i) !derive from Ldyns331
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,costa,sinta
    integer,intent(in) :: i
    real(8) :: rbc,ra
    real(8):: z2a2,z3a2,z3a3,rbc2,rbc6,ra2,ra5
    real(8) :: tmp

    z2a2=z2a*z2a
    z3a2=z3a*z3a
    z3a3=z3a*z3a2
    rbc=sqrt(z2a2+z3a2)
    rbc2=rbc*rbc
    rbc6=rbc2*rbc2*rbc2
    ra=sqrt(z1a*z1a+rbc2)
    ra2=ra*ra
    ra5=ra2*ra2*ra

    lds332(0,i)=-2.0_8*Ct3*z3a*(3.0_8*z2a2-z3a2)/(3.0_8*pi*rbc6)*costa
    lds332(1,i)=Ct3*z1a*z3a3/(pi*ra5*rbc2)*costa
    tmp=Ct*z3a/(pi*rbc2)*costa
    lds332(2,i)=tmp*(1.0_8-2.0_8*p2+2.0_8*p2*z3a2/rbc2)
    lds332(3,i)=tmp*(1.0_8-2.0_8*z3a2/rbc2)
  end subroutine Ldyns332

  subroutine Ldyns333(z1a,z2a,z3a,costa,sinta,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,costa,sinta ! actually, cos & sin are not needed
    integer,intent(in) :: i
    real(8) :: rbc,ra
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5
    real(8) :: tmp

    z2a2=z2a*z2a
    z3a2=z3a*z3a
    rbc=sqrt(z2a2+z3a2)
    rbc2=rbc*rbc
    rbc6=rbc2*rbc2*rbc2
    ra= sqrt(z1a*z1a+rbc2)
    ra2=ra*ra
    ra5=ra2*ra2*ra

    lds333(0,i)=-2.0_8*Ct3*z2a*(3.0_8*z3a2-z2a2)/(3.0_8*pi*rbc6)
    lds333(1,i)=-Ct3*z1a*z2a*z3a2/(pi*ra5*rbc2)
    tmp=-Ct*z2a/(pi*rbc2)
    lds333(2,i)=tmp*(2.0_8*(1-p2)+2.0_8*p2*z3a2/rbc2)
    lds333(3,i)=tmp*(-2.0_8*z3a2/rbc2)
    lds333(4,i)=1/(2.0_8*pi*p)
  end subroutine Ldyns333

  ! (function Ldynf: multiplications of the g and spatial dependent coefficient to calculate L0)
  real(8) function Ldynf111(z1a,z2a,z3a,t,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,t
    integer,intent(in) :: i
    real(8) :: rbc,ra 
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2,sinta3,g1t,g1l,g2t,g2l,g4t,g4l
    real(8) :: L111

    rbc=sqrt(z2a*z2a+z3a*z3a)
    rbc2=rbc*rbc
    ra=sqrt(z1a*z1a+rbc2)
    call Calg142(z1a,z2a,t,rbc,ra,g1l,g2l,g4l,g1t,g2t,g4t)
    L111 = lds111(0,i)*( g1l - g1t )&
         + lds111(1,i)*( g4l - g4t )&
         + lds111(2,i)*g2l&
         + lds111(3,i)*g2t
    Ldynf111 = L111
  end function Ldynf111

  real(8) function Ldynf112(z1a,z2a,z3a,t,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,t
    integer,intent(in) :: i
    real(8) :: rbc,ra 
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2,g1t,g1l,g2t,g2l,g4t,g4l 
    real(8) :: L112

    rbc=sqrt(z2a*z2a+z3a*z3a)
    rbc2=rbc*rbc
    ra=sqrt(z1a*z1a+rbc2)
    call Calg142(z1a,z2a,t,rbc,ra,g1l,g2l,g4l,g1t,g2t,g4t)
    L112 = lds112(0,i)*( g1l - g1t )&
         + lds112(1,i)*( g4l - g4t )&
         + lds112(2,i)*g2l&
         + lds112(3,i)*g2t
    Ldynf112 = L112
  end function Ldynf112

  real(8) function Ldynf113(z1a,z2a,z3a,t,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,t
    integer,intent(in) :: i
    real(8)::rbc,ra 
    real(8)::z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2,g1t,g1l,g2t,g2l,g3t,g3l,g4t,g4l,g5t,g5l
    real(8)::L113

    rbc=sqrt(z2a*z2a+z3a*z3a)
    rbc2=rbc*rbc
    ra=sqrt(z1a*z1a+rbc2)
    call Calg142(z1a,z2a,t,rbc,ra,g1l,g2l,g4l,g1t,g2t,g4t)
    call Calg5(z1a,z2a,t,rbc,ra,g5l,g5t)
    call Calg3(z1a,z2a,t,rbc,ra,g3l,g3t)
    L113 = lds113(0,i)*( g1l - g1t )&
         + lds113(1,i)*( g4l - g4t )&
         + lds113(2,i)*g2l&
         + lds113(3,i)*g2t&
         + lds113(4,i)*g5l&
         + lds113(5,i)*g5t&
         + lds113(6,i)*g3l
    Ldynf113 = L113
  end function Ldynf113

  real(8) function Ldynf121(z1a,z2a,z3a,t,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,t
    integer,intent(in) :: i
    real(8) :: rbc,ra 
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2,g1t,g1l,g2t,g2l,g4t,g4l
    real(8) :: L121

    rbc=sqrt(z2a*z2a+z3a*z3a)
    rbc2=rbc*rbc
    ra=sqrt(z1a*z1a+rbc2)
    call Calg142(z1a,z2a,t,rbc,ra,g1l,g2l,g4l,g1t,g2t,g4t)
    L121 = lds121(0,i)*( g1l - g1t )&
         + lds121(1,i)*( g4l - g4t )&
         + lds121(2,i)*g2l&
         + lds121(3,i)*g2t
    Ldynf121 = L121
  end function Ldynf121

  real(8) function Ldynf122(z1a,z2a,z3a,t,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,t
    integer,intent(in) :: i
    real(8) :: rbc,ra
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2,g1t,g1l,g2t,g2l,g4t,g4l
    real(8) :: L122

    rbc=sqrt(z2a*z2a+z3a*z3a)
    rbc2=rbc*rbc
    ra=sqrt(z1a*z1a+rbc2)
    call Calg142(z1a,z2a,t,rbc,ra,g1l,g2l,g4l,g1t,g2t,g4t)
    L122 = lds122(0,i)*( g1l - g1t )&
         + lds122(1,i)*( g4l - g4t )&
         + lds122(2,i)*g2l&
         + lds122(3,i)*g2t
    Ldynf122 = L122
  end function Ldynf122

  real(8) function Ldynf123(z1a,z2a,z3a,t,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,t
    integer,intent(in) :: i
    real(8) :: rbc,ra
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2,g1t,g1l,g2t,g2l,g4t,g4l,g5t,g5l
    real(8) :: L123
      
    rbc=sqrt(z2a*z2a+z3a*z3a)
    rbc2=rbc*rbc
    ra=sqrt(z1a*z1a+rbc2)
    call Calg142(z1a,z2a,t,rbc,ra,g1l,g2l,g4l,g1t,g2t,g4t)
    call Calg5(z1a,z2a,t,rbc,ra,g5l,g5t)
    L123 = lds123(0,i)*( g1l - g1t )&
         + lds123(1,i)*( g4l - g4t )&
         + lds123(2,i)*g2l&
         + lds123(3,i)*g2t&
         + lds123(4,i)*g5l&
         + lds123(5,i)*g5t
    Ldynf123 = L123
  end function Ldynf123

  real(8) function Ldynf131(z1a,z2a,z3a,t,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,t
    integer,intent(in) :: i
    real(8) :: rbc,ra
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2,g1t,g1l,g2t,g2l,g3t,g3l,g4t,g4l,g5t,g5l
    real(8) :: L131

    rbc=sqrt(z2a*z2a+z3a*z3a)
    rbc2=rbc*rbc
    ra=sqrt(z1a*z1a+rbc2)
    call Calg142(z1a,z2a,t,rbc,ra,g1l,g2l,g4l,g1t,g2t,g4t)
    call Calg5(z1a,z2a,t,rbc,ra,g5l,g5t)
    call Calg3(z1a,z2a,t,rbc,ra,g3l,g3t)
    L131 = lds131(0,i)*( g1l - g1t )&
         + lds131(1,i)*( g4l - g4t )&
         + lds131(2,i)*g2l&
         + lds131(3,i)*g2t&
         + lds131(4,i)*g5l&
         + lds131(5,i)*g5t&
         + lds131(6,i)*g3t
    Ldynf131 = L131
  end function Ldynf131

  real(8) function Ldynf132(z1a,z2a,z3a,t,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,t
    integer,intent(in) :: i
    real(8) :: rbc,ra 
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2,g1t,g1l,g2t,g2l,g4t,g4l,g5t,g5l
    real(8) :: L132
      
    rbc=sqrt(z2a*z2a+z3a*z3a)
    rbc2=rbc*rbc
    ra=sqrt(z1a*z1a+rbc2)
    call Calg142(z1a,z2a,t,rbc,ra,g1l,g2l,g4l,g1t,g2t,g4t)
    call Calg5(z1a,z2a,t,rbc,ra,g5l,g5t)
    L132 = lds132(0,i)*( g1l - g1t )&
         + lds132(1,i)*( g4l - g4t )&
         + lds132(2,i)*g2l&
         + lds132(3,i)*g2t&
         + lds132(4,i)*g5l&
         + lds132(5,i)*g5t
    Ldynf132 = L132
  end function Ldynf132

  real(8) function Ldynf133(z1a,z2a,z3a,t,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,t
    integer,intent(in) :: i
    real(8) :: rbc,ra
    real(8) :: z2a2,z3a2,z3a3,rbc2,rbc6,ra2,ra5,sinta2,costa2,g1t,g1l,g2t,g2l,g4t,g4l
    real(8) :: L133

    rbc=sqrt(z2a*z2a+z3a*z3a)
    rbc2=rbc*rbc
    ra=sqrt(z1a*z1a+rbc2)
    call Calg142(z1a,z2a,t,rbc,ra,g1l,g2l,g4l,g1t,g2t,g4t)
    L133 = lds133(0,i)*( g1l - g1t )&
         + lds133(1,i)*( g4l - g4t )&
         + lds133(2,i)*g2l&
         + lds133(3,i)*g2t
    Ldynf133 = L133
  end function Ldynf133

  real(8) function Ldynf211(z1a,z2a,z3a,t,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,t
    integer,intent(in) :: i
    real(8) :: rbc,ra 
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2,g1t,g1l,g2t,g2l,g4t,g4l
    real(8) :: L211

    rbc=sqrt(z2a*z2a+z3a*z3a)
    rbc2=rbc*rbc
    ra=sqrt(z1a*z1a+rbc2)
    call Calg142(z1a,z2a,t,rbc,ra,g1l,g2l,g4l,g1t,g2t,g4t)
    L211 = lds211(0,i)*( g1l - g1t )&
         + lds211(1,i)*( g4l - g4t )&
         + lds211(2,i)*g2l&
         + lds211(3,i)*g2t
    Ldynf211 = L211
  end function Ldynf211

  real(8) function Ldynf212(z1a,z2a,z3a,t,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,t
    integer,intent(in) :: i
    real(8) :: rbc,ra
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2,g1t,g1l,g2t,g2l,g4t,g4l
    real(8) :: L212

    rbc=sqrt(z2a*z2a+z3a*z3a)
    rbc2=rbc*rbc
    ra=sqrt(z1a*z1a+rbc2)
    call Calg142(z1a,z2a,t,rbc,ra,g1l,g2l,g4l,g1t,g2t,g4t)
    L212 = lds212(0,i)*( g1l - g1t )&
         + lds212(1,i)*( g4l - g4t )&
         + lds212(2,i)*g2l&
         + lds212(3,i)*g2t
    Ldynf212 = L212
  end function Ldynf212

  real(8) function Ldynf213(z1a,z2a,z3a,t,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,t
    integer,intent(in) :: i
    real(8) :: rbc,ra
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2,g1t,g1l,g2t,g2l,g4t,g4l,g5t,g5l
    real(8) :: L213
      
    rbc=sqrt(z2a*z2a+z3a*z3a)
    rbc2=rbc*rbc
    ra=sqrt(z1a*z1a+rbc2)
    call Calg142(z1a,z2a,t,rbc,ra,g1l,g2l,g4l,g1t,g2t,g4t)
    call Calg5(z1a,z2a,t,rbc,ra,g5l,g5t)
    L213 = lds213(0,i)*( g1l - g1t )&
         + lds213(1,i)*( g4l - g4t )&
         + lds213(2,i)*g2l&
         + lds213(3,i)*g2t&
         + lds213(4,i)*g5l&
         + lds213(5,i)*g5t
    Ldynf213 = L213
  end function Ldynf213
  
  real(8) function Ldynf221(z1a,z2a,z3a,t,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,t
    integer,intent(in) :: i
    real(8) :: rbc,ra 
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2,g1t,g1l,g2t,g2l,g4t,g4l 
    real(8) :: L221

    rbc=sqrt(z2a*z2a+z3a*z3a)
    rbc2=rbc*rbc
    ra=sqrt(z1a*z1a+rbc2)
    call Calg142(z1a,z2a,t,rbc,ra,g1l,g2l,g4l,g1t,g2t,g4t)
    L221 = lds221(0,i)*( g1l - g1t )&
         + lds221(1,i)*( g4l - g4t )&
         + lds221(2,i)*g2l&
         + lds221(3,i)*g2t
    Ldynf221 = L221
  end function Ldynf221

  real(8) function Ldynf222(z1a,z2a,z3a,t,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,t
    integer,intent(in) :: i
    real(8) :: rbc,ra
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2,sinta3,g1t,g1l,g2t,g2l,g4t,g4l
    real(8) :: L222

    rbc=sqrt(z2a*z2a+z3a*z3a)
    rbc2=rbc*rbc
    ra=sqrt(z1a*z1a+rbc2)
    call Calg142(z1a,z2a,t,rbc,ra,g1l,g2l,g4l,g1t,g2t,g4t)

    L222 = lds222(0,i)*( g1l - g1t )&
         + lds222(1,i)*( g4l - g4t )&
         + lds222(2,i)*g2l&
         + lds222(3,i)*g2t
    Ldynf222 = L222
  end function Ldynf222

  real(8) function Ldynf223(z1a,z2a,z3a,t,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,t
    integer,intent(in) :: i
    real(8) :: rbc,ra 
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2,g1t,g1l,g2t,g2l,g3t,g3l,g4t,g4l,g5t,g5l
    real(8) :: L223

    rbc=sqrt(z2a*z2a+z3a*z3a)
    rbc2=rbc*rbc
    ra=sqrt(z1a*z1a+rbc2)
    call Calg142(z1a,z2a,t,rbc,ra,g1l,g2l,g4l,g1t,g2t,g4t)
    call Calg5(z1a,z2a,t,rbc,ra,g5l,g5t)
    call Calg3(z1a,z2a,t,rbc,ra,g3l,g3t)
    L223 = lds223(0,i)*( g1l - g1t )&
         + lds223(1,i)*( g4l - g4t )&
         + lds223(2,i)*g2l&
         + lds223(3,i)*g2t&
         + lds223(4,i)*g5l&
         + lds223(5,i)*g5t&
         + lds223(6,i)*g3l
    Ldynf223 = L223
  end function Ldynf223

  real(8) function Ldynf231(z1a,z2a,z3a,t,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,t
    integer,intent(in) :: i
    real(8) :: rbc,ra 
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2,g1t,g1l,g2t,g2l,g4t,g4l,g5t,g5l
    real(8) :: L231
   
    rbc=sqrt(z2a*z2a+z3a*z3a)
    rbc2=rbc*rbc
    ra=sqrt(z1a*z1a+rbc2)
    call Calg142(z1a,z2a,t,rbc,ra,g1l,g2l,g4l,g1t,g2t,g4t)
    call Calg5(z1a,z2a,t,rbc,ra,g5l,g5t)
    L231 = lds231(0,i)*( g1l - g1t )&
         + lds231(1,i)*( g4l - g4t )&
         + lds231(2,i)*g2l&
         + lds231(3,i)*g2t&
         + lds231(4,i)*g5l&
         + lds231(5,i)*g5t
    Ldynf231 = L231
  end function Ldynf231

  real(8) function Ldynf232(z1a,z2a,z3a,t,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,t
    integer,intent(in) :: i
    real(8) :: rbc,ra 
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2,g1t,g1l,g2t,g2l,g3t,g3l,g4t,g4l,g5t,g5l
    real(8) :: L232

    rbc=sqrt(z2a*z2a+z3a*z3a)
    rbc2=rbc*rbc
    ra=sqrt(z1a*z1a+rbc2)
    call Calg142(z1a,z2a,t,rbc,ra,g1l,g2l,g4l,g1t,g2t,g4t)
    call Calg5(z1a,z2a,t,rbc,ra,g5l,g5t)
    call Calg3(z1a,z2a,t,rbc,ra,g3l,g3t)
    L232 = lds232(0,i)*( g1l - g1t )&
         + lds232(1,i)*( g4l - g4t )&
         + lds232(2,i)*g2l&
         + lds232(3,i)*g2t&
         + lds232(4,i)*g5l&
         + lds232(5,i)*g5t&
         + lds232(6,i)*g3t
    Ldynf232 = L232
  end function Ldynf232

  real(8) function Ldynf233(z1a,z2a,z3a,t,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,t
    integer,intent(in) :: i
    real(8) :: rbc,ra 
    real(8) :: z2a2,z3a2,z3a3,rbc2,rbc6,ra2,ra5,sinta2,costa2,g1t,g1l,g2t,g2l,g4t,g4l
    real(8) :: L233

    rbc=sqrt(z2a*z2a+z3a*z3a)
    rbc2=rbc*rbc
    ra=sqrt(z1a*z1a+rbc2)
    call Calg142(z1a,z2a,t,rbc,ra,g1l,g2l,g4l,g1t,g2t,g4t)
    L233 = lds233(0,i)*( g1l - g1t )&
         + lds233(1,i)*( g4l - g4t )&
         + lds233(2,i)*g2l&
         + lds233(3,i)*g2t
    Ldynf233 = L233
  end function Ldynf233

  real(8) function Ldynf311(z1a,z2a,z3a,t,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,t
    integer,intent(in) :: i
    real(8) :: rbc,ra
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2,g1t,g1l,g2t,g2l,g3t,g3l,g4t,g4l,g5t,g5l
    real(8) :: L311

    rbc=sqrt(z2a*z2a+z3a*z3a)
    rbc2=rbc*rbc
    ra=sqrt(z1a*z1a+rbc2)
    call Calg142(z1a,z2a,t,rbc,ra,g1l,g2l,g4l,g1t,g2t,g4t)
    call Calg5(z1a,z2a,t,rbc,ra,g5l,g5t)
    call Calg3(z1a,z2a,t,rbc,ra,g3l,g3t)
    L311 = lds311(0,i)*( g1l - g1t )&
         + lds311(1,i)*( g4l - g4t )&
         + lds311(2,i)*g2l&
         + lds311(3,i)*g2t&
         + lds311(4,i)*g5l&
         + lds311(5,i)*g5t&
         + lds311(6,i)*g3t
    Ldynf311 = L311
  end function Ldynf311

  real(8) function Ldynf312(z1a,z2a,z3a,t,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,t
    integer,intent(in) :: i
    real(8) :: rbc,ra 
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2,g1t,g1l,g2t,g2l,g4t,g4l,g5t,g5l
    real(8) :: L312
      
    rbc=sqrt(z2a*z2a+z3a*z3a)
    rbc2=rbc*rbc
    ra=sqrt(z1a*z1a+rbc2)
    call Calg142(z1a,z2a,t,rbc,ra,g1l,g2l,g4l,g1t,g2t,g4t)
    call Calg5(z1a,z2a,t,rbc,ra,g5l,g5t)
    L312 = lds312(0,i)*( g1l - g1t )&
         + lds312(1,i)*( g4l - g4t )&
         + lds312(2,i)*g2l&
         + lds312(3,i)*g2t&
         + lds312(4,i)*g5l&
         + lds312(5,i)*g5t
    Ldynf312 = L312
  end function Ldynf312

  real(8) function Ldynf313(z1a,z2a,z3a,t,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,t
    integer,intent(in) :: i
    real(8) :: rbc,ra
    real(8) :: z2a2,z3a2,z3a3,rbc2,rbc6,ra2,ra5,sinta2,costa2,g1t,g1l,g2t,g2l,g4t,g4l
    real(8) :: L313

    rbc=sqrt(z2a*z2a+z3a*z3a)
    rbc2=rbc*rbc
    ra=sqrt(z1a*z1a+rbc2)
    call Calg142(z1a,z2a,t,rbc,ra,g1l,g2l,g4l,g1t,g2t,g4t)
    L313 = lds313(0,i)*( g1l - g1t )&
         + lds313(1,i)*( g4l - g4t )&
         + lds313(2,i)*g2l&
         + lds313(3,i)*g2t
    Ldynf313 = L313
  end function Ldynf313

  real(8) function Ldynf321(z1a,z2a,z3a,t,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,t
    integer,intent(in) :: i
    real(8) :: rbc,ra 
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2,g1t,g1l,g2t,g2l,g4t,g4l,g5t,g5l
    real(8) :: L321
   
    rbc=sqrt(z2a*z2a+z3a*z3a)
    rbc2=rbc*rbc
    ra=sqrt(z1a*z1a+rbc2)
    call Calg142(z1a,z2a,t,rbc,ra,g1l,g2l,g4l,g1t,g2t,g4t)
    call Calg5(z1a,z2a,t,rbc,ra,g5l,g5t)
    L321 = lds321(0,i)*( g1l - g1t )&
         + lds321(1,i)*( g4l - g4t )&
         + lds321(2,i)*g2l&
         + lds321(3,i)*g2t&
         + lds321(4,i)*g5l&
         + lds321(5,i)*g5t
    Ldynf321 = L321
  end function Ldynf321

  real(8) function Ldynf322(z1a,z2a,z3a,t,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,t
    integer,intent(in) :: i
    real(8) :: rbc,ra 
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2,g1t,g1l,g2t,g2l,g3t,g3l,g4t,g4l,g5t,g5l
    real(8) :: L322

    rbc=sqrt(z2a*z2a+z3a*z3a)
    rbc2=rbc*rbc
    ra=sqrt(z1a*z1a+rbc2)
    call Calg142(z1a,z2a,t,rbc,ra,g1l,g2l,g4l,g1t,g2t,g4t)
    call Calg5(z1a,z2a,t,rbc,ra,g5l,g5t)
    call Calg3(z1a,z2a,t,rbc,ra,g3l,g3t)
    L322 = lds322(0,i)*( g1l - g1t )&
         + lds322(1,i)*( g4l - g4t )&
         + lds322(2,i)*g2l&
         + lds322(3,i)*g2t&
         + lds322(4,i)*g5l&
         + lds322(5,i)*g5t&
         + lds322(6,i)*g3t
    Ldynf322 = L322
  end function Ldynf322

  real(8) function Ldynf323(z1a,z2a,z3a,t,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,t
    integer,intent(in) :: i
    real(8) :: rbc,ra 
    real(8) :: z2a2,z3a2,z3a3,rbc2,rbc6,ra2,ra5,sinta2,costa2,g1t,g1l,g2t,g2l,g4t,g4l
    real(8) :: L323

    rbc=sqrt(z2a*z2a+z3a*z3a)
    rbc2=rbc*rbc
    ra=sqrt(z1a*z1a+rbc2)
    call Calg142(z1a,z2a,t,rbc,ra,g1l,g2l,g4l,g1t,g2t,g4t)
    L323 = lds323(0,i)*( g1l - g1t )&
         + lds323(1,i)*( g4l - g4t )&
         + lds323(2,i)*g2l&
         + lds323(3,i)*g2t
    Ldynf323 = L323
  end function Ldynf323

  real(8) function Ldynf331(z1a,z2a,z3a,t,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,t
    integer,intent(in) :: i
    real(8) :: rbc,ra 
    real(8) :: z2a2,z3a2,z3a3,rbc2,rbc6,ra2,ra5,sinta2,costa2,g1t,g1l,g2t,g2l,g4t,g4l
    real(8) :: L331

    rbc=sqrt(z2a*z2a+z3a*z3a)
    rbc2=rbc*rbc
    ra=sqrt(z1a*z1a+rbc2)
    call Calg142(z1a,z2a,t,rbc,ra,g1l,g2l,g4l,g1t,g2t,g4t)
    L331 = lds331(0,i)*( g1l - g1t )&
         + lds331(1,i)*( g4l - g4t )&
         + lds331(2,i)*g2l&
         + lds331(3,i)*g2t
    Ldynf331 = L331
  end function Ldynf331

  real(8) function Ldynf332(z1a,z2a,z3a,t,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,t
    integer,intent(in) :: i
    real(8) :: rbc,ra 
    real(8) :: z2a2,z3a2,z3a3,rbc2,rbc6,ra2,ra5,sinta2,costa2,g1t,g1l,g2t,g2l,g4t,g4l
    real(8) :: L332

    rbc=sqrt(z2a*z2a+z3a*z3a)
    rbc2=rbc*rbc
    ra=sqrt(z1a*z1a+rbc2)
    call Calg142(z1a,z2a,t,rbc,ra,g1l,g2l,g4l,g1t,g2t,g4t)
    L332 = lds332(0,i)*( g1l - g1t )&
         + lds332(1,i)*( g4l - g4t )&
         + lds332(2,i)*g2l&
         + lds332(3,i)*g2t
    Ldynf332 = L332
  end function Ldynf332

  real(8) function Ldynf333(z1a,z2a,z3a,t,i)
    implicit none
    real(8),intent(in) :: z1a,z2a,z3a,t
    integer,intent(in) :: i
    real(8) :: rbc,ra
    real(8) :: z2a2,z3a2,rbc2,rbc6,ra2,ra5,sinta2,costa2,g1t,g1l,g2t,g2l,g3t,g3l,g4t,g4l
    real(8) :: L333
      
    rbc=sqrt(z2a*z2a+z3a*z3a)
    rbc2=rbc*rbc
    ra=sqrt(z1a*z1a+rbc2)
    call Calg142(z1a,z2a,t,rbc,ra,g1l,g2l,g4l,g1t,g2t,g4t)
    call Calg3(z1a,z2a,t,rbc,ra,g3l,g3t)
    L333 = lds333(0,i)*( g1l - g1t )&
         + lds333(1,i)*( g4l - g4t )&
         + lds333(2,i)*g2l&
         + lds333(3,i)*g2t&
         + lds333(4,i)*g3l
    Ldynf333 = L333
  end function Ldynf333

  !(subroutine Kdyns: compute six spatial-dependent coefficients to obtain six L0 which make up L)
  !( 27 subroutine are exactly the same code except for the subscripts)
  
  subroutine Kdyns111(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3
    real(8) :: z1a,zi1,zii1,z2,z3,costa,sinta
 
    call Rotcoordnt(xr1,xr2,xr3,xsb1,xsb2,xsc1,xsc2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(0) = zii1 ; zi1s(0) = zi1 ; z2s(0) = z2 ; z3s(0) = z3
    call Ldyns111(zii1,z2,z3,costa,sinta,0)
    call Ldyns111(zi1,z2,z3,costa,sinta,1)
 
    call Rotcoordnt(xr1,xr2,xr3,xsc1,xsc2,xsa1,xsa2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(1) = zii1 ; zi1s(1) = zi1 ; z2s(1) = z2 ; z3s(1) = z3
    call Ldyns111(zii1,z2,z3,costa,sinta,2)
    call Ldyns111(zi1,z2,z3,costa,sinta,3)
 
    call Rotcoordnt(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(2) = zii1 ; zi1s(2) = zi1 ; z2s(2) = z2 ; z3s(2) = z3
    call Ldyns111(zii1,z2,z3,costa,sinta,4) 
    call Ldyns111(zi1,z2,z3,costa,sinta,5) 
  end subroutine Kdyns111

  subroutine Kdyns112(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3
    real(8) :: z1a,zi1,zii1,z2,z3,costa,sinta
 
    call Rotcoordnt(xr1,xr2,xr3,xsb1,xsb2,xsc1,xsc2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(0) = zii1 ; zi1s(0) = zi1 ; z2s(0) = z2 ; z3s(0) = z3
    call Ldyns112(zii1,z2,z3,costa,sinta,0)
    call Ldyns112(zi1,z2,z3,costa,sinta,1)
 
    call Rotcoordnt(xr1,xr2,xr3,xsc1,xsc2,xsa1,xsa2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(1) = zii1 ; zi1s(1) = zi1 ; z2s(1) = z2 ; z3s(1) = z3
    call Ldyns112(zii1,z2,z3,costa,sinta,2)
    call Ldyns112(zi1,z2,z3,costa,sinta,3)
 
    call Rotcoordnt(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(2) = zii1 ; zi1s(2) = zi1 ; z2s(2) = z2 ; z3s(2) = z3
    call Ldyns112(zii1,z2,z3,costa,sinta,4) 
    call Ldyns112(zi1,z2,z3,costa,sinta,5) 
  end subroutine Kdyns112

  subroutine Kdyns113(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3
    real(8) :: z1a,zi1,zii1,z2,z3,costa,sinta
 
    call Rotcoordnt(xr1,xr2,xr3,xsb1,xsb2,xsc1,xsc2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(0) = zii1 ; zi1s(0) = zi1 ; z2s(0) = z2 ; z3s(0) = z3
    call Ldyns113(zii1,z2,z3,costa,sinta,0)
    call Ldyns113(zi1,z2,z3,costa,sinta,1)
 
    call Rotcoordnt(xr1,xr2,xr3,xsc1,xsc2,xsa1,xsa2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(1) = zii1 ; zi1s(1) = zi1 ; z2s(1) = z2 ; z3s(1) = z3
    call Ldyns113(zii1,z2,z3,costa,sinta,2)
    call Ldyns113(zi1,z2,z3,costa,sinta,3)
 
    call Rotcoordnt(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(2) = zii1 ; zi1s(2) = zi1 ; z2s(2) = z2 ; z3s(2) = z3
    call Ldyns113(zii1,z2,z3,costa,sinta,4) 
    call Ldyns113(zi1,z2,z3,costa,sinta,5) 
  end subroutine Kdyns113

  subroutine Kdyns121(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3
    real(8) :: z1a,zi1,zii1,z2,z3,costa,sinta
 
    call Rotcoordnt(xr1,xr2,xr3,xsb1,xsb2,xsc1,xsc2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(0) = zii1 ; zi1s(0) = zi1 ; z2s(0) = z2 ; z3s(0) = z3
    call Ldyns121(zii1,z2,z3,costa,sinta,0)
    call Ldyns121(zi1,z2,z3,costa,sinta,1)
 
    call Rotcoordnt(xr1,xr2,xr3,xsc1,xsc2,xsa1,xsa2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(1) = zii1 ; zi1s(1) = zi1 ; z2s(1) = z2 ; z3s(1) = z3
    call Ldyns121(zii1,z2,z3,costa,sinta,2)
    call Ldyns121(zi1,z2,z3,costa,sinta,3)
 
    call Rotcoordnt(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(2) = zii1 ; zi1s(2) = zi1 ; z2s(2) = z2 ; z3s(2) = z3
    call Ldyns121(zii1,z2,z3,costa,sinta,4) 
    call Ldyns121(zi1,z2,z3,costa,sinta,5) 
  end subroutine Kdyns121

  subroutine Kdyns122(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3
    real(8) :: z1a,zi1,zii1,z2,z3,costa,sinta
 
    call Rotcoordnt(xr1,xr2,xr3,xsb1,xsb2,xsc1,xsc2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(0) = zii1 ; zi1s(0) = zi1 ; z2s(0) = z2 ; z3s(0) = z3
    call Ldyns122(zii1,z2,z3,costa,sinta,0)
    call Ldyns122(zi1,z2,z3,costa,sinta,1)
 
    call Rotcoordnt(xr1,xr2,xr3,xsc1,xsc2,xsa1,xsa2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(1) = zii1 ; zi1s(1) = zi1 ; z2s(1) = z2 ; z3s(1) = z3
    call Ldyns122(zii1,z2,z3,costa,sinta,2)
    call Ldyns122(zi1,z2,z3,costa,sinta,3)
 
    call Rotcoordnt(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(2) = zii1 ; zi1s(2) = zi1 ; z2s(2) = z2 ; z3s(2) = z3
    call Ldyns122(zii1,z2,z3,costa,sinta,4) 
    call Ldyns122(zi1,z2,z3,costa,sinta,5) 
  end subroutine Kdyns122

  subroutine Kdyns123(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3
    real(8) :: z1a,zi1,zii1,z2,z3,costa,sinta
 
    call Rotcoordnt(xr1,xr2,xr3,xsb1,xsb2,xsc1,xsc2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(0) = zii1 ; zi1s(0) = zi1 ; z2s(0) = z2 ; z3s(0) = z3
    call Ldyns123(zii1,z2,z3,costa,sinta,0)
    call Ldyns123(zi1,z2,z3,costa,sinta,1)
 
    call Rotcoordnt(xr1,xr2,xr3,xsc1,xsc2,xsa1,xsa2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(1) = zii1 ; zi1s(1) = zi1 ; z2s(1) = z2 ; z3s(1) = z3
    call Ldyns123(zii1,z2,z3,costa,sinta,2)
    call Ldyns123(zi1,z2,z3,costa,sinta,3)
 
    call Rotcoordnt(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(2) = zii1 ; zi1s(2) = zi1 ; z2s(2) = z2 ; z3s(2) = z3
    call Ldyns123(zii1,z2,z3,costa,sinta,4) 
    call Ldyns123(zi1,z2,z3,costa,sinta,5) 
  end subroutine Kdyns123

  subroutine Kdyns131(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3
    real(8) :: z1a,zi1,zii1,z2,z3,costa,sinta
 
    call Rotcoordnt(xr1,xr2,xr3,xsb1,xsb2,xsc1,xsc2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(0) = zii1 ; zi1s(0) = zi1 ; z2s(0) = z2 ; z3s(0) = z3
    call Ldyns131(zii1,z2,z3,costa,sinta,0)
    call Ldyns131(zi1,z2,z3,costa,sinta,1)
 
    call Rotcoordnt(xr1,xr2,xr3,xsc1,xsc2,xsa1,xsa2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(1) = zii1 ; zi1s(1) = zi1 ; z2s(1) = z2 ; z3s(1) = z3
    call Ldyns131(zii1,z2,z3,costa,sinta,2)
    call Ldyns131(zi1,z2,z3,costa,sinta,3)
 
    call Rotcoordnt(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(2) = zii1 ; zi1s(2) = zi1 ; z2s(2) = z2 ; z3s(2) = z3
    call Ldyns131(zii1,z2,z3,costa,sinta,4) 
    call Ldyns131(zi1,z2,z3,costa,sinta,5) 
  end subroutine Kdyns131

  subroutine Kdyns132(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3
    real(8) :: z1a,zi1,zii1,z2,z3,costa,sinta
 
    call Rotcoordnt(xr1,xr2,xr3,xsb1,xsb2,xsc1,xsc2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(0) = zii1 ; zi1s(0) = zi1 ; z2s(0) = z2 ; z3s(0) = z3
    call Ldyns132(zii1,z2,z3,costa,sinta,0)
    call Ldyns132(zi1,z2,z3,costa,sinta,1)
 
    call Rotcoordnt(xr1,xr2,xr3,xsc1,xsc2,xsa1,xsa2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(1) = zii1 ; zi1s(1) = zi1 ; z2s(1) = z2 ; z3s(1) = z3
    call Ldyns132(zii1,z2,z3,costa,sinta,2)
    call Ldyns132(zi1,z2,z3,costa,sinta,3)
 
    call Rotcoordnt(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(2) = zii1 ; zi1s(2) = zi1 ; z2s(2) = z2 ; z3s(2) = z3
    call Ldyns132(zii1,z2,z3,costa,sinta,4) 
    call Ldyns132(zi1,z2,z3,costa,sinta,5) 
  end subroutine Kdyns132

  subroutine Kdyns133(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3
    real(8) :: z1a,zi1,zii1,z2,z3,costa,sinta
 
    call Rotcoordnt(xr1,xr2,xr3,xsb1,xsb2,xsc1,xsc2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(0) = zii1 ; zi1s(0) = zi1 ; z2s(0) = z2 ; z3s(0) = z3
    call Ldyns133(zii1,z2,z3,costa,sinta,0)
    call Ldyns133(zi1,z2,z3,costa,sinta,1)
 
    call Rotcoordnt(xr1,xr2,xr3,xsc1,xsc2,xsa1,xsa2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(1) = zii1 ; zi1s(1) = zi1 ; z2s(1) = z2 ; z3s(1) = z3
    call Ldyns133(zii1,z2,z3,costa,sinta,2)
    call Ldyns133(zi1,z2,z3,costa,sinta,3)
 
    call Rotcoordnt(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(2) = zii1 ; zi1s(2) = zi1 ; z2s(2) = z2 ; z3s(2) = z3
    call Ldyns133(zii1,z2,z3,costa,sinta,4) 
    call Ldyns133(zi1,z2,z3,costa,sinta,5) 
  end subroutine Kdyns133

  subroutine Kdyns211(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3
    real(8) :: z1a,zi1,zii1,z2,z3,costa,sinta
 
    call Rotcoordnt(xr1,xr2,xr3,xsb1,xsb2,xsc1,xsc2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(0) = zii1 ; zi1s(0) = zi1 ; z2s(0) = z2 ; z3s(0) = z3
    call Ldyns211(zii1,z2,z3,costa,sinta,0)
    call Ldyns211(zi1,z2,z3,costa,sinta,1)
 
    call Rotcoordnt(xr1,xr2,xr3,xsc1,xsc2,xsa1,xsa2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(1) = zii1 ; zi1s(1) = zi1 ; z2s(1) = z2 ; z3s(1) = z3
    call Ldyns211(zii1,z2,z3,costa,sinta,2)
    call Ldyns211(zi1,z2,z3,costa,sinta,3)
 
    call Rotcoordnt(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(2) = zii1 ; zi1s(2) = zi1 ; z2s(2) = z2 ; z3s(2) = z3
    call Ldyns211(zii1,z2,z3,costa,sinta,4) 
    call Ldyns211(zi1,z2,z3,costa,sinta,5) 
  end subroutine Kdyns211

  subroutine Kdyns212(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3
    real(8) :: z1a,zi1,zii1,z2,z3,costa,sinta
 
    call Rotcoordnt(xr1,xr2,xr3,xsb1,xsb2,xsc1,xsc2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(0) = zii1 ; zi1s(0) = zi1 ; z2s(0) = z2 ; z3s(0) = z3
    call Ldyns212(zii1,z2,z3,costa,sinta,0)
    call Ldyns212(zi1,z2,z3,costa,sinta,1)
 
    call Rotcoordnt(xr1,xr2,xr3,xsc1,xsc2,xsa1,xsa2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(1) = zii1 ; zi1s(1) = zi1 ; z2s(1) = z2 ; z3s(1) = z3
    call Ldyns212(zii1,z2,z3,costa,sinta,2)
    call Ldyns212(zi1,z2,z3,costa,sinta,3)
 
    call Rotcoordnt(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(2) = zii1 ; zi1s(2) = zi1 ; z2s(2) = z2 ; z3s(2) = z3
    call Ldyns212(zii1,z2,z3,costa,sinta,4) 
    call Ldyns212(zi1,z2,z3,costa,sinta,5) 
  end subroutine Kdyns212

  subroutine Kdyns213(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3
    real(8) :: z1a,zi1,zii1,z2,z3,costa,sinta
 
    call Rotcoordnt(xr1,xr2,xr3,xsb1,xsb2,xsc1,xsc2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(0) = zii1 ; zi1s(0) = zi1 ; z2s(0) = z2 ; z3s(0) = z3
    call Ldyns213(zii1,z2,z3,costa,sinta,0)
    call Ldyns213(zi1,z2,z3,costa,sinta,1)
 
    call Rotcoordnt(xr1,xr2,xr3,xsc1,xsc2,xsa1,xsa2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(1) = zii1 ; zi1s(1) = zi1 ; z2s(1) = z2 ; z3s(1) = z3
    call Ldyns213(zii1,z2,z3,costa,sinta,2)
    call Ldyns213(zi1,z2,z3,costa,sinta,3)
 
    call Rotcoordnt(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(2) = zii1 ; zi1s(2) = zi1 ; z2s(2) = z2 ; z3s(2) = z3
    call Ldyns213(zii1,z2,z3,costa,sinta,4) 
    call Ldyns213(zi1,z2,z3,costa,sinta,5) 
  end subroutine Kdyns213

  subroutine Kdyns221(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3
    real(8) :: z1a,zi1,zii1,z2,z3,costa,sinta
 
    call Rotcoordnt(xr1,xr2,xr3,xsb1,xsb2,xsc1,xsc2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(0) = zii1 ; zi1s(0) = zi1 ; z2s(0) = z2 ; z3s(0) = z3
    call Ldyns221(zii1,z2,z3,costa,sinta,0)
    call Ldyns221(zi1,z2,z3,costa,sinta,1)
 
    call Rotcoordnt(xr1,xr2,xr3,xsc1,xsc2,xsa1,xsa2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(1) = zii1 ; zi1s(1) = zi1 ; z2s(1) = z2 ; z3s(1) = z3
    call Ldyns221(zii1,z2,z3,costa,sinta,2)
    call Ldyns221(zi1,z2,z3,costa,sinta,3)
 
    call Rotcoordnt(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(2) = zii1 ; zi1s(2) = zi1 ; z2s(2) = z2 ; z3s(2) = z3
    call Ldyns221(zii1,z2,z3,costa,sinta,4) 
    call Ldyns221(zi1,z2,z3,costa,sinta,5) 
  end subroutine Kdyns221

  subroutine Kdyns222(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3
    real(8) :: z1a,zi1,zii1,z2,z3,costa,sinta
 
    call Rotcoordnt(xr1,xr2,xr3,xsb1,xsb2,xsc1,xsc2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(0) = zii1 ; zi1s(0) = zi1 ; z2s(0) = z2 ; z3s(0) = z3
    call Ldyns222(zii1,z2,z3,costa,sinta,0)
    call Ldyns222(zi1,z2,z3,costa,sinta,1)
 
    call Rotcoordnt(xr1,xr2,xr3,xsc1,xsc2,xsa1,xsa2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(1) = zii1 ; zi1s(1) = zi1 ; z2s(1) = z2 ; z3s(1) = z3
    call Ldyns222(zii1,z2,z3,costa,sinta,2)
    call Ldyns222(zi1,z2,z3,costa,sinta,3)
 
    call Rotcoordnt(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(2) = zii1 ; zi1s(2) = zi1 ; z2s(2) = z2 ; z3s(2) = z3
    call Ldyns222(zii1,z2,z3,costa,sinta,4) 
    call Ldyns222(zi1,z2,z3,costa,sinta,5) 
  end subroutine Kdyns222

  subroutine Kdyns223(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3
    real(8) :: z1a,zi1,zii1,z2,z3,costa,sinta
 
    call Rotcoordnt(xr1,xr2,xr3,xsb1,xsb2,xsc1,xsc2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(0) = zii1 ; zi1s(0) = zi1 ; z2s(0) = z2 ; z3s(0) = z3
    call Ldyns223(zii1,z2,z3,costa,sinta,0)
    call Ldyns223(zi1,z2,z3,costa,sinta,1)
 
    call Rotcoordnt(xr1,xr2,xr3,xsc1,xsc2,xsa1,xsa2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(1) = zii1 ; zi1s(1) = zi1 ; z2s(1) = z2 ; z3s(1) = z3
    call Ldyns223(zii1,z2,z3,costa,sinta,2)
    call Ldyns223(zi1,z2,z3,costa,sinta,3)
 
    call Rotcoordnt(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(2) = zii1 ; zi1s(2) = zi1 ; z2s(2) = z2 ; z3s(2) = z3
    call Ldyns223(zii1,z2,z3,costa,sinta,4) 
    call Ldyns223(zi1,z2,z3,costa,sinta,5) 
  end subroutine Kdyns223

  subroutine Kdyns231(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3
    real(8) :: z1a,zi1,zii1,z2,z3,costa,sinta
 
    call Rotcoordnt(xr1,xr2,xr3,xsb1,xsb2,xsc1,xsc2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(0) = zii1 ; zi1s(0) = zi1 ; z2s(0) = z2 ; z3s(0) = z3
    call Ldyns231(zii1,z2,z3,costa,sinta,0)
    call Ldyns231(zi1,z2,z3,costa,sinta,1)
 
    call Rotcoordnt(xr1,xr2,xr3,xsc1,xsc2,xsa1,xsa2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(1) = zii1 ; zi1s(1) = zi1 ; z2s(1) = z2 ; z3s(1) = z3
    call Ldyns231(zii1,z2,z3,costa,sinta,2)
    call Ldyns231(zi1,z2,z3,costa,sinta,3)
 
    call Rotcoordnt(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(2) = zii1 ; zi1s(2) = zi1 ; z2s(2) = z2 ; z3s(2) = z3
    call Ldyns231(zii1,z2,z3,costa,sinta,4) 
    call Ldyns231(zi1,z2,z3,costa,sinta,5) 
  end subroutine Kdyns231

  subroutine Kdyns232(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3
    real(8) :: z1a,zi1,zii1,z2,z3,costa,sinta
 
    call Rotcoordnt(xr1,xr2,xr3,xsb1,xsb2,xsc1,xsc2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(0) = zii1 ; zi1s(0) = zi1 ; z2s(0) = z2 ; z3s(0) = z3
    call Ldyns232(zii1,z2,z3,costa,sinta,0)
    call Ldyns232(zi1,z2,z3,costa,sinta,1)
 
    call Rotcoordnt(xr1,xr2,xr3,xsc1,xsc2,xsa1,xsa2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(1) = zii1 ; zi1s(1) = zi1 ; z2s(1) = z2 ; z3s(1) = z3
    call Ldyns232(zii1,z2,z3,costa,sinta,2)
    call Ldyns232(zi1,z2,z3,costa,sinta,3)
 
    call Rotcoordnt(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(2) = zii1 ; zi1s(2) = zi1 ; z2s(2) = z2 ; z3s(2) = z3
    call Ldyns232(zii1,z2,z3,costa,sinta,4) 
    call Ldyns232(zi1,z2,z3,costa,sinta,5) 
  end subroutine Kdyns232

  subroutine Kdyns233(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3
    real(8) :: z1a,zi1,zii1,z2,z3,costa,sinta
 
    call Rotcoordnt(xr1,xr2,xr3,xsb1,xsb2,xsc1,xsc2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(0) = zii1 ; zi1s(0) = zi1 ; z2s(0) = z2 ; z3s(0) = z3
    call Ldyns233(zii1,z2,z3,costa,sinta,0)
    call Ldyns233(zi1,z2,z3,costa,sinta,1)
 
    call Rotcoordnt(xr1,xr2,xr3,xsc1,xsc2,xsa1,xsa2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(1) = zii1 ; zi1s(1) = zi1 ; z2s(1) = z2 ; z3s(1) = z3
    call Ldyns233(zii1,z2,z3,costa,sinta,2)
    call Ldyns233(zi1,z2,z3,costa,sinta,3)
 
    call Rotcoordnt(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(2) = zii1 ; zi1s(2) = zi1 ; z2s(2) = z2 ; z3s(2) = z3
    call Ldyns233(zii1,z2,z3,costa,sinta,4) 
    call Ldyns233(zi1,z2,z3,costa,sinta,5) 
  end subroutine Kdyns233

  subroutine Kdyns311(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3
    real(8) :: z1a,zi1,zii1,z2,z3,costa,sinta
 
    call Rotcoordnt(xr1,xr2,xr3,xsb1,xsb2,xsc1,xsc2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(0) = zii1 ; zi1s(0) = zi1 ; z2s(0) = z2 ; z3s(0) = z3
    call Ldyns311(zii1,z2,z3,costa,sinta,0)
    call Ldyns311(zi1,z2,z3,costa,sinta,1)
 
    call Rotcoordnt(xr1,xr2,xr3,xsc1,xsc2,xsa1,xsa2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(1) = zii1 ; zi1s(1) = zi1 ; z2s(1) = z2 ; z3s(1) = z3
    call Ldyns311(zii1,z2,z3,costa,sinta,2)
    call Ldyns311(zi1,z2,z3,costa,sinta,3)
 
    call Rotcoordnt(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(2) = zii1 ; zi1s(2) = zi1 ; z2s(2) = z2 ; z3s(2) = z3
    call Ldyns311(zii1,z2,z3,costa,sinta,4) 
    call Ldyns311(zi1,z2,z3,costa,sinta,5) 
  end subroutine Kdyns311

  subroutine Kdyns312(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3
    real(8) :: z1a,zi1,zii1,z2,z3,costa,sinta
 
    call Rotcoordnt(xr1,xr2,xr3,xsb1,xsb2,xsc1,xsc2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(0) = zii1 ; zi1s(0) = zi1 ; z2s(0) = z2 ; z3s(0) = z3
    call Ldyns312(zii1,z2,z3,costa,sinta,0)
    call Ldyns312(zi1,z2,z3,costa,sinta,1)
 
    call Rotcoordnt(xr1,xr2,xr3,xsc1,xsc2,xsa1,xsa2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(1) = zii1 ; zi1s(1) = zi1 ; z2s(1) = z2 ; z3s(1) = z3
    call Ldyns312(zii1,z2,z3,costa,sinta,2)
    call Ldyns312(zi1,z2,z3,costa,sinta,3)
 
    call Rotcoordnt(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(2) = zii1 ; zi1s(2) = zi1 ; z2s(2) = z2 ; z3s(2) = z3
    call Ldyns312(zii1,z2,z3,costa,sinta,4) 
    call Ldyns312(zi1,z2,z3,costa,sinta,5) 
  end subroutine Kdyns312

  subroutine Kdyns313(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3
    real(8) :: z1a,zi1,zii1,z2,z3,costa,sinta
 
    call Rotcoordnt(xr1,xr2,xr3,xsb1,xsb2,xsc1,xsc2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(0) = zii1 ; zi1s(0) = zi1 ; z2s(0) = z2 ; z3s(0) = z3
    call Ldyns313(zii1,z2,z3,costa,sinta,0)
    call Ldyns313(zi1,z2,z3,costa,sinta,1)
 
    call Rotcoordnt(xr1,xr2,xr3,xsc1,xsc2,xsa1,xsa2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(1) = zii1 ; zi1s(1) = zi1 ; z2s(1) = z2 ; z3s(1) = z3
    call Ldyns313(zii1,z2,z3,costa,sinta,2)
    call Ldyns313(zi1,z2,z3,costa,sinta,3)
 
    call Rotcoordnt(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(2) = zii1 ; zi1s(2) = zi1 ; z2s(2) = z2 ; z3s(2) = z3
    call Ldyns313(zii1,z2,z3,costa,sinta,4) 
    call Ldyns313(zi1,z2,z3,costa,sinta,5) 
  end subroutine Kdyns313

  subroutine Kdyns321(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3
    real(8) :: z1a,zi1,zii1,z2,z3,costa,sinta
 
    call Rotcoordnt(xr1,xr2,xr3,xsb1,xsb2,xsc1,xsc2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(0) = zii1 ; zi1s(0) = zi1 ; z2s(0) = z2 ; z3s(0) = z3
    call Ldyns321(zii1,z2,z3,costa,sinta,0)
    call Ldyns321(zi1,z2,z3,costa,sinta,1)
 
    call Rotcoordnt(xr1,xr2,xr3,xsc1,xsc2,xsa1,xsa2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(1) = zii1 ; zi1s(1) = zi1 ; z2s(1) = z2 ; z3s(1) = z3
    call Ldyns321(zii1,z2,z3,costa,sinta,2)
    call Ldyns321(zi1,z2,z3,costa,sinta,3)
 
    call Rotcoordnt(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(2) = zii1 ; zi1s(2) = zi1 ; z2s(2) = z2 ; z3s(2) = z3
    call Ldyns321(zii1,z2,z3,costa,sinta,4) 
    call Ldyns321(zi1,z2,z3,costa,sinta,5) 
  end subroutine Kdyns321

  subroutine Kdyns322(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3
    real(8) :: z1a,zi1,zii1,z2,z3,costa,sinta
 
    call Rotcoordnt(xr1,xr2,xr3,xsb1,xsb2,xsc1,xsc2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(0) = zii1 ; zi1s(0) = zi1 ; z2s(0) = z2 ; z3s(0) = z3
    call Ldyns322(zii1,z2,z3,costa,sinta,0)
    call Ldyns322(zi1,z2,z3,costa,sinta,1)
 
    call Rotcoordnt(xr1,xr2,xr3,xsc1,xsc2,xsa1,xsa2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(1) = zii1 ; zi1s(1) = zi1 ; z2s(1) = z2 ; z3s(1) = z3
    call Ldyns322(zii1,z2,z3,costa,sinta,2)
    call Ldyns322(zi1,z2,z3,costa,sinta,3)
 
    call Rotcoordnt(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(2) = zii1 ; zi1s(2) = zi1 ; z2s(2) = z2 ; z3s(2) = z3
    call Ldyns322(zii1,z2,z3,costa,sinta,4) 
    call Ldyns322(zi1,z2,z3,costa,sinta,5) 
  end subroutine Kdyns322

  subroutine Kdyns323(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3
    real(8) :: z1a,zi1,zii1,z2,z3,costa,sinta
 
    call Rotcoordnt(xr1,xr2,xr3,xsb1,xsb2,xsc1,xsc2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(0) = zii1 ; zi1s(0) = zi1 ; z2s(0) = z2 ; z3s(0) = z3
    call Ldyns323(zii1,z2,z3,costa,sinta,0)
    call Ldyns323(zi1,z2,z3,costa,sinta,1)
 
    call Rotcoordnt(xr1,xr2,xr3,xsc1,xsc2,xsa1,xsa2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(1) = zii1 ; zi1s(1) = zi1 ; z2s(1) = z2 ; z3s(1) = z3
    call Ldyns323(zii1,z2,z3,costa,sinta,2)
    call Ldyns323(zi1,z2,z3,costa,sinta,3)
 
    call Rotcoordnt(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(2) = zii1 ; zi1s(2) = zi1 ; z2s(2) = z2 ; z3s(2) = z3
    call Ldyns323(zii1,z2,z3,costa,sinta,4) 
    call Ldyns323(zi1,z2,z3,costa,sinta,5) 
  end subroutine Kdyns323

  subroutine Kdyns331(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3
    real(8) :: z1a,zi1,zii1,z2,z3,costa,sinta
 
    call Rotcoordnt(xr1,xr2,xr3,xsb1,xsb2,xsc1,xsc2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(0) = zii1 ; zi1s(0) = zi1 ; z2s(0) = z2 ; z3s(0) = z3
    call Ldyns331(zii1,z2,z3,costa,sinta,0)
    call Ldyns331(zi1,z2,z3,costa,sinta,1)
 
    call Rotcoordnt(xr1,xr2,xr3,xsc1,xsc2,xsa1,xsa2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(1) = zii1 ; zi1s(1) = zi1 ; z2s(1) = z2 ; z3s(1) = z3
    call Ldyns331(zii1,z2,z3,costa,sinta,2)
    call Ldyns331(zi1,z2,z3,costa,sinta,3)
 
    call Rotcoordnt(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(2) = zii1 ; zi1s(2) = zi1 ; z2s(2) = z2 ; z3s(2) = z3
    call Ldyns331(zii1,z2,z3,costa,sinta,4) 
    call Ldyns331(zi1,z2,z3,costa,sinta,5) 
  end subroutine Kdyns331

  subroutine Kdyns332(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3
    real(8) :: z1a,zi1,zii1,z2,z3,costa,sinta
 
    call Rotcoordnt(xr1,xr2,xr3,xsb1,xsb2,xsc1,xsc2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(0) = zii1 ; zi1s(0) = zi1 ; z2s(0) = z2 ; z3s(0) = z3
    call Ldyns332(zii1,z2,z3,costa,sinta,0)
    call Ldyns332(zi1,z2,z3,costa,sinta,1)
 
    call Rotcoordnt(xr1,xr2,xr3,xsc1,xsc2,xsa1,xsa2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(1) = zii1 ; zi1s(1) = zi1 ; z2s(1) = z2 ; z3s(1) = z3
    call Ldyns332(zii1,z2,z3,costa,sinta,2)
    call Ldyns332(zi1,z2,z3,costa,sinta,3)
 
    call Rotcoordnt(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(2) = zii1 ; zi1s(2) = zi1 ; z2s(2) = z2 ; z3s(2) = z3
    call Ldyns332(zii1,z2,z3,costa,sinta,4) 
    call Ldyns332(zi1,z2,z3,costa,sinta,5) 
  end subroutine Kdyns332

  subroutine Kdyns333(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3
    real(8) :: z1a,zi1,zii1,z2,z3,costa,sinta
 
    call Rotcoordnt(xr1,xr2,xr3,xsb1,xsb2,xsc1,xsc2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(0) = zii1 ; zi1s(0) = zi1 ; z2s(0) = z2 ; z3s(0) = z3
    call Ldyns333(zii1,z2,z3,costa,sinta,0)
    call Ldyns333(zi1,z2,z3,costa,sinta,1)
 
    call Rotcoordnt(xr1,xr2,xr3,xsc1,xsc2,xsa1,xsa2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(1) = zii1 ; zi1s(1) = zi1 ; z2s(1) = z2 ; z3s(1) = z3
    call Ldyns333(zii1,z2,z3,costa,sinta,2)
    call Ldyns333(zi1,z2,z3,costa,sinta,3)
 
    call Rotcoordnt(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xs3,costa,sinta,zi1,zii1,z2,z3)
    zii1s(2) = zii1 ; zi1s(2) = zi1 ; z2s(2) = z2 ; z3s(2) = z3
    call Ldyns333(zii1,z2,z3,costa,sinta,4) 
    call Ldyns333(zi1,z2,z3,costa,sinta,5) 
  end subroutine Kdyns333

  !(subroutine Kdynf: compute stress Green function L created by constant slip rate of unit magnitude)
  real(8) function Kdynf111(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt
    real(8) :: L111,L111t

    L111t = Ldynf111(zii1s(0),z2s(0),z3s(0),t,0)&
         - Ldynf111(zi1s(0),z2s(0),z3s(0),t,1)&
         + Ldynf111(zii1s(1),z2s(1),z3s(1),t,2)&
         - Ldynf111(zi1s(1),z2s(1),z3s(1),t,3)&
         + Ldynf111(zii1s(2),z2s(2),z3s(2),t,4)&
         - Ldynf111(zi1s(2),z2s(2),z3s(2),t,5)

    L111 = Ldynf111(zii1s(0),z2s(0),z3s(0),t-delt,0)&
         - Ldynf111(zi1s(0),z2s(0),z3s(0),t-delt,1)&
         + Ldynf111(zii1s(1),z2s(1),z3s(1),t-delt,2)&
         - Ldynf111(zi1s(1),z2s(1),z3s(1),t-delt,3)&
         + Ldynf111(zii1s(2),z2s(2),z3s(2),t-delt,4)&
         - Ldynf111(zi1s(2),z2s(2),z3s(2),t-delt,5)

    L111 = L111t - L111
    Kdynf111 = L111
  end function Kdynf111

  real(8) function Kdynf112(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt
    real(8) :: L112,L112t

    L112t = Ldynf112(zii1s(0),z2s(0),z3s(0),t,0)&
          - Ldynf112(zi1s(0),z2s(0),z3s(0),t,1)&
          + Ldynf112(zii1s(1),z2s(1),z3s(1),t,2)&
          - Ldynf112(zi1s(1),z2s(1),z3s(1),t,3)&
          + Ldynf112(zii1s(2),z2s(2),z3s(2),t,4)&
          - Ldynf112(zi1s(2),z2s(2),z3s(2),t,5)

    L112 = Ldynf112(zii1s(0),z2s(0),z3s(0),t-delt,0)&
          - Ldynf112(zi1s(0),z2s(0),z3s(0),t-delt,1)&
          + Ldynf112(zii1s(1),z2s(1),z3s(1),t-delt,2)&
          - Ldynf112(zi1s(1),z2s(1),z3s(1),t-delt,3)&
          + Ldynf112(zii1s(2),z2s(2),z3s(2),t-delt,4)&
          - Ldynf112(zi1s(2),z2s(2),z3s(2),t-delt,5)
      
    L112 = L112t - L112
    Kdynf112 = L112
  end function Kdynf112

  real(8) function Kdynf113(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt
    real(8) :: L113,L113t

    L113t = Ldynf113(zii1s(0),z2s(0),z3s(0),t,0)&
          - Ldynf113(zi1s(0),z2s(0),z3s(0),t,1)&
          + Ldynf113(zii1s(1),z2s(1),z3s(1),t,2)&
          - Ldynf113(zi1s(1),z2s(1),z3s(1),t,3)&
          + Ldynf113(zii1s(2),z2s(2),z3s(2),t,4)&
          - Ldynf113(zi1s(2),z2s(2),z3s(2),t,5)&
          + Caa3*Hevix( z2s(0) )*Hevix( z2s(1) )*Hevix( z2s(2) )*Heviw( t-abs(z3s(0))/Cl )

    L113 = Ldynf113(zii1s(0),z2s(0),z3s(0),t-delt,0)&
         - Ldynf113(zi1s(0),z2s(0),z3s(0),t-delt,1)&
         + Ldynf113(zii1s(1),z2s(1),z3s(1),t-delt,2)&
         - Ldynf113(zi1s(1),z2s(1),z3s(1),t-delt,3)&
         + Ldynf113(zii1s(2),z2s(2),z3s(2),t-delt,4)&
         - Ldynf113(zi1s(2),z2s(2),z3s(2),t-delt,5)&
         + Caa3*Hevix( z2s(0) )*Hevix( z2s(1) )*Hevix( z2s(2) )*Heviw( t-delt-abs(z3s(0))/Cl )

    L113 = L113t - L113
    Kdynf113 = L113
  end function Kdynf113

  real(8) function Kdynf121(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt
    real(8)::L121,L121t

    L121t = Ldynf121(zii1s(0),z2s(0),z3s(0),t,0)&
          - Ldynf121(zi1s(0),z2s(0),z3s(0),t,1)&
          + Ldynf121(zii1s(1),z2s(1),z3s(1),t,2)&
          - Ldynf121(zi1s(1),z2s(1),z3s(1),t,3)&
          + Ldynf121(zii1s(2),z2s(2),z3s(2),t,4)&
          - Ldynf121(zi1s(2),z2s(2),z3s(2),t,5)
    
    L121 = Ldynf121(zii1s(0),z2s(0),z3s(0),t-delt,0)&
         - Ldynf121(zi1s(0),z2s(0),z3s(0),t-delt,1)&
         + Ldynf121(zii1s(1),z2s(1),z3s(1),t-delt,2)&
         - Ldynf121(zi1s(1),z2s(1),z3s(1),t-delt,3)&
         + Ldynf121(zii1s(2),z2s(2),z3s(2),t-delt,4)&
         - Ldynf121(zi1s(2),z2s(2),z3s(2),t-delt,5)

    L121 = L121t - L121
    Kdynf121 = L121
  end function Kdynf121

  real(8) function Kdynf122(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt
    real(8) :: L122,L122t

    L122 = Ldynf122(zii1s(0),z2s(0),z3s(0),t,0)&
         - Ldynf122(zi1s(0),z2s(0),z3s(0),t,1)&
         + Ldynf122(zii1s(1),z2s(1),z3s(1),t,2)&
         - Ldynf122(zi1s(1),z2s(1),z3s(1),t,3)&
         + Ldynf122(zii1s(2),z2s(2),z3s(2),t,4)&
         - Ldynf122(zi1s(2),z2s(2),z3s(2),t,5)

    L122 = Ldynf122(zii1s(0),z2s(0),z3s(0),t-delt,0)&
         - Ldynf122(zi1s(0),z2s(0),z3s(0),t-delt,1)&
         + Ldynf122(zii1s(1),z2s(1),z3s(1),t-delt,2)&
         - Ldynf122(zi1s(1),z2s(1),z3s(1),t-delt,3)&
         + Ldynf122(zii1s(2),z2s(2),z3s(2),t-delt,4)&
         - Ldynf122(zi1s(2),z2s(2),z3s(2),t-delt,5)

    L122 = L122t - L122
    Kdynf122 = L122
  end function Kdynf122


  real(8) function Kdynf123(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt
    real(8) :: L123,L123t

    L123t = Ldynf123(zii1s(0),z2s(0),z3s(0),t,0)&
          - Ldynf123(zi1s(0),z2s(0),z3s(0),t,1)&
          + Ldynf123(zii1s(1),z2s(1),z3s(1),t,2)&
          - Ldynf123(zi1s(1),z2s(1),z3s(1),t,3)&
          + Ldynf123(zii1s(2),z2s(2),z3s(2),t,4)&
          - Ldynf123(zi1s(2),z2s(2),z3s(2),t,5)

    L123 = Ldynf123(zii1s(0),z2s(0),z3s(0),t-delt,0)&
         - Ldynf123(zi1s(0),z2s(0),z3s(0),t-delt,1)&
         + Ldynf123(zii1s(1),z2s(1),z3s(1),t-delt,2)&
         - Ldynf123(zi1s(1),z2s(1),z3s(1),t-delt,3)&
         + Ldynf123(zii1s(2),z2s(2),z3s(2),t-delt,4)&
         - Ldynf123(zi1s(2),z2s(2),z3s(2),t-delt,5)

    L123 = L123t - L123
    Kdynf123 = L123
  end function Kdynf123

  real(8) function Kdynf131(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt
    real(8) :: L131,L131t
    
    L131t = Ldynf131(zii1s(0),z2s(0),z3s(0),t,0)&
          - Ldynf131(zi1s(0),z2s(0),z3s(0),t,1)&
          + Ldynf131(zii1s(1),z2s(1),z3s(1),t,2)&
          - Ldynf131(zi1s(1),z2s(1),z3s(1),t,3)&
          + Ldynf131(zii1s(2),z2s(2),z3s(2),t,4)&
          - Ldynf131(zi1s(2),z2s(2),z3s(2),t,5)

    L131 = Ldynf131(zii1s(0),z2s(0),z3s(0),t-delt,0)&
         - Ldynf131(zi1s(0),z2s(0),z3s(0),t-delt,1)&
         + Ldynf131(zii1s(1),z2s(1),z3s(1),t-delt,2)&
         - Ldynf131(zi1s(1),z2s(1),z3s(1),t-delt,3)&
         + Ldynf131(zii1s(2),z2s(2),z3s(2),t-delt,4)&
         - Ldynf131(zi1s(2),z2s(2),z3s(2),t-delt,5)

    L131 = L131t - L131
    Kdynf131 = L131
  end function Kdynf131

  real(8) function Kdynf132(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt
    real(8) :: L132,L132t
    
    L132t = Ldynf132(zii1s(0),z2s(0),z3s(0),t,0)&
          - Ldynf132(zi1s(0),z2s(0),z3s(0),t,1)&
          + Ldynf132(zii1s(1),z2s(1),z3s(1),t,2)&
          - Ldynf132(zi1s(1),z2s(1),z3s(1),t,3)&
          + Ldynf132(zii1s(2),z2s(2),z3s(2),t,4)&
          - Ldynf132(zi1s(2),z2s(2),z3s(2),t,5)

    L132 = Ldynf132(zii1s(0),z2s(0),z3s(0),t-delt,0)&
         - Ldynf132(zi1s(0),z2s(0),z3s(0),t-delt,1)&
         + Ldynf132(zii1s(1),z2s(1),z3s(1),t-delt,2)&
         - Ldynf132(zi1s(1),z2s(1),z3s(1),t-delt,3)&
         + Ldynf132(zii1s(2),z2s(2),z3s(2),t-delt,4)&
         - Ldynf132(zi1s(2),z2s(2),z3s(2),t-delt,5)

    L132 = L132t - L132
    Kdynf132 = L132
  end function Kdynf132

  real(8) function Kdynf133(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt
    real(8) :: L133,L133t
    
    L133t = Ldynf133(zii1s(0),z2s(0),z3s(0),t,0)&
          - Ldynf133(zi1s(0),z2s(0),z3s(0),t,1)&
          + Ldynf133(zii1s(1),z2s(1),z3s(1),t,2)&
          - Ldynf133(zi1s(1),z2s(1),z3s(1),t,3)&
          + Ldynf133(zii1s(2),z2s(2),z3s(2),t,4)&
          - Ldynf133(zi1s(2),z2s(2),z3s(2),t,5)

    L133 = Ldynf133(zii1s(0),z2s(0),z3s(0),t-delt,0)&
         - Ldynf133(zi1s(0),z2s(0),z3s(0),t-delt,1)&
         + Ldynf133(zii1s(1),z2s(1),z3s(1),t-delt,2)&
         - Ldynf133(zi1s(1),z2s(1),z3s(1),t-delt,3)&
         + Ldynf133(zii1s(2),z2s(2),z3s(2),t-delt,4)&
         - Ldynf133(zi1s(2),z2s(2),z3s(2),t-delt,5)

    L133 = L133t - L133
    Kdynf133 = L133
  end function Kdynf133

  real(8) function Kdynf211(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt
    real(8) :: L211,L211t

    L211t = Ldynf211(zii1s(0),z2s(0),z3s(0),t,0)&
          - Ldynf211(zi1s(0),z2s(0),z3s(0),t,1)&
          + Ldynf211(zii1s(1),z2s(1),z3s(1),t,2)&
          - Ldynf211(zi1s(1),z2s(1),z3s(1),t,3)&
          + Ldynf211(zii1s(2),z2s(2),z3s(2),t,4)&
          - Ldynf211(zi1s(2),z2s(2),z3s(2),t,5)

    L211 = Ldynf211(zii1s(0),z2s(0),z3s(0),t-delt,0)&
         - Ldynf211(zi1s(0),z2s(0),z3s(0),t-delt,1)&
         + Ldynf211(zii1s(1),z2s(1),z3s(1),t-delt,2)&
         - Ldynf211(zi1s(1),z2s(1),z3s(1),t-delt,3)&
         + Ldynf211(zii1s(2),z2s(2),z3s(2),t-delt,4)&
         - Ldynf211(zi1s(2),z2s(2),z3s(2),t-delt,5)

    L211 = L211t - L211
    Kdynf211 = L211
  end function Kdynf211

  real(8) function Kdynf212(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt
    real(8) :: L212,L212t

    L212t = Ldynf212(zii1s(0),z2s(0),z3s(0),t,0)&
          - Ldynf212(zi1s(0),z2s(0),z3s(0),t,1)&
          + Ldynf212(zii1s(1),z2s(1),z3s(1),t,2)&
          - Ldynf212(zi1s(1),z2s(1),z3s(1),t,3)&
          + Ldynf212(zii1s(2),z2s(2),z3s(2),t,4)&
          - Ldynf212(zi1s(2),z2s(2),z3s(2),t,5)

    L212 = Ldynf212(zii1s(0),z2s(0),z3s(0),t-delt,0)&
         - Ldynf212(zi1s(0),z2s(0),z3s(0),t-delt,1)&
         + Ldynf212(zii1s(1),z2s(1),z3s(1),t-delt,2)&
         - Ldynf212(zi1s(1),z2s(1),z3s(1),t-delt,3)&
         + Ldynf212(zii1s(2),z2s(2),z3s(2),t-delt,4)&
         - Ldynf212(zi1s(2),z2s(2),z3s(2),t-delt,5)

    L212 = L212t - L212
    Kdynf212 = L212
  end function Kdynf212

  real(8) function Kdynf213(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt
    real(8) :: L213,L213t

    L213t = Ldynf213(zii1s(0),z2s(0),z3s(0),t,0)&
          - Ldynf213(zi1s(0),z2s(0),z3s(0),t,1)&
          + Ldynf213(zii1s(1),z2s(1),z3s(1),t,2)&
          - Ldynf213(zi1s(1),z2s(1),z3s(1),t,3)&
          + Ldynf213(zii1s(2),z2s(2),z3s(2),t,4)&
          - Ldynf213(zi1s(2),z2s(2),z3s(2),t,5)

    L213 = Ldynf213(zii1s(0),z2s(0),z3s(0),t-delt,0)&
         - Ldynf213(zi1s(0),z2s(0),z3s(0),t-delt,1)&
         + Ldynf213(zii1s(1),z2s(1),z3s(1),t-delt,2)&
         - Ldynf213(zi1s(1),z2s(1),z3s(1),t-delt,3)&
         + Ldynf213(zii1s(2),z2s(2),z3s(2),t-delt,4)&
         - Ldynf213(zi1s(2),z2s(2),z3s(2),t-delt,5)

    L213 = L213t - L213
    Kdynf213 = L213
  end function Kdynf213

  real(8) function Kdynf221(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt
    real(8) :: L221,L221t

    L221t = Ldynf221(zii1s(0),z2s(0),z3s(0),t,0)&
          - Ldynf221(zi1s(0),z2s(0),z3s(0),t,1)&
          + Ldynf221(zii1s(1),z2s(1),z3s(1),t,2)&
          - Ldynf221(zi1s(1),z2s(1),z3s(1),t,3)&
          + Ldynf221(zii1s(2),z2s(2),z3s(2),t,4)&
          - Ldynf221(zi1s(2),z2s(2),z3s(2),t,5)

    L221 = Ldynf221(zii1s(0),z2s(0),z3s(0),t-delt,0)&
         - Ldynf221(zi1s(0),z2s(0),z3s(0),t-delt,1)&
         + Ldynf221(zii1s(1),z2s(1),z3s(1),t-delt,2)&
         - Ldynf221(zi1s(1),z2s(1),z3s(1),t-delt,3)&
         + Ldynf221(zii1s(2),z2s(2),z3s(2),t-delt,4)&
         - Ldynf221(zi1s(2),z2s(2),z3s(2),t-delt,5)

    L221 = L221t - L221
    Kdynf221 = L221
  end function Kdynf221

  real(8) function Kdynf222(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt
    real(8) :: L222,L222t

    L222t = Ldynf222(zii1s(0),z2s(0),z3s(0),t,0)&
          - Ldynf222(zi1s(0),z2s(0),z3s(0),t,1)&
          + Ldynf222(zii1s(1),z2s(1),z3s(1),t,2)&
          - Ldynf222(zi1s(1),z2s(1),z3s(1),t,3)&
          + Ldynf222(zii1s(2),z2s(2),z3s(2),t,4)&
          - Ldynf222(zi1s(2),z2s(2),z3s(2),t,5)
      
    L222 = Ldynf222(zii1s(0),z2s(0),z3s(0),t-delt,0)&
         - Ldynf222(zi1s(0),z2s(0),z3s(0),t-delt,1)&
         + Ldynf222(zii1s(1),z2s(1),z3s(1),t-delt,2)&
         - Ldynf222(zi1s(1),z2s(1),z3s(1),t-delt,3)&
         + Ldynf222(zii1s(2),z2s(2),z3s(2),t-delt,4)&
         - Ldynf222(zi1s(2),z2s(2),z3s(2),t-delt,5)

    L222 = L222t - L222
    Kdynf222 = L222
  end function Kdynf222

  real(8) function Kdynf223(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt
    real(8) :: L223,L223t

    L223 = Ldynf223(zii1s(0),z2s(0),z3s(0),t,0)&
         - Ldynf223(zi1s(0),z2s(0),z3s(0),t,1)&
         + Ldynf223(zii1s(1),z2s(1),z3s(1),t,2)&
         - Ldynf223(zi1s(1),z2s(1),z3s(1),t,3)&
         + Ldynf223(zii1s(2),z2s(2),z3s(2),t,4)&
         - Ldynf223(zi1s(2),z2s(2),z3s(2),t,5)&
         + Caa3*Hevix( z2s(0) )*Hevix( z2s(1) )*Hevix( z2s(2) )*Heviw( t-abs(z3s(0))/Cl )

    L223 = Ldynf223(zii1s(0),z2s(0),z3s(0),t-delt,0)&
         - Ldynf223(zi1s(0),z2s(0),z3s(0),t-delt,1)&
         + Ldynf223(zii1s(1),z2s(1),z3s(1),t-delt,2)&
         - Ldynf223(zi1s(1),z2s(1),z3s(1),t-delt,3)&
         + Ldynf223(zii1s(2),z2s(2),z3s(2),t-delt,4)&
         - Ldynf223(zi1s(2),z2s(2),z3s(2),t-delt,5)&
         + Caa3*Hevix( z2s(0) )*Hevix( z2s(1) )*Hevix( z2s(2) )*Heviw( t-delt-abs(z3s(0))/Cl )

    L223 = L223t - L223
    Kdynf223 = L223
  end function Kdynf223

  real(8) function Kdynf231(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt
    real(8) :: L231,L231t

    L231t = Ldynf231(zii1s(0),z2s(0),z3s(0),t,0)&
          - Ldynf231(zi1s(0),z2s(0),z3s(0),t,1)&
          + Ldynf231(zii1s(1),z2s(1),z3s(1),t,2)&
          - Ldynf231(zi1s(1),z2s(1),z3s(1),t,3)&
          + Ldynf231(zii1s(2),z2s(2),z3s(2),t,4)&
          - Ldynf231(zi1s(2),z2s(2),z3s(2),t,5)

    L231 = Ldynf231(zii1s(0),z2s(0),z3s(0),t-delt,0)&
         - Ldynf231(zi1s(0),z2s(0),z3s(0),t-delt,1)&
         + Ldynf231(zii1s(1),z2s(1),z3s(1),t-delt,2)&
         - Ldynf231(zi1s(1),z2s(1),z3s(1),t-delt,3)&
         + Ldynf231(zii1s(2),z2s(2),z3s(2),t-delt,4)&
         - Ldynf231(zi1s(2),z2s(2),z3s(2),t-delt,5)

    L231 = L231t - L231
    Kdynf231 = L231
  end function Kdynf231

  real(8) function Kdynf232(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt
    real(8) :: L232,L232t

    L232t = Ldynf232(zii1s(0),z2s(0),z3s(0),t,0)&
          - Ldynf232(zi1s(0),z2s(0),z3s(0),t,1)&
          + Ldynf232(zii1s(1),z2s(1),z3s(1),t,2)&
          - Ldynf232(zi1s(1),z2s(1),z3s(1),t,3)&
          + Ldynf232(zii1s(2),z2s(2),z3s(2),t,4)&
          - Ldynf232(zi1s(2),z2s(2),z3s(2),t,5)

    L232 = Ldynf232(zii1s(0),z2s(0),z3s(0),t-delt,0)&
         - Ldynf232(zi1s(0),z2s(0),z3s(0),t-delt,1)&
         + Ldynf232(zii1s(1),z2s(1),z3s(1),t-delt,2)&
         - Ldynf232(zi1s(1),z2s(1),z3s(1),t-delt,3)&
         + Ldynf232(zii1s(2),z2s(2),z3s(2),t-delt,4)&
         - Ldynf232(zi1s(2),z2s(2),z3s(2),t-delt,5)

    L232 = L232t - L232
    Kdynf232 = L232
  end function Kdynf232

  real(8) function Kdynf233(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt
    real(8) :: L233,L233t

    L233t = Ldynf233(zii1s(0),z2s(0),z3s(0),t,0)&
          - Ldynf233(zi1s(0),z2s(0),z3s(0),t,1)&
          + Ldynf233(zii1s(1),z2s(1),z3s(1),t,2)&
          - Ldynf233(zi1s(1),z2s(1),z3s(1),t,3)&
          + Ldynf233(zii1s(2),z2s(2),z3s(2),t,4)&
          - Ldynf233(zi1s(2),z2s(2),z3s(2),t,5)

    L233 = Ldynf233(zii1s(0),z2s(0),z3s(0),t-delt,0)&
         - Ldynf233(zi1s(0),z2s(0),z3s(0),t-delt,1)&
         + Ldynf233(zii1s(1),z2s(1),z3s(1),t-delt,2)&
         - Ldynf233(zi1s(1),z2s(1),z3s(1),t-delt,3)&
         + Ldynf233(zii1s(2),z2s(2),z3s(2),t-delt,4)&
         - Ldynf233(zi1s(2),z2s(2),z3s(2),t-delt,5)

    L233 = L233t - L233
    Kdynf233 = L233
  end function Kdynf233

  real(8) function Kdynf311(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt
    real(8) :: L311,L311t

    L311t = Ldynf311(zii1s(0),z2s(0),z3s(0),t,0)&
          - Ldynf311(zi1s(0),z2s(0),z3s(0),t,1)&
          + Ldynf311(zii1s(1),z2s(1),z3s(1),t,2)&
          - Ldynf311(zi1s(1),z2s(1),z3s(1),t,3)&
          + Ldynf311(zii1s(2),z2s(2),z3s(2),t,4)&
          - Ldynf311(zi1s(2),z2s(2),z3s(2),t,5)&
          + Hevix( z2s(0) )*Hevix( z2s(1) )*Hevix( z2s(2) )*Heviw( t-abs(z3s(0))/Ct )

    L311 = Ldynf311(zii1s(0),z2s(0),z3s(0),t-delt,0)&
         - Ldynf311(zi1s(0),z2s(0),z3s(0),t-delt,1)&
         + Ldynf311(zii1s(1),z2s(1),z3s(1),t-delt,2)&
         - Ldynf311(zi1s(1),z2s(1),z3s(1),t-delt,3)&
         + Ldynf311(zii1s(2),z2s(2),z3s(2),t-delt,4)&
         - Ldynf311(zi1s(2),z2s(2),z3s(2),t-delt,5)&
         + Hevix( z2s(0) )*Hevix( z2s(1) )*Hevix( z2s(2) )*Heviw( t-delt-abs(z3s(0))/Ct )

    L311 = L311t - L311
    Kdynf311 = L311
  end function Kdynf311

  real(8) function Kdynf312(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt
    real(8) :: L312,L312t

    L312t = Ldynf312(zii1s(0),z2s(0),z3s(0),t,0)&
          - Ldynf312(zi1s(0),z2s(0),z3s(0),t,1)&
          + Ldynf312(zii1s(1),z2s(1),z3s(1),t,2)&
          - Ldynf312(zi1s(1),z2s(1),z3s(1),t,3)&
          + Ldynf312(zii1s(2),z2s(2),z3s(2),t,4)&
          - Ldynf312(zi1s(2),z2s(2),z3s(2),t,5)
    
    L312 = Ldynf312(zii1s(0),z2s(0),z3s(0),t-delt,0)&
         - Ldynf312(zi1s(0),z2s(0),z3s(0),t-delt,1)&
         + Ldynf312(zii1s(1),z2s(1),z3s(1),t-delt,2)&
         - Ldynf312(zi1s(1),z2s(1),z3s(1),t-delt,3)&
         + Ldynf312(zii1s(2),z2s(2),z3s(2),t-delt,4)&
         - Ldynf312(zi1s(2),z2s(2),z3s(2),t-delt,5)

    L312 = L312t - L312
    Kdynf312 = L312
  end function Kdynf312
    
  real(8) function Kdynf313(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt
    real(8) :: L313,L313t

    L313t = Ldynf313(zii1s(0),z2s(0),z3s(0),t,0)&
          - Ldynf313(zi1s(0),z2s(0),z3s(0),t,1)&
          + Ldynf313(zii1s(1),z2s(1),z3s(1),t,2)&
          - Ldynf313(zi1s(1),z2s(1),z3s(1),t,3)&
          + Ldynf313(zii1s(2),z2s(2),z3s(2),t,4)&
          - Ldynf313(zi1s(2),z2s(2),z3s(2),t,5)

    L313 = Ldynf313(zii1s(0),z2s(0),z3s(0),t-delt,0)&
         - Ldynf313(zi1s(0),z2s(0),z3s(0),t-delt,1)&
         + Ldynf313(zii1s(1),z2s(1),z3s(1),t-delt,2)&
         - Ldynf313(zi1s(1),z2s(1),z3s(1),t-delt,3)&
         + Ldynf313(zii1s(2),z2s(2),z3s(2),t-delt,4)&
         - Ldynf313(zi1s(2),z2s(2),z3s(2),t-delt,5)

    L313 = L313t - L313
    Kdynf313 = L313
  end function Kdynf313

  real(8) function Kdynf321(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt
    real(8) :: L321,L321t

    L321t = Ldynf321(zii1s(0),z2s(0),z3s(0),t,0)&
          - Ldynf321(zi1s(0),z2s(0),z3s(0),t,1)&
          + Ldynf321(zii1s(1),z2s(1),z3s(1),t,2)&
          - Ldynf321(zi1s(1),z2s(1),z3s(1),t,3)&
          + Ldynf321(zii1s(2),z2s(2),z3s(2),t,4)&
          - Ldynf321(zi1s(2),z2s(2),z3s(2),t,5)

    L321 = Ldynf321(zii1s(0),z2s(0),z3s(0),t-delt,0)&
         - Ldynf321(zi1s(0),z2s(0),z3s(0),t-delt,1)&
         + Ldynf321(zii1s(1),z2s(1),z3s(1),t-delt,2)&
         - Ldynf321(zi1s(1),z2s(1),z3s(1),t-delt,3)&
         + Ldynf321(zii1s(2),z2s(2),z3s(2),t-delt,4)&
         - Ldynf321(zi1s(2),z2s(2),z3s(2),t-delt,5)
    
    L321 = L321t - L321
    Kdynf321 = L321
  end function Kdynf321

  real(8) function Kdynf322(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt
    real(8) :: L322,L322t

    L322 = Ldynf322(zii1s(0),z2s(0),z3s(0),t,0)&
         - Ldynf322(zi1s(0),z2s(0),z3s(0),t,1)&
         + Ldynf322(zii1s(1),z2s(1),z3s(1),t,2)&
         - Ldynf322(zi1s(1),z2s(1),z3s(1),t,3)&
         + Ldynf322(zii1s(2),z2s(2),z3s(2),t,4)&
         - Ldynf322(zi1s(2),z2s(2),z3s(2),t,5)&
         + Hevix( z2s(0) )*Hevix( z2s(1) )*Hevix( z2s(2) )*Heviw( t-abs(z3s(0))/Ct )

    L322 = Ldynf322(zii1s(0),z2s(0),z3s(0),t-delt,0)&
         - Ldynf322(zi1s(0),z2s(0),z3s(0),t-delt,1)&
         + Ldynf322(zii1s(1),z2s(1),z3s(1),t-delt,2)&
         - Ldynf322(zi1s(1),z2s(1),z3s(1),t-delt,3)&
         + Ldynf322(zii1s(2),z2s(2),z3s(2),t-delt,4)&
         - Ldynf322(zi1s(2),z2s(2),z3s(2),t-delt,5)&
         + Hevix( z2s(0) )*Hevix( z2s(1) )*Hevix( z2s(2) )*Heviw( t-delt-abs(z3s(0))/Ct )
    
    L322 = L322t - L322
    Kdynf322 = L322
  end function Kdynf322

  real(8) function Kdynf323(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt
    real(8) :: L323,L323t
    
    L323t = Ldynf323(zii1s(0),z2s(0),z3s(0),t,0)&
          - Ldynf323(zi1s(0),z2s(0),z3s(0),t,1)&
          + Ldynf323(zii1s(1),z2s(1),z3s(1),t,2)&
          - Ldynf323(zi1s(1),z2s(1),z3s(1),t,3)&
          + Ldynf323(zii1s(2),z2s(2),z3s(2),t,4)&
          - Ldynf323(zi1s(2),z2s(2),z3s(2),t,5)

    L323 = Ldynf323(zii1s(0),z2s(0),z3s(0),t-delt,0)&
         - Ldynf323(zi1s(0),z2s(0),z3s(0),t-delt,1)&
         + Ldynf323(zii1s(1),z2s(1),z3s(1),t-delt,2)&
         - Ldynf323(zi1s(1),z2s(1),z3s(1),t-delt,3)&
         + Ldynf323(zii1s(2),z2s(2),z3s(2),t-delt,4)&
         - Ldynf323(zi1s(2),z2s(2),z3s(2),t-delt,5)

    L323 = L323t - L323
    Kdynf323 = L323
  end function Kdynf323

  real(8) function Kdynf331(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt
    real(8) :: L331,L331t

    L331t = Ldynf331(zii1s(0),z2s(0),z3s(0),t,0)&
          - Ldynf331(zi1s(0),z2s(0),z3s(0),t,1)&
          + Ldynf331(zii1s(1),z2s(1),z3s(1),t,2)&
          - Ldynf331(zi1s(1),z2s(1),z3s(1),t,3)&
          + Ldynf331(zii1s(2),z2s(2),z3s(2),t,4)&
          - Ldynf331(zi1s(2),z2s(2),z3s(2),t,5)

    L331 = Ldynf331(zii1s(0),z2s(0),z3s(0),t-delt,0)&
         - Ldynf331(zi1s(0),z2s(0),z3s(0),t-delt,1)&
         + Ldynf331(zii1s(1),z2s(1),z3s(1),t-delt,2)&
         - Ldynf331(zi1s(1),z2s(1),z3s(1),t-delt,3)&
         + Ldynf331(zii1s(2),z2s(2),z3s(2),t-delt,4)&
         - Ldynf331(zi1s(2),z2s(2),z3s(2),t-delt,5)

    L331 = L331t - L331
    Kdynf331 = L331
  end function Kdynf331

  real(8) function Kdynf332(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt
    real(8) :: L332,L332t

    L332t = Ldynf332(zii1s(0),z2s(0),z3s(0),t,0)&
          - Ldynf332(zi1s(0),z2s(0),z3s(0),t,1)&
          + Ldynf332(zii1s(1),z2s(1),z3s(1),t,2)&
          - Ldynf332(zi1s(1),z2s(1),z3s(1),t,3)&
          + Ldynf332(zii1s(2),z2s(2),z3s(2),t,4)&
          - Ldynf332(zi1s(2),z2s(2),z3s(2),t,5)

    L332 = Ldynf332(zii1s(0),z2s(0),z3s(0),t-delt,0)&
         - Ldynf332(zi1s(0),z2s(0),z3s(0),t-delt,1)&
         + Ldynf332(zii1s(1),z2s(1),z3s(1),t-delt,2)&
         - Ldynf332(zi1s(1),z2s(1),z3s(1),t-delt,3)&
         + Ldynf332(zii1s(2),z2s(2),z3s(2),t-delt,4)&
         - Ldynf332(zi1s(2),z2s(2),z3s(2),t-delt,5)

    L332 = L332t - L332
    Kdynf332 = L332
  end function Kdynf332

  real(8) function Kdynf333(xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt)
    implicit none
    real(8),intent(in) :: xr1,xr2,xr3,xsa1,xsa2,xsb1,xsb2,xsc1,xsc2,xs3,t,delt
    real(8) :: L333,L333t
     
    L333t = Ldynf333(zii1s(0),z2s(0),z3s(0),t,0)&
          - Ldynf333(zi1s(0),z2s(0),z3s(0),t,1)&
          + Ldynf333(zii1s(1),z2s(1),z3s(1),t,2)&
          - Ldynf333(zi1s(1),z2s(1),z3s(1),t,3)&
          + Ldynf333(zii1s(2),z2s(2),z3s(2),t,4)&
          - Ldynf333(zi1s(2),z2s(2),z3s(2),t,5)&
          + 1.0_8/p*Hevix( z2s(0) )*Hevix( z2s(1) )*Hevix( z2s(2) )*Heviw( t-abs(z3s(0))/Cl )

    L333 = Ldynf333(zii1s(0),z2s(0),z3s(0),t-delt,0)&
         - Ldynf333(zi1s(0),z2s(0),z3s(0),t-delt,1)&
         + Ldynf333(zii1s(1),z2s(1),z3s(1),t-delt,2)&
         - Ldynf333(zi1s(1),z2s(1),z3s(1),t-delt,3)&
         + Ldynf333(zii1s(2),z2s(2),z3s(2),t-delt,4)&
         - Ldynf333(zi1s(2),z2s(2),z3s(2),t-delt,5)&
         + 1.0_8/p*Hevix( z2s(0) )*Hevix( z2s(1) )*Hevix( z2s(2) )*Heviw( t-delt-abs(z3s(0))/Cl )
      
    L333 = L333t - L333
    Kdynf333 = L333
  end function Kdynf333

end module kernel_module
