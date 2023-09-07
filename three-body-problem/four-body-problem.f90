program coreografias_quatro
  implicit none
  integer i
  real(8)VXA,VXB,VXC,VXD,VYA,VYB,VYC,VYD,XA,XB,XC,XD,YA,YB,YC,YD,&
  gm,dm,h,N
  dm=0.0d0;h=1.0d-5;gm=1.0d0+dm; N = 1.0d1
  open(9,file='condicoes.dat');
  open(10,file='A.dat');open(20,file='B.dat');open(30,file='C.dat');open(40,file='D.dat')
  read(9,*)xa,ya,vxa,vya;
  read(9,*)xb,yb,vxb,vyb;
  read(9,*)xc,yc,vxc,vyc
  read(9,*)xd,yd,vxd,vyd
  write(10,*)xa,ya;write(20,*)xb,yb;write(30,*)xc,yc;write(40,*)XD,YD;

  do i = 0, int(N/h)
    call rk4(VXA,VXB,VXC,VYA,VYB,VYC,XA,XB,XC,YA,YB,YC,VXD,VYD,XD,YD)
    write(10,*)xa,ya;write(20,*)xb,yb;write(30,*)xc,yc;write(40,*)XD,YD;
  end do
  close(9);close(10);close(20);close(30);close(40)

contains
  real(8)function fx_(x1,y1,x2,y2)
    implicit none
    real(8) r,x1,y1,x2,y2
    r = sqrt((x1-x2)**2.0d0 + (y1-y2)**2.0d0)
    fx_ = -gm * ((x1-x2)/((r)**3.0d0))
  end function
  real(8)function fy_(x1,y1,x2,y2)
    implicit none
    real(8) r,x1,y1,x2,y2
    r = sqrt((x1-x2)**2.0d0 + (y1-y2)**2.0d0)
    fy_ = -gm * ((y1-y2)/((r)**3.0d0))
  end function
  subroutine rk4(VXA,VXB,VXC,VYA,VYB,VYC,XA,XB,XC,YA,YB,YC,VXD,VYD,XD,YD)
    implicit none
    real(8)VXA,VXB,VXC,VYA,VYB,VYC,XA,XB,XC,YA,YB,YC,XD,YD,VXD,VYD,&
    VXAI,VXBI,VXCI,VYAI,VYBI,VYCI,XAI,XBI,XCI,YAI,YBI,YCI,YDI,XDI,VXDI,VYDI,&
    KVXA(4),KVXB(4),KVXC(4),KVXD(4),KVYA(4),KVYB(4),KVYC(4),KVYD(4),&
    KXA(4),KXB(4),KXC(4),KXD(4),KYA(4),KYB(4),KYC(4),KYD(4)


    KXA(1)=H*VXA;KXB(1)=H*VXB;KXC(1)=H*VXC;KXD(1)=H*VXD;
    KYA(1)=H*VYA;KYB(1)=H*VYB;KYC(1)=H*VYC;KYD(1)=H*VYD;
    KVXA(1)=H*(FX_(XA,YA,XB,YB)+FX_(XA,YA,XC,YC)+FX_(XA,YA,XD,YD));
    KVXB(1)=H*(FX_(XB,YB,XA,YA)+FX_(XB,YB,XC,YC)+FX_(XB,YB,XD,YD));
    KVXC(1)=H*(FX_(XC,YC,XA,YA)+FX_(XC,YC,XB,YB)+FX_(XC,YC,XD,YD));
    KVXD(1)=H*(FX_(XD,YD,XA,YA)+FX_(XD,YD,XB,YB)+FX_(XD,YD,XC,YC));
    KVYA(1)=H*(FY_(XA,YA,XB,YB)+FY_(XA,YA,XC,YC)+FY_(XA,YA,XD,YD));
    KVYB(1)=H*(FY_(XB,YB,XA,YA)+FY_(XB,YB,XC,YC)+FY_(XB,YB,XD,YD));
    KVYC(1)=H*(FY_(XC,YC,XA,YA)+FY_(XC,YC,XB,YB)+FY_(XC,YC,XD,YD));
    KVYD(1)=H*(FY_(XD,YD,XA,YA)+FY_(XD,YD,XB,YB)+FY_(XD,YD,XC,YC));

    XAI=XA+KXA(1)/2.0D0;XBI=XB+KXB(1)/2.0D0;XCI=XC+KXC(1)/2.0D0;XDI=XD+KXD(1)/2.0D0;
    YAI=YA+KYA(1)/2.0D0;YBI=YB+KYB(1)/2.0D0;YCI=YC+KYC(1)/2.0D0;YDI=YD+KYD(1)/2.0D0;
    VXAI=VXA+KVXA(1)/2.0D0;VXBI=VXB+KVXB(1)/2.0D0;VXCI=VXC+KVXC(1)/2.0D0;VXDI=VXD+KVXD(1)/2.0D0;
    VYAI=VYA+KVYA(1)/2.0D0;VYBI=VYB+KVYB(1)/2.0D0;VYCI=VYC+KVYC(1)/2.0D0;VYDI=VYD+KVYD(1)/2.0D0;

    KXA(2)=H*VXA;KXB(2)=H*VXB;KXC(2)=H*VXC;KXD(2)=H*VXD;
    KYA(2)=H*VYA;KYB(2)=H*VYB;KYC(2)=H*VYC;KYD(2)=H*VYD;
    KVXA(2)=H*(FX_(XA,YA,XB,YB)+FX_(XA,YA,XC,YC)+FX_(XA,YA,XD,YD));
    KVXB(2)=H*(FX_(XB,YB,XA,YA)+FX_(XB,YB,XC,YC)+FX_(XB,YB,XD,YD));
    KVXC(2)=H*(FX_(XC,YC,XA,YA)+FX_(XC,YC,XB,YB)+FX_(XC,YC,XD,YD));
    KVXD(2)=H*(FX_(XD,YD,XA,YA)+FX_(XD,YD,XB,YB)+FX_(XD,YD,XC,YC));
    KVYA(2)=H*(FY_(XA,YA,XB,YB)+FY_(XA,YA,XC,YC)+FY_(XA,YA,XD,YD));
    KVYB(2)=H*(FY_(XB,YB,XA,YA)+FY_(XB,YB,XC,YC)+FY_(XB,YB,XD,YD));
    KVYC(2)=H*(FY_(XC,YC,XA,YA)+FY_(XC,YC,XB,YB)+FY_(XC,YC,XD,YD));
    KVYD(2)=H*(FY_(XD,YD,XA,YA)+FY_(XD,YD,XB,YB)+FY_(XD,YD,XC,YC));

    XAI=XA+KXA(2)/2.0D0;XBI=XB+KXB(2)/2.0D0;XCI=XC+KXC(2)/2.0D0;XDI=XD+KXD(2)/2.0D0;
    YAI=YA+KYA(2)/2.0D0;YBI=YB+KYB(2)/2.0D0;YCI=YC+KYC(2)/2.0D0;YDI=YD+KYD(2)/2.0D0;
    VXAI=VXA+KVXA(2)/2.0D0;VXBI=VXB+KVXB(2)/2.0D0;VXCI=VXC+KVXC(2)/2.0D0;VXDI=VXD+KVXD(2)/2.0D0;
    VYAI=VYA+KVYA(2)/2.0D0;VYBI=VYB+KVYB(2)/2.0D0;VYCI=VYC+KVYC(2)/2.0D0;VYDI=VYD+KVYD(2)/2.0D0;

    KXA(3)=H*VXA;KXB(3)=H*VXB;KXC(3)=H*VXC;KXD(3)=H*VXD;
    KYA(3)=H*VYA;KYB(3)=H*VYB;KYC(3)=H*VYC;KYD(3)=H*VYD;
    KVXA(3)=H*(FX_(XA,YA,XB,YB)+FX_(XA,YA,XC,YC)+FX_(XA,YA,XD,YD));
    KVXB(3)=H*(FX_(XB,YB,XA,YA)+FX_(XB,YB,XC,YC)+FX_(XB,YB,XD,YD));
    KVXC(3)=H*(FX_(XC,YC,XA,YA)+FX_(XC,YC,XB,YB)+FX_(XC,YC,XD,YD));
    KVXD(3)=H*(FX_(XD,YD,XA,YA)+FX_(XD,YD,XB,YB)+FX_(XD,YD,XC,YC));
    KVYA(3)=H*(FY_(XA,YA,XB,YB)+FY_(XA,YA,XC,YC)+FY_(XA,YA,XD,YD));
    KVYB(3)=H*(FY_(XB,YB,XA,YA)+FY_(XB,YB,XC,YC)+FY_(XB,YB,XD,YD));
    KVYC(3)=H*(FY_(XC,YC,XA,YA)+FY_(XC,YC,XB,YB)+FY_(XC,YC,XD,YD));
    KVYD(3)=H*(FY_(XD,YD,XA,YA)+FY_(XD,YD,XB,YB)+FY_(XD,YD,XC,YC));

    XAI=XA+KXA(3);XBI=XB+KXB(3);XCI=XC+KXC(3);XDI=XD+KXD(3);
    YAI=YA+KYA(3);YBI=YB+KYB(3);YCI=YC+KYC(3);YDI=YD+KYD(3);
    VXAI=VXA+KVXA(3);VXBI=VXB+KVXB(3);VXCI=VXC+KVXC(3);VXDI=VXD+KVXD(3);
    VYAI=VYA+KVYA(3);VYBI=VYB+KVYB(3);VYCI=VYC+KVYC(3);VYDI=VYD+KVYD(3);

    KXA(4)=H*VXA;KXB(4)=H*VXB;KXC(4)=H*VXC;KXD(4)=H*VXD;
    KYA(4)=H*VYA;KYB(4)=H*VYB;KYC(4)=H*VYC;KYD(4)=H*VYD;
    KVXA(4)=H*(FX_(XA,YA,XB,YB)+FX_(XA,YA,XC,YC)+FX_(XA,YA,XD,YD));
    KVXB(4)=H*(FX_(XB,YB,XA,YA)+FX_(XB,YB,XC,YC)+FX_(XB,YB,XD,YD));
    KVXC(4)=H*(FX_(XC,YC,XA,YA)+FX_(XC,YC,XB,YB)+FX_(XC,YC,XD,YD));
    KVXD(4)=H*(FX_(XD,YD,XA,YA)+FX_(XD,YD,XB,YB)+FX_(XD,YD,XC,YC));
    KVYA(4)=H*(FY_(XA,YA,XB,YB)+FY_(XA,YA,XC,YC)+FY_(XA,YA,XD,YD));
    KVYB(4)=H*(FY_(XB,YB,XA,YA)+FY_(XB,YB,XC,YC)+FY_(XB,YB,XD,YD));
    KVYC(4)=H*(FY_(XC,YC,XA,YA)+FY_(XC,YC,XB,YB)+FY_(XC,YC,XD,YD));
    KVYD(4)=H*(FY_(XD,YD,XA,YA)+FY_(XD,YD,XB,YB)+FY_(XD,YD,XC,YC));

    XA=XA+(1.0d0/6.0d0)*(KXA(1)+KXA(2)*2.0D0+KXA(3)*2.0D0+KXA(4));
    XB=XB+(1.0d0/6.0d0)*(KXB(1)+KXB(2)*2.0D0+KXB(3)*2.0D0+KXB(4));
    XC=XC+(1.0d0/6.0d0)*(KXC(1)+KXC(2)*2.0D0+KXC(3)*2.0D0+KXC(4));
    XD=XD+(1.0d0/6.0d0)*(KXD(1)+KXD(2)*2.0D0+KXD(3)*2.0D0+KXD(4));
    YA=YA+(1.0d0/6.0d0)*(KYA(1)+KYA(2)*2.0D0+KYA(3)*2.0D0+KYA(4));
    YB=YB+(1.0d0/6.0d0)*(KYB(1)+KYB(2)*2.0D0+KYB(3)*2.0D0+KYB(4));
    YC=YC+(1.0d0/6.0d0)*(KYC(1)+KYC(2)*2.0D0+KYC(3)*2.0D0+KYC(4));
    YD=YD+(1.0d0/6.0d0)*(KYD(1)+KYD(2)*2.0D0+KYD(3)*2.0D0+KYD(4));
    VXA=VXA+(1.0d0/6.0d0)*(KVXA(1)+KVXA(2)*2.0D0+KVXA(3)*2.0D0+KVXA(4));
    VXB=VXB+(1.0d0/6.0d0)*(KVXB(1)+KVXB(2)*2.0D0+KVXB(3)*2.0D0+KVXB(4));
    VXC=VXC+(1.0d0/6.0d0)*(KVXC(1)+KVXC(2)*2.0D0+KVXC(3)*2.0D0+KVXC(4));
    VXD=VXD+(1.0d0/6.0d0)*(KVXD(1)+KVXD(2)*2.0D0+KVXD(3)*2.0D0+KVXD(4));
    VYA=VYA+(1.0d0/6.0d0)*(KVYA(1)+KVYA(2)*2.0D0+KVYA(3)*2.0D0+KVYA(4));
    VYB=VYB+(1.0d0/6.0d0)*(KVYB(1)+KVYB(2)*2.0D0+KVYB(3)*2.0D0+KVYB(4));
    VYC=VYC+(1.0d0/6.0d0)*(KVYC(1)+KVYC(2)*2.0D0+KVYC(3)*2.0D0+KVYC(4));
    VYD=VYD+(1.0d0/6.0d0)*(KVYD(1)+KVYD(2)*2.0D0+KVYD(3)*2.0D0+KVYD(4));
  end subroutine

end program
