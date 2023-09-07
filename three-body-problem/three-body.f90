program coreografias
  implicit none
  integer i
  real(8)VXA,VXB,VXC,VYA,VYB,VYC,XA,XB,XC,YA,YB,YC,&
  gm,dm,h,N,x0
  dm=0.0d0;h=1.0d-5;gm=1.0d0+dm; N = 1.0d1
  open(9,file='condicoes.dat');
  open(10,file='A.dat');open(20,file='B.dat');open(30,file='C.dat')
  open(41,file='AB.dat');open(42,file='BC.dat');open(43,file='AC.dat')
  read(9,*)xa,ya,vxa,vya;
  read(9,*)xb,yb,vxb,vyb;
  read(9,*)xc,yc,vxc,vyc
  write(10,*)xa,ya;write(20,*)xb,yb;write(30,*)xc,yc;

  do i = 0, int(N/h)
    call rk4(VXA,VXB,VXC,VYA,VYB,VYC,XA,XB,XC,YA,YB,YC)
    write(10,*)xa,ya;write(20,*)xb,yb;write(30,*)xc,yc;
  end do
  close(9);close(10);close(20);close(30);

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
  subroutine rk4(VXA,VXB,VXC,VYA,VYB,VYC,XA,XB,XC,YA,YB,YC)
    implicit none
    real(8)VXA,VXB,VXC,VYA,VYB,VYC,XA,XB,XC,YA,YB,YC,&
    VXAI,VXBI,VXCI,VYAI,VYBI,VYCI,XAI,XBI,XCI,YAI,YBI,YCI,&
    KVXA(4),KVXB(4),KVXC(4),KVYA(4),KVYB(4),KVYC(4),&
    KXA(4),KXB(4),KXC(4),KYA(4),KYB(4),KYC(4)


    KXA(1)=H*VXA;KXB(1)=H*VXB;KXC(1)=H*VXC;
    KYA(1)=H*VYA;KYB(1)=H*VYB;KYC(1)=H*VYC;
    KVXA(1)=H*(FX_(XA,YA,XB,YB)+FX_(XA,YA,XC,YC));
    KVXB(1)=H*(FX_(XB,YB,XA,YA)+FX_(XB,YB,XC,YC));
    KVXC(1)=H*(FX_(XC,YC,XA,YA)+FX_(XC,YC,XB,YB));
    KVYA(1)=H*(FY_(XA,YA,XB,YB)+FY_(XA,YA,XC,YC));
    KVYB(1)=H*(FY_(XB,YB,XA,YA)+FY_(XB,YB,XC,YC));
    KVYC(1)=H*(FY_(XC,YC,XA,YA)+FY_(XC,YC,XB,YB));

    XAI=XA+KXA(1)/2.0D0;XBI=XB+KXB(1)/2.0D0;XCI=XC+KXC(1)/2.0D0;
    YAI=YA+KYA(1)/2.0D0;YBI=YB+KYB(1)/2.0D0;YCI=YC+KYC(1)/2.0D0;
    VXAI=VXA+KVXA(1)/2.0D0;VXBI=VXB+KVXB(1)/2.0D0;VXCI=VXC+KVXC(1)/2.0D0;
    VYAI=VYA+KVYA(1)/2.0D0;VYBI=VYB+KVYB(1)/2.0D0;VYCI=VYC+KVYC(1)/2.0D0;

    KXA(2)=H*VXAI;KXB(2)=H*VXBI;KXC(2)=H*VXCI;
    KYA(2)=H*VYAI;KYB(2)=H*VYBI;KYC(2)=H*VYCI;
    KVXA(2)=H*(FX_(XAI,YAI,XBI,YBI)+FX_(XAI,YAI,XCI,YCI));
    KVXB(2)=H*(FX_(XBI,YBI,XAI,YAI)+FX_(XBI,YBI,XCI,YCI));
    KVXC(2)=H*(FX_(XCI,YCI,XAI,YAI)+FX_(XCI,YCI,XBI,YBI));
    KVYA(2)=H*(FY_(XAI,YAI,XBI,YBI)+FY_(XAI,YAI,XCI,YCI));
    KVYB(2)=H*(FY_(XBI,YBI,XAI,YAI)+FY_(XBI,YBI,XCI,YCI));
    KVYC(2)=H*(FY_(XCI,YCI,XAI,YAI)+FY_(XCI,YCI,XBI,YBI));

    XAI=XA+KXA(2)/2.0D0;XBI=XB+KXB(2)/2.0D0;XCI=XC+KXC(2)/2.0D0;
    YAI=YA+KYA(2)/2.0D0;YBI=YB+KYB(2)/2.0D0;YCI=YC+KYC(2)/2.0D0;
    VXAI=VXA+KVXA(2)/2.0D0;VXBI=VXB+KVXB(2)/2.0D0;VXCI=VXC+KVXC(2)/2.0D0;
    VYAI=VYA+KVYA(2)/2.0D0;VYBI=VYB+KVYB(2)/2.0D0;VYCI=VYC+KVYC(2)/2.0D0;

    KXA(3)=H*VXAI;KXB(3)=H*VXBI;KXC(3)=H*VXCI;
    KYA(3)=H*VYAI;KYB(3)=H*VYBI;KYC(3)=H*VYCI;
    KVXA(3)=H*(FX_(XAI,YAI,XBI,YBI)+FX_(XAI,YAI,XCI,YCI));
    KVXB(3)=H*(FX_(XBI,YBI,XAI,YAI)+FX_(XBI,YBI,XCI,YCI));
    KVXC(3)=H*(FX_(XCI,YCI,XAI,YAI)+FX_(XCI,YCI,XBI,YBI));
    KVYA(3)=H*(FY_(XAI,YAI,XBI,YBI)+FY_(XAI,YAI,XCI,YCI));
    KVYB(3)=H*(FY_(XBI,YBI,XAI,YAI)+FY_(XBI,YBI,XCI,YCI));
    KVYC(3)=H*(FY_(XCI,YCI,XAI,YAI)+FY_(XCI,YCI,XBI,YBI));

    XAI=XA+KXA(3);XBI=XB+KXB(3);XCI=XC+KXC(3);
    YAI=YA+KYA(3);YBI=YB+KYB(3);YCI=YC+KYC(3);
    VXAI=VXA+KVXA(3);VXBI=VXB+KVXB(3);VXCI=VXC+KVXC(3);
    VYAI=VYA+KVYA(3);VYBI=VYB+KVYB(3);VYCI=VYC+KVYC(3);

    KXA(4)=H*VXAI;KXB(4)=H*VXBI;KXC(4)=H*VXCI;
    KYA(4)=H*VYAI;KYB(4)=H*VYBI;KYC(4)=H*VYCI;
    KVXA(4)=H*(FX_(XAI,YAI,XBI,YBI)+FX_(XAI,YAI,XCI,YCI));
    KVXB(4)=H*(FX_(XBI,YBI,XAI,YAI)+FX_(XBI,YBI,XCI,YCI));
    KVXC(4)=H*(FX_(XCI,YCI,XAI,YAI)+FX_(XCI,YCI,XBI,YBI));
    KVYA(4)=H*(FY_(XAI,YAI,XBI,YBI)+FY_(XAI,YAI,XCI,YCI));
    KVYB(4)=H*(FY_(XBI,YBI,XAI,YAI)+FY_(XBI,YBI,XCI,YCI));
    KVYC(4)=H*(FY_(XCI,YCI,XAI,YAI)+FY_(XCI,YCI,XBI,YBI));

    XA=XA+(1.0d0/6.0d0)*(KXA(1)+KXA(2)*2.0D0+KXA(3)*2.0D0+KXA(4));
    XB=XB+(1.0d0/6.0d0)*(KXB(1)+KXB(2)*2.0D0+KXB(3)*2.0D0+KXB(4));
    XC=XC+(1.0d0/6.0d0)*(KXC(1)+KXC(2)*2.0D0+KXC(3)*2.0D0+KXC(4));
    YA=YA+(1.0d0/6.0d0)*(KYA(1)+KYA(2)*2.0D0+KYA(3)*2.0D0+KYA(4));
    YB=YB+(1.0d0/6.0d0)*(KYB(1)+KYB(2)*2.0D0+KYB(3)*2.0D0+KYB(4));
    YC=YC+(1.0d0/6.0d0)*(KYC(1)+KYC(2)*2.0D0+KYC(3)*2.0D0+KYC(4));
    VXA=VXA+(1.0d0/6.0d0)*(KVXA(1)+KVXA(2)*2.0D0+KVXA(3)*2.0D0+KVXA(4));
    VXB=VXB+(1.0d0/6.0d0)*(KVXB(1)+KVXB(2)*2.0D0+KVXB(3)*2.0D0+KVXB(4));
    VXC=VXC+(1.0d0/6.0d0)*(KVXC(1)+KVXC(2)*2.0D0+KVXC(3)*2.0D0+KVXC(4));
    VYA=VYA+(1.0d0/6.0d0)*(KVYA(1)+KVYA(2)*2.0D0+KVYA(3)*2.0D0+KVYA(4));
    VYB=VYB+(1.0d0/6.0d0)*(KVYB(1)+KVYB(2)*2.0D0+KVYB(3)*2.0D0+KVYB(4));
    VYC=VYC+(1.0d0/6.0d0)*(KVYC(1)+KVYC(2)*2.0D0+KVYC(3)*2.0D0+KVYC(4));
  end subroutine

end program
