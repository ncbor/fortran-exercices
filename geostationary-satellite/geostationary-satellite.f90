program satelite
  implicit none
  integer i, N
  real(8)gmt, gml, rs, rsl,rl, xl, yl, xs, ys, vxs, vys, ts, tl, h, t, phi, phii
  gmt=6.6743d-11*5.9736d24;gml=gmt*1.23d-2;rl=3.844d8;ts=8.64d4;tl=27.32d0*ts;h=1.0d-2
  xs=4.2244d7;ys=0.0d0;vxs=0.0d0;vys=2*acos(-1.0d0)*xs/ts;
  open(10,file='p_satelite.dat');open(20,file='p_lua.dat')
  open(41,file='r_satelite.dat');open(42,file='dr_satelite.dat')
  open(51,file='phi_satelite.dat');open(52,file='dphi_satelite.dat')

  write(10,*)xs, ys;write(15,*)vxs, vys;write(20,*)rl, 0.0d0
  i=0;xl = rl*cos(2*acos(-1.0d0)*i/tl);yl = rl*sin(2*acos(-1.0d0)*i/tl)
  rs=rs_(xs,ys);rsl=rsl_(xs,ys,xl,yl)

  do i=1,int(4.5d0*tl)
    t=real(i,8)/ts
    call rungek4(xs, vxs, ys, vys, rs, rsl)
    xl = rl*cos(2*acos(-1.0d0)*i/tl)
    yl = rl*sin(2*acos(-1.0d0)*i/tl)
    phi=mod(atan2(ys,xs)+acos(-1.0d0)*2.0d0,acos(-1.0d0)*2.0d0)
    phii =mod(acos(-1.0d0)*2.0d0*real(i,8)/ts,acos(-1.0d0)*2.0d0)
    write(10,*)xs, ys;write(20,*)xl, yl
    write(41,*)t,rs_(xs,ys);write(42,*)t,rs_(xs,ys)-4.2244d7;
    write(51,*)t,phi;
    if(phi.lt.phii)then
      write(52,*)t,-(phi-phii+acos(-1.0d0)*2)
    else
      write(52,*)t,phii-phi
    endif
  end do

  close(10);close(20);close(51);close(52);close(41);close(42)
 !
contains
  !
  real(8)function rs_(xs,ys);implicit none;real(8)xs,ys
    rs_=sqrt(xs**2.0d0 + ys**2.0d0)
  end function
  real(8)function rsl_(xs,ys,xl,yl);implicit none;real(8)xs,ys,xl,yl
    rsl_ = sqrt((xs-xl)**2.0d0 + (ys-yl)**2.0d0)
  end function
  real(8)function fx(xs, ys, xl, yl);implicit none;real(8)xs,ys,xl,yl
    fx=-gml*(xs-xl)/rsl_(xs, ys, xl, yl)**3.0d0 -gmt*xs/rs_(xs,ys)**3.0d0
  end function
  real(8)function fy(xs, ys, xl, yl);implicit none;real(8)xs,ys,xl,yl
    fy=-gml*(ys-yl)/rsl_(xs, ys, xl, yl)**3.0d0 -gmt*ys/rs_(xs,ys)**3.0d0
  end function
  subroutine rungek4(xs, vxs, ys, vys, rs, rsl)
    implicit none
    real(8)xs, vxs, ys, vys, rs, rsl, kx(4), kxv(4), ky(4), kyv(4)

    kx(1)=h*vxs;ky(1)=h*vys
    kxv(1)=h*fx(xs, ys, xl, yl);kyv(1)=h*fy(xs, ys, xl, yl)

    kx(2)=h*(vxs+kxv(1)/2.0d0);ky(2)=h*(vys+kyv(1)/2.0d0)
    kxv(2)=h*fx(xs+kx(1)/2.0d0, ys+ky(1)/2.0d0, xl, yl);kyv(2)=h*fy(xs+kx(1)/2.0d0, ys+ky(1)/2.0d0, xl, yl)

    kx(3)=h*(vxs+kxv(2)/2.0d0);ky(3)=h*(vys+kyv(2)/2.0d0)
    kxv(3)=h*fx(xs+kx(2)/2.0d0, ys+ky(2)/2.0d0, xl, yl);kyv(3)=h*fy(xs+kx(2)/2.0d0, ys+ky(2)/2.0d0, xl, yl)

    kx(4)=h*(vxs+kxv(3));ky(4)=h*(vys+kyv(3))
    kxv(4)=h*fx(xs+kx(3), ys+ky(3), xl, yl);kyv(4)=h*fy(xs+kx(3), ys+ky(3), xl, yl)

    xs = xs + (kx(1)+kx(2)*2.0d0+kx(3)*2.0d0+kx(4))/6.0d0
    ys = ys + (ky(1)+ky(2)*2.0d0+ky(3)*2.0d0+ky(4))/6.0d0
    vxs = vxs +(kxv(1)+kxv(2)*2.0d0+kxv(3)*2.0d0+kxv(4))/6.0d0
    vys = vys +(kyv(1)+kyv(2)*2.0d0+kyv(3)*2.0d0+kyv(4))/6.0d0

    rs = rs_(xs,ys)
    rsl= rsl_(xs,ys,vxs,vys)
  end subroutine

end program
