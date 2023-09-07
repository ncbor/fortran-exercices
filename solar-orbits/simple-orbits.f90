program orbitas
  implicit none
  character(5) xdat, vdat
  character(6) nome
  integer i, j, l, N
  real(8)x, y, vx, vy,r, a, count, G, m, gms, h
  gms=4*(acos(-1.0d0))**2;h=1.0d-5;N=1.0d5
  xdat='x.dat';vdat='v.dat';G=2.0d-30*(acos(-1.0d0))**2
  open(17,file='kepleres.dat')
  write(17,*)'Planeta   Semieixo Maior(UA)    VelocidadeInicial (UA/ano)     Periodo (anos)             Razão Kepler'

  open(9,file='planetas.dat')
  do l = 10, 80, 10
    if (l.gt.45)then
      h=1.0d-3;N=1.0d7
    end if
    read(9,*)nome, m, a
    open(l,file=nome//xdat);open(l+1,file=nome//vdat)
    y=0.0d0;vx=0.0d0;
    x=a
    vy=sqrt(gms/a)
    write(l,*)x, y
    write(l+1,*)vx, vy

    do i=1,N
      r = sqrt(x**2 + y**2)
      call rungek4(x, vx, y, vy, r)
      write(l,*)x, y
      write(l+1,*)vx, vy
      if(y.lt.0)then
        do j=1,N
          r = sqrt(x**2 + y**2)
          call rungek4(x, vx, y, vy, r)
          write(l,*)x, y
          write(l+1,*)vx, vy
          if(y.gt.0)then
            write(17,*)nome, a, sqrt(gms/a), (i+j)*h, ((i+j)*h)**2/a**3
            exit
          end if
        end do
        exit
      end if
    end do
  end do
contains

  subroutine rungek4(x, vx, y, vy, r)
    implicit none
    real(8)x, vx, y, vy, r, kx(4), kvx(4), ky(4), kvy(4)
    kx(1)=h*vx
    ky(1)=h*vy
    kvx(1)=-h*gms*x/r**3
    kvy(1)=-h*gms*y/r**3
    kx(2)=h*(vx+kvx(1)/2.0d0)
    ky(2)=h*(vy+kvy(1)/2.0d0)
    r=sqrt((x+kx(1)/2.0d0)**2+(y+ky(1)/2.0d0)**2)
    kvx(2)=-h*gms*(x+kx(1)/2.0d0)/r**3
    kvy(2)=-h*gms*(y+ky(1)/2.0d0)/r**3
    kx(3)=h*(vx+kvx(2)/2.0d0)
    ky(3)=h*(vy+kvy(2)/2.0d0)
    r=sqrt((x+kx(2)/2.0d0)**2+(y+ky(2)/2.0d0)**2)
    kvx(3)=-h*gms*(x+kx(2)/2.0d0)/r**3
    kvy(3)=-h*gms*(y+ky(2)/2.0d0)/r**3
    kx(4)=h*(vx+kvx(3))
    ky(4)=h*(vy+kvy(3))
    r=sqrt((x+kx(3))**2+(y+ky(3))**2)
    kvx(4)=-h*gms*(x+kx(3))/r**3
    kvy(4)=-h*gms*(y+ky(3))/r**3
    x = x + (kx(1)+kx(2)*2.0d0+kx(3)*2.0d0+kx(4))/6.0d0
    y = y + (ky(1)+ky(2)*2.0d0+ky(3)*2.0d0+ky(4))/6.0d0
    vx = vx +(kvx(1)+kvx(2)*2.0d0+kvx(3)*2.0d0+kvx(4))/6.0d0
    vy = vy +(kvy(1)+kvy(2)*2.0d0+kvy(3)*2.0d0+kvy(4))/6.0d0
  end subroutine
  
end program
