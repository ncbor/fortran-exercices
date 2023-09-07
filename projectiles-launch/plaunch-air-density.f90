program lancamento_dens
  implicit none
  integer i
  real(8),  parameter :: pi = 4 * atan (1.0d0)
  real(8) x, y, vr, vx, vy, delta, theta, g, t, res, coef, dens
  delta=0.01d0;vr=700.0d0;g=9.8d0;coef=4.0d-5;
  res= coef*sqrt(vr**2)

  open(10,file='dados1.dat');open(20,file='dados2.dat');open(30,file='dados3.dat')
  open(40,file='dados4.dat');open(50,file='dados5.dat');open(60,file='dados6MAX.dat')

  write(*,*)'       Theta(rad)                 Tempo(s)'
  do i = 10, 60, 10
    t=delta
    x=0.0d0;y=0.0d0

    if(i.eq.10)then
      theta=pi/2.5d0
    else if(i.eq.20)then
      theta=pi/3.0d0
    else if(i.eq.30)then
      theta=pi/4.0d0
    else if(i.eq.40)then
      theta=pi/6.0d0
    else if(i.eq.50)then
      theta=pi/12.0d0
    else if(i.eq.60)then
      theta=0.76000632679224944d0
    end if

    write(i,*)x, y

    vx = vr*cos(theta); vy = vr*sin(theta)
    x = vx*delta; y = vy*delta

      do
        dens=(1-6.5d-3*y/300d0)**2.5d0
        t=t+delta
        res = dens*coef*sqrt((vx**2)+(vy**2))
        vx = vx - res*delta*vx
        vy = vy - g*delta - res*delta*vy
        x = x + vx*delta
        y = y + vy*delta

        write(i,*)x, y
        if(y.lt.0)exit
      end do
      write(*,*)theta, t
  end do

  close(10);close(20);close(30);close(40);close(50);close(60)

end program
