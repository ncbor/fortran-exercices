program lancamento_vacuo
  implicit none
  integer i
  real(8),  parameter :: pi = 4 * atan (1.0d0)
  real(8) x, y, vr, vx, vy, delta, theta, g, t, alcance
  delta=0.01d0;vr=700.0d0;g=9.8d0

  open(10,file='dados1.dat');open(20,file='dados2.dat');open(30,file='dados3.dat')
  open(40,file='dados4.dat');open(50,file='dados5.dat')

  write(*,*)
  write(*,*)"    Theta(rad)                Alcançe(m)                Tempo(s)"
  do i = 10, 50, 10
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
    end if

    alcance = (vr**2 * sin(2*theta))/g
    t = alcance/(vr*cos(theta))

    write(*,*)theta, alcance, t
    write(i,*)x, y

    vx = vr*cos(theta); vy = vr*sin(theta) - g*delta
    x = vx*delta; y = vy*delta

      do
        x = x + vx*delta
        vy = vy - g*delta
        y = y + vy*delta

        if(y.lt.0)exit
        write(i,*)x, y
      end do

      write(i,*)alcance, 0.0d0
  end do
  write(*,*)
  close(10);close(20);close(30);close(40);close(50)

end program
