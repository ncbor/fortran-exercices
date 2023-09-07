program derivada_numerica
  implicit none

  integer ios
  real(8) df, dt, dsim, d2sim, h, ef, et, esim, e2sim, da, d2a

  open(10,file="input.dat")
  open(20,file="derivadas.dat")
  open(30,file="erros.dat")
  open(40,file="logdf.dat")
  open(50,file="logdsim.dat")
  open(60,file="logd2sim.dat")

  da = ( 0.5d0 * exp(0.5d0) * sin(1.0d0/3.0d0) )      + ( cos(1.0d0/3.0d0) * exp(0.5d0) ) / 3.0d0
  d2a =( exp(0.5d0)*sin(1.0d0/3.0d0)*(5.0d0/36.0d0) ) + ( cos(1.0d0/3.0d0) * exp(0.5d0) ) / 3.0d0

  do
    read(10,*,iostat=ios) h
    if (ios/=0) EXIT
    df = ((F(1.0d0+h) - F(1.0d0)))/h
    dt = (F(1.0d0) - F(1.0d0-h))/h
    dsim = 0.5d0 * (df + dt)
    d2sim = (F(1.0d0+h) - 2.0d0*F(1.0d0) + F(1.0d0-h))/(h**2.0d0)
    write(20,*)df, dt, dsim, d2sim
    ef = abs(df - da)
    et = abs(dt - da)
    esim = abs(dsim - da)
    e2sim = abs(d2sim - d2a)
    write(30,*) ef, et, esim, e2sim
    write(40,*) dlog10(h), dlog10(ef)
    write(50,*) dlog10(h), dlog10(esim)
    write(60,*) dlog10(h), dlog10(e2sim)
  end do

  close(10)
  close(20)
  close(30)
  close(40)
  close(50)
  close(60)

contains

real(8) function F(x)
  implicit none

  real(8) x

  F = EXP(X/2.0d0)*SIN(X/3.0d0)

end function

end program
