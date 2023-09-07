program bhaskara
  implicit none

  !declarando as variaveis
  real(8) a, b, c, delta, x1, x2
  integer ios

  !recebendo valores do usuario
  print *, "Insira os coeficientes a, b, e c, de uma equação polinomial do segundo grau:"
  read(*,*, iostat=ios) a, b, c

  !admnistrando erros no input
  if (ios.ne.0) then

    print *, "Isso não são numeros reais."

  else

    !calculando delta
    delta = b**2 - 4*a*c

    !fazendo o bhaskara de fato com suas 3 possibilidaes
    if (delta>0) then !delta positivo

      x1 = ((-b + sqrt(delta))/(2*a))
      x2 = ((-b - sqrt(delta))/(2*a))
      print *, "As raizes reais são:"
      print *, x1, x2

    else if (delta<0) then !delta negativo

      print *, "Não há raizes reais."

    else

      !delta igual a 0
      x1 = -b/(2*a)
      print *, "A equação possui apenas uma raiz real:"
      print *, x1

    end if
  end if

end program bhaskara
