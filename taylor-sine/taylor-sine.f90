program seno_staylor
  implicit none

  !Declarando as variaveis
  real(8),  parameter :: pi = 4 * atan (1.0d0)
  real(8) x, res, comp
  integer i, n

  !Abrindo e formatando o arquivo de saída
  open(10, file='senos.dat')
  write(10,*)"  Ângulo X                sin(x)                    taylor_sin(x)          abs(taylor_sin - sin)/ sin      exp ordem"
  write(10,'(/)')

  !Abrindo um loop para executar o programa de fato
  do i=1, 10
    !Gerando um número aleatório entre 0 e 2pi
    call Random_Number(x)
    x = x * 2 * pi
    !Calculando o seno por taylor e a ordem de expansão N
    call seno_taylor(x, res, n)
    !Comparando o valor calculado com o valor real de seno
    comp = abs(res-sin(x))/sin(x)
    !Escrevendo os dados no arquivo de saída
    write(10,*) x, sin(x), res, comp, n

  enddo

  close(10)
end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!SUBROUTINE PARA O CÁLCULO DE FATO!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine seno_taylor(angulo, resultado, ordem)
  implicit none

  !Declarando as variaveis
  real(8) angulo, soma, taylor, resultado
  integer i, ordem, k

  soma = 0; i = 0

  !Loop para cálculo de cada termo
  do
    !Ordem atual
    k = 1 + 2*i
    !Cálculo do termo em questão
    taylor = ((-1)**i) * (angulo**(k)) / real( fatorial( k ) , 16 )
    !Limite para a precisão
    if (abs(taylor).lt.1d-06) EXIT
    !Fazendo o somatório
    soma = soma + taylor
    !Próximo índice
    i = i+1

  end do

  !Atribuindo valores aos argumentos da subroutine
  resultado = soma
  ordem = k

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!PARA CALCULAR O FATORIAL!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real function fatorial(x)
    implicit none

  integer i, x
  fatorial = 1.0d0
  do i = 2, x
    fatorial = fatorial * i

  end do

  end function
!!!
!!!!
!!!!!

end subroutine
