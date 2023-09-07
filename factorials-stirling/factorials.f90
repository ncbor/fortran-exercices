program fatoriais
  implicit none

  !declarando as variaveis
  integer i
  real*8 aux

  !atribuindo valor a aux para poder ser utilizado na linha 20
  aux = 1d0

  !criando arquivo com o nome em questão
  open(10, file='fatoriais_1-20.dat')

  !esrevendo o que representa cada coluna
  write(10,*)"          n    n!"

  !para tratar cada N em sua devia linha
  do i=1, 30
    !fazendo o cálculo de cada fatorial
    aux = aux*real(i,8)
    !escrevendo cada valor em sua devida linha
    write(10,*)i, aux
  enddo

  !finalizando o programa

  close(10)
end program
