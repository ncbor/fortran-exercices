program comparacao_stirling
  implicit none

  !declarando as variaveis e o valor de pi
  integer j
  real*8 i, s, fat
  real(8),  parameter :: pi = 4 * atan (1.0d0)

  !criando o arquivo e formatando a tabela
  open(10,file ="stirling_2-30.dat")
  write(10,*)"          n   S                         n!"

  !trabalhando com cada N de cada vez
  do j=2, 30
    !aqui foi utilizado "i" para facilitar a escrita
    i = real(j,8)

    !calculo do stirling de fato
    s = i*log(i) - i + ( ( log(2*pi*i) ) / 2.0d0 )

    !calculo do fatorial através do log dele
    fat = exp(s)
    !escrevendo os valores em suas respectivas linhas
    write(10,*)j, s, fat
  enddo
  !finalizando o programa
  close(10)
end program
