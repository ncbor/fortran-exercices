program log_fat
  implicit none

  !declarando as variaveis
  integer i
  real*8 resultado

  !abrindo arquivo e formatando a tabela
  open(10, file='log_fat_2-30.dat')
  write(10,*)"          n   LOG(n!)"

  !acossiando cada fatorial e log a uma linha dedicada
  do i=2, 30
    !chamando a subrotina para calcular o valor do fatorial de i
    call fatorial(i, resultado)
    !escrevendo o ln do fatorial de i em sua própria linha
    write(10,*)i, log(resultado)
  enddo


  !finalizando o programa
  close(10)

end program

 !****************iniciando subrotina ************************!
!************************************************************!

subroutine fatorial(x, resultado)
  implicit none

  !declarando as variaveis
  integer x, i
  real*8 resultado

  !calculando o fatorial de X
  resultado = 1d0
  do i=2, x
    resultado = resultado * real(i,8)
  enddo

end subroutine
