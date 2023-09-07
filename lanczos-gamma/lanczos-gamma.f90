program gamma_function
  implicit none

  !declarando as variaveis
  integer i
  real(8) x

  !criando e formatando o arquivo do output
  open(10, file="fatoriais_gamma.dat")
  write(10, *)"          x          x-1    Gamma(x) = x-1!"

  !usando do para mexer em um numero de cada vez, em ordem
  do i=5, 20
    x = real(i,8) !mantenho I inteiro mas utilizo X real na função para formatação
    write(10,*)i+1, i, lanczos(x) !escrevendo no arquivo de output
  enddo

  !computando os 2 ultimos casos solicitados pelo exercício
  x = 1.5
  write(10,'(/,A, F3.1, A, F18.16)')"           ",x, "             ",lanczos(x-1)
  x = 2.5
  write(10,'(A, F3.1, A, F18.16)')"           ",x,"             ", lanczos(x-1)

  close(10)

contains

  !___________________________________________________________!
 !*   CRIANDO FUNÇÃO PARA O CÁLCULO  DA FORMULA DE LANCZOS  *!
!‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾!

real(8) function lanczos(x)
  implicit none

  real(8),  parameter :: PI  = 4 * atan (1.0d0) !Definindo PI

  !declarando as variaveis
  integer aux
  real(8) gamma, coeficiente, x
  real(8), dimension(0:6) :: c

  !definindo as constantes
  c(0)=  1.000000000190015; c(1)=  76.18009172947146
  c(2)= -86.50532032941677; c(3)=  24.01409824083091
  c(4)= -1.231739572450155; c(5)=  0.1208650973866179D-02; c(6)= -0.5395239384953D-05

  coeficiente = c(0)!para poupar um IF na hora de somar os coeficientes
  gamma = 5.0 !definindo o gamma de acordo com o requisitado

  !calculando a soma dos c(n)/(x+n)
  do aux = 1,6

    coeficiente = coeficiente + ((c(aux))/(x+real(aux, 8)))
  end do

  !de fato executando a fórmula de Lanczos
  lanczos = ((x+gamma+0.5d0)**(x+0.5d0)) * exp(-(x+gamma+0.5d0)) * sqrt(2*PI) * coeficiente

end function
!!!!!!!!!!!
!!!!!!!!!!!
end program
