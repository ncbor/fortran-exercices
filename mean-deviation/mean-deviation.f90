program media_desvio
  implicit none

  !Declarando as variáveis
  real(8)  aux, media, media_da, desvio, DesvioAprox
  integer ios, n, i
  character(40) arg

  i = 0; n = 0; ios  = 0; media = 0; media_da = 0; desvio = 0; DesvioAprox = 0

  !Para que o input do arquivo seja através do terminal
  call getarg(1, arg)
  open(10, file= arg)

  !Leitura da quantidade de linhas
  do
    read(10,*,iostat=ios) aux
    if (ios/=0) EXIT
    n = n + 1

    !comecando a media
    media = media + aux

    !comecando a media para o desvio aproximado
    media_da = media_da + (aux**2)

  end do

  !De fato calculando as médias
  media = media / real(n,8)

  media_da = media_da / real(n,8)

  !Terminando DesvioAprox
  DesvioAprox = sqrt(media_da - (media**2))

  !Necessário para ler o arquivo de novo
  rewind(10)

  !Cálculo do desvio padrão
  do i=1, N
    read(10,*) aux
    desvio = desvio + (aux - media) ** 2
  end do

  desvio = sqrt(desvio/(n-1))

  !Fechando o arquivo utilizado
  close(10)

  !Imprimindo os resultados no terminal
  print*, "Média              ", media
  print*, "Desvio Padrão      ", desvio
  print*, "Desvio Aproximado  ", DesvioAprox

end program
