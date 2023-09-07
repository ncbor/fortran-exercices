program comparador
  implicit none

  !declarando variaveis
  integer n, i
  real(8) s, nfat, logn

  !abrindo arquivos de leitura
  open(10, file="stirling_2-30.dat")
  open(20, file="log_fat_2-30.dat")
  !criando arquivo de output
  open(30, file="tabela_final.dat")

  !pulando a linha sem valores numéricos de cada arquivo lido
  read(10,*)
  read(20,*)

  !formatando o arquivo de saída para a tabela em questão
  write(30,"(A)")" N      ln(n)                   S                      (ln(n) - S) / ln(n)"

  !uma linha para cada N
  do i = 1, 29

    !lendo os arquivos e atribuindo as variaveis
    read(10,*)n, s, nfat
    read(20,*)n, logn

    !escrevendo no arquivo de saída
    write(30,'(I2, A, F19.16,A, F19.16, A,F19.16)')n,"     ", logn,"     ", s, "     ",((logn-s)/logn)
  enddo

  !fechando os arquivos
  close(10); close(20); close(30)

endprogram
