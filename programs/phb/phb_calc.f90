program bbb
implicit none
real*8, allocatable :: t(:), c(:)
real*8 :: read_c, dum
integer :: i, j, nline, n_file
character(50) :: argv, acp, dnr
character(50) :: infile, outfile

call getarg(1,acp)
call getarg(2,dnr)
call getarg(3,argv)
read(argv,*) n_file

write(argv, '(I4.4)') 1
infile = 'Data_phb_'//trim(acp)//'_'//trim(dnr)//'.'//trim(argv)

open(15,file=infile,status='old')
nline = 0
do
  read(15,*,end=888)
  nline = nline + 1
end do
888 close(15)

allocate(c(nline))
allocate(t(nline))

c(:) = 0d0

do i = 1, n_file
  write(argv,'(I4.4)') i 
  infile = 'Data_phb_'//trim(acp)//'_'//trim(dnr)//'.'//trim(argv)
  open(15,file=infile,status='old')
  do j = 1, nline
    read(15,*) t(j), read_c
    c(j) = c(j) + read_c
  end do
  close(15)
end do


outfile = 'Data_phb_'//trim(acp)//'_'//trim(dnr)//'.xvg'
open(17, file=outfile)
do i = 1, nline
  write(17,'(2f20.7)') t(i) , c(i)/c(1)
end do
close(17)










end program bbb
