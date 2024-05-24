program phb_tcf
use info
use calc
use xdr, only: xtcfile
implicit none
type(xtcfile) :: xtc

!======================================================================!
!============================== variable ==============================!
!======================================================================!
! parameter
integer, parameter :: ON = 1
integer, parameter :: OFF = 0
real*8, parameter :: logn = 1.3d0 
real*8, parameter ::pi = dacos(-1d0)
! MDinfo
integer :: nstep
real*8 :: dt
character(50) :: trjfile
! MOLinfo
type(MOL_atom),allocatable :: mol(:)
character(5),allocatable :: order(:)
! GRPinfo
type(group_atom) :: grp(2)
integer, allocatable :: list1(:), list2(:)
! get_boxdim
real :: box_dim(3)
! allocatable
real, allocatable :: pos(:,:)
integer, allocatable :: n_count(:,:)
integer, allocatable :: h(:,:,:,:)
integer*8, allocatable :: ht(:)
! others
integer :: i, j, k, l
integer :: i_cite, i_tag, j_tag
integer :: logstep, kstep, interval
integer :: step, num, npos1, npos2, id_1, id_2, id_h(2)
integer :: ini_pos(2), posatm, ncite
integer :: chain_id1, chain_id2, flag
integer :: n_ini, n_end
real*8 :: v_oo(3), r_oo, v_oh(3), r_oh, beta
real*8 :: phb0, phb, t, a
character(50) :: outfile,arg1, arg2, argv
!======================================================================!
print *, '*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*'
print*, ''
print *, '                      == pmf =='
print *, ''

call  mdinfo(trjfile,nstep,dt,1)
call  molinfo(mol, order)
num = COMMAND_ARGUMENT_COUNT()
logstep = int(log(real(nstep))/log(logn))

if ( num == 0 ) call warn_arg()
call getarg(1,arg1)
call getarg(2,arg2)
print*, 'donor side    : ', trim(arg1)
print*, 'accepcor side : ', trim(arg2)
print *, ''


call getarg(3,argv)
read(argv,*) i
write(argv,'(I4.4)') i
outfile='Data_phb_'//trim(arg1)//'_'//trim(arg2)//'.'//trim(argv)
call getarg(4,argv)
read(argv,*) n_ini
call getarg(5,argv)
read(argv,*) n_end


grp%gname = trim(arg1)
call grpinfo(order, mol, grp(1))
grp%gname = trim(arg2)
call grpinfo(order, mol, grp(2))
allocate(list1(size(grp(1)%index_list)))
allocate(list2(size(grp(2)%index_list)))
list1 = grp(1)%index_list
list2 = grp(2)%index_list
call inipos(mol,mol(grp(1)%order),ini_pos(1))
call inipos(mol,mol(grp(2)%order),ini_pos(2))
call xtc % init(trjfile)
call xtc % read
call get_boxdim(xtc%box, box_dim)
allocate(pos(3,xtc % natoms))
allocate(ht(0:logstep))
ht(:) = 0


npos1 = mol(grp(1)%order)%nmol*grp(1)%natoms
npos2 = mol(grp(2)%order)%nmol*grp(2)%natoms

ncite = 0
if ( arg1 == 'ow') then
  call water_info(mol(grp(1)%order))
  ncite = 2
  id_h(1) = 1
  id_h(2) = 2
else
  call bondinfo(mol(grp(1)%order))
  do i = 1, grp(1)%natoms
    do j = 1, mol(grp(1)%order)%nbonds(grp(1)%list(i))
      num = mol(grp(1)%order)%list_bonds(j,grp(1)%list(i))
      if (mol(grp(1)%order)%flag_heavy(num) == OFF) then
        ncite = ncite + 1
        id_h(i) = num - grp(1)%list(i)
      end if
    end do
  end do
end if

allocate(h(10,n_ini:n_end,0:nstep,ncite))
allocate(n_count(n_ini:n_end,ncite))
h=0


num = nstep
a=1/dble(num)
if(nstep == 0) num = 1
print*, 'reading trajctory and calculating length...'
do step = 0, nstep
!print*, step
  if (mod(step*10,num) == 0) print *, int(dble(step*100)*a),'%'
  pos(:,:) = xtc % pos(:,:)
  n_count(:,:) = 0
!$omp parallel shared(step, pos, n_ini, n_end, npos2,ncite, box_dim, list1,list2,id_h) &
!$omp & private(j, i_cite, id_1, id_2, v_oo, r_oo,v_oh,r_oh,beta) &
!$omp & reduction(+:n_count)
!$omp do collapse(2)
  do i = n_ini, n_end ! donor side
    do j = 1,  npos2 ! acceptor side
    id_1 = list1(i)
      id_2 = list2(j)
      v_oo(:) = pos(:,id_2) - pos(:,id_1)
      call length(box_dim, v_oo, r_oo)
      if ( r_oo > 0.0d0 .and. r_oo <= 0.35d0) then
        do i_cite = 1, ncite
          v_oh = pos(:,id_1 + id_h(i_cite)) - pos(:,id_1)
          call length(box_dim, v_oh, r_oh)
          beta = calc_degree(v_oo,v_oh,r_oo,r_oh,ON)
          if ( beta < 30.0d0 ) then
            n_count(i,i_cite) = n_count(i,i_cite) + 1
         !   print*, n_count(i,i_cite), i, step, i_cite
            h(n_count(i,i_cite), i, step, i_cite) = id_2
          end if
        end do
      end if
    end do
  end do
!$omp end do
!$omp end parallel
  call xtc % read
  call get_boxdim(xtc%box, box_dim)
end do

call xtc % close
print *, 'END Checking Hbond !!!'
! making TCF.....
print *, 'Start Counting H-Bond Time ...'

! t = 0
print*, ncite, npos1, npos2
phb0 = 0d0
!$omp parallel shared(nstep, n_ini, n_end, h) &
!$omp & private(j, i_cite, i_tag) &
!$omp & reduction(+: phb0)
!$omp do
do step = 0, nstep
  do i = n_ini, n_end
    do i_cite = 1, ncite
      do i_tag = 1, 10
        j = h(i_tag,i,step,i_cite)
        if (j /= 0 ) then
          phb0 = phb0 + 1d0
        else
          exit
        end if
      end do
    end do
  end do
end do
!$omp end do
!$omp end parallel
print*, phb0
! t = t
do step = 0, logstep
  interval = int(logn**step)
  kstep = nstep - interval
!$omp parallel shared(step, n_ini, n_end, h, interval, kstep) &
!$omp & private(j, k, i_tag, j_tag, i_cite) &
!$omp & reduction(+: ht)
!$omp do
  do k = 0, kstep
    do i = n_ini, n_end
      do i_cite = 1, ncite
        do i_tag = 1, 10
          j = h(i_tag,i,k,i_cite)
          if (j /= 0 ) then
            do j_tag = 1, 10
              if (h(j_tag,i,k+interval,i_cite) == j ) then
                ht(step) = ht(step) + 1
                exit
              end if
            end do
          else
            exit
          end if
        end do
      end do
    end do
  end do
!$omp end do
!$omp end parallel
end do


print*, 'outputting...'
open(15,file=outfile)



t=0d0
print*, phb0
phb0 = phb0/(dble(nstep + 1))

write(15,'(2f20.7)') t, phb0

do i = 0, logstep
  interval = int(logn**i)
  t = interval * dt
  phb = dble(ht(i))/(dble(nstep - interval))
  write(15,'(2f20.7)') t, phb
end do


close(15)

print *,''
print *, '             == MISSION COMPLETE !!! =='
print *, ''
print *, '*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*'

contains


subroutine warn_arg()
  implicit none
  print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  print *, '|                 !! ERROR !!                  |'
  print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  print *, ''
  print *, '! This program needs 2 arguments !'
  print *, ''
  print *, './phb a b'
  print *, ''
  print *, 'a : atom1 (donar)'
  print *, 'b : atom2 (acceptor)'
  print *, ''
  print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  print *, '|               !! YARINAOSHI !!                |'
  print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  stop

end subroutine warn_arg



end program phb_tcf
