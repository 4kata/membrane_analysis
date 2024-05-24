program tcf_for_stay_region
use info
use calc
use xdr, only: xtcfile
implicit none

type(xtcfile) :: xtc

!======================================================================!
!============================== variable ==============================!
!======================================================================!
! parameter
integer, parameter :: ON = 1, UP = 1
integer, parameter :: OFF = 0, DOWN = 0
integer, parameter :: n_region = 3
real*8, parameter :: logn = 1.1d0
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
integer, allocatable :: list_c(:), list_d(:)
! get_boxdim
real :: box_dim(3)
! allocatable
real*8, allocatable :: pos(:,:)
real*8, allocatable :: cutoff_length(:)
real*8, allocatable :: c0(:),ct_up(:), ct_down(:),ct_stay(:)
integer*1, allocatable :: state(:,:), flag(:,:)
real*8, allocatable :: sum_state(:)
real*8, allocatable :: c_up(:,:), c_down(:,:), c_stay(:,:)
integer, allocatable ::id_leaf(:,:)
integer ::leaf
! others
integer :: i, j, k, l
integer :: nhist, ig, id, id_min
integer :: step, num, npos_c, npos_d
integer :: logstep, interval, runstep
integer :: i_state, e_state
real*8 :: z, r_max, dv, v(2), t
real*8 :: posg, newpos(3)
real*8 :: zg(0:1), sum_box(2)
integer :: n_leaf(0:1)
real*8 :: num_div
character(50) :: outfile
character(50) :: argv, argw
character(100) :: fmt
!======================================================================!
print *, '*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*'
print*, ''
print *, '                   == zp 0 = p =='
print *, ''

!MDinfo
call  mdinfo(trjfile,nstep,dt,1)
call  molinfo(mol, order)
if ( logn /= 1d0) then
  logstep = int(log(real(nstep))/log(logn))
else
  logstep = nstep
end if


!MOLinfo
num = COMMAND_ARGUMENT_COUNT()
if ( num == 0 ) call warn_arg()
do i = 1, 2
  call getarg(i,argv)
  grp(i)%gname = trim(argv)
  call grpinfo(order, mol, grp(i))
end do
allocate(cutoff_length(n_region - 1))
outfile='Data_tcf_transition_'//trim(argv)//'.xvg'
do i = 3 , 2 + n_region - 1
  call getarg(i,argv)
  read(argv,*) cutoff_length(i-2)
  write(*,'(a6,i2,a3,i2,a7,a3)') 'region', i-2, ' - ', i-1,trim(argv), ' nm'
end do

!Bondinfo
do i = 1, size(mol)
  if ( mol(i) % flag_com == 1 ) then
    call bondinfo(mol(i))
  end if
end do

!Open xtc
call xtc % init(trjfile)
call xtc % read
call get_boxdim(xtc%box, box_dim)

r_max = box_dim(3)/2d0
allocate(pos(3,xtc % natoms))
sum_box(:)=0
npos_c = size(grp(1) % index_list)
allocate(id_leaf(npos_c,0:1))
allocate(list_c(size(grp(1) % index_list)))
list_c = grp(1) % index_list
npos_d = size(grp(2) % index_list)
allocate(list_d(size(grp(2) % index_list)))
list_d = grp(2) % index_list

allocate(state(npos_d,nstep+1))
allocate(c0(n_region))
allocate(ct_up(1:n_region-1))
allocate(ct_down(2:n_region))
allocate(ct_stay(1:n_region))
allocate(flag(npos_d,nstep+1))
allocate(sum_state(n_region))
allocate(c_up(1:n_region-1,0:logstep))
allocate(c_down(2:n_region,0:logstep))
allocate(c_stay(1:n_region,0:logstep))

sum_state(1:n_region) = 0d0
!sum_state(n_region) = 1d0 !(nstep + 1)*npos_d
num_div = 1d0/dble(nstep+1) 
print*, num_div
print*, sum_state
state(:,:) = n_region

!Read trj
num = nstep
if(nstep == 0) num = 1
print*, 'reading trajctory and calculating length...'
do step = 1,nstep + 1
  if (mod(step*10,num) == 0) print *, int(dble(step*100)/dble(num)),'%'
  pos(:,:) = xtc % pos(:,:)
  call centering(mol, pos, posg,box_dim(3))
  if ( posg /= posg ) then
    posg=0
  end if
  n_leaf(:) = 0
  id_leaf(:,:) = 0

  do i = 1, npos_c
    id = list_c(i)
    call z_change(pos(1:3,id), box_dim(3), posg, newpos)
    call check_leaf(newpos(3), leaf)
    pos(3,id) = newpos(3)
    n_leaf(leaf) = n_leaf(leaf) + 1
    id_leaf(n_leaf(leaf),leaf) = id
  end do
!$omp parallel shared(pos,npos_d,npos_c,box_dim,nhist,list_d,id_leaf,posg,state,cutoff_length,num_div) &
!$omp & private(j,id,id_min,v,z,ig,newpos,leaf) &
!$omp & reduction(+:sum_state)
!$omp do  
  do i = 1, npos_d
    id = list_d(i)
    call z_change(pos(1:3,id), box_dim(3), posg, newpos)
    call check_leaf(newpos(3), leaf)
    z = 9999999
    do j = 1, n_leaf(leaf)
      v(:) = pos(1:2,id) - pos(1:2,id_leaf(j,leaf))
      if ( z > sqrt(sum(v**2)) ) then
        z = sqrt(sum(v**2))
        id_min = id_leaf(j,leaf)
      end if
    end do
    z = abs(newpos(3)) - abs(pos(3,id_min))
    do j = 1, n_region - 1
      if  ( z < cutoff_length(j) ) then
        state(i,step) = j
        sum_state(j) = sum_state(j) + 1d0*num_div !dble(((nstep + 1)*npos_d))
        sum_state(n_region) = sum_state(n_region) + 1d0*num_div !dble((nstep + 1)*npos_d)
        !if ( sum_state(3) < 0) then
        !print*, num_div, nstep+1
        !stop

        !end if
        exit
      end if
    end do
  end do
!$omp end do
!$omp end parallel

  sum_box(1:2) = sum_box(1:2) + box_dim(1:2)
  call xtc % read
  call get_boxdim(xtc%box, box_dim)
end do
sum_box(1:2) = sum_box(1:2)/dble(nstep + 1)
print*, sum_state
sum_state(n_region) = dble(npos_d) - sum_state(n_region)
print*, 'total : ', sum(sum_state)

call xtc % close

!==========================================================
print*, 'tcf part...'
c_up(:,:) = 0d0
c_down(:,:) = 0d0
c_stay(:,:) = 0d0
flag(:,:) = -1

do step = 0, logstep
!  if (mod(step*10,num) == 0) print *, int(dble(step*100)/dble(num)),'%'
  interval = calc_interval(logn, step)
  runstep = nstep + 1 - interval
  num_div = 1d0/dble(runstep)
!$omp parallel shared(npos_d,runstep,interval,state,flag,num_div) &
!$omp & private(j,i_state,e_state) &
!$omp & reduction(+:c_up,c_down,c_stay)
!$omp do  
  do i = 1, runstep
    do j = 1, npos_d
      i_state = state(j,i)
      if ( flag(j,i) ==  -1 ) then
        e_state = state(j,i + interval)
        if ( i_state > e_state ) then
          c_down(i_state,step) = c_down(i_state,step) + 1d0*num_div
          flag(j,i) = OFF
        else if (i_state < e_state ) then
          c_up(i_state,step) = c_up(i_state,step) + 1d0*num_div
          flag(j,i) = ON
        !  !if ( i_state == 1) print*, c_up(i_state,logstep)
        else if (i_state == e_state) then
          c_stay(i_state,step) = c_stay(i_state,step) + 1d0*num_div
        end if
      else if ( flag(j,i) == ON ) then
        c_up(i_state,step) = c_up(i_state,step) + 1d0*num_div
      else if ( flag(j,i) == OFF ) then
        c_down(i_state,step) = c_down(i_state,step) + 1d0*num_div
      end if
    end do
  end do
!$omp end do
!$omp end parallel
end do
!stop

!==========================================================
!!!!!thinking
print*, 'outputting...'
open(15,file=outfile)
fmt='t '
do i = 1, n_region
  write(argv,'(I2.2)') i
  write(argw,'(I2.2)') i
  fmt = trim(fmt)//' '//trim(argv)//'=>'//trim(argw)//' '
end do
do i = 1, n_region-1
  write(argv,'(I2.2)') i
  write(argw,'(I2.2)') i+1
  fmt = trim(fmt)//' '//trim(argv)//'=>'//trim(argw)//' '
end do
do i = 1, n_region-1
  write(argv,'(I2.2)') i
  write(argw,'(I2.2)') i+1
  fmt = trim(fmt)//' '//trim(argw)//'=>'//trim(argv)//' '
end do
write(15,'(A)') trim(fmt)
write(fmt, '(a, i0, a)') '(f15.5,',2*(n_region-1)+n_region, 'f13.8)'
c0(:) = dble(sum_state(:))!/dble((nstep+1)*npos_d)
do step = 0, logstep
  interval = calc_interval(logn, step)
  t = interval * dt
  runstep = nstep + 1 - interval
  ct_up(:) = c_up(:,step)
  ct_down(:) = c_down(:,step)
  ct_stay(:) = c_stay(:,step)
  ct_up(1:n_region-1) = ct_up(1:n_region-1)/c0(1:n_region-1)
  ct_down(2:n_region) = ct_down(2:n_region)/c0(2:n_region)
  ct_stay(:) = ct_stay(:)/c0(:)
  write(15,fmt) t, ct_stay(:),ct_up(:), ct_down(:)
    
end do
  write(*,fmt) t, ct_stay(:),ct_up(:), ct_down(:)
print *,''
print *, '             == MISSION COMPLETE !!! =='
print *, ''
print *, '*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*'

contains

subroutine z_change(pos, box_dim, posg, newpos)
  implicit none
  real, intent(in) :: box_dim
  real*8, intent(in) :: pos(3), posg
  real*8, intent(out) :: newpos(3)
  real*8 :: z

  newpos(:)= pos(:)

  if ( posg >= box_dim/2d0 ) then
    z = pos(3) !- box_dim(3)*nint(pos(3, posatm)/box_dim(3))
    if ( abs(z) > box_dim) z = z - box_dim*nint(z/box_dim)
    z = z - posg
    if ( abs(z) > box_dim/2d0 ) z = z + box_dim
  else
    z = pos(3) !- box_dim(3)*nint(pos(3, posatm)/box_dim(3))
    if ( abs(z) > box_dim) z = z - box_dim*nint(z/box_dim)
    z = z - posg
    if ( abs(z) > box_dim/2d0 ) z = z - box_dim
  end if

  !z = abs(z)
  newpos(3) = z
endsubroutine z_change


subroutine check_leaf(pos, leaf)
  implicit none
  real*8, intent(in) :: pos
  integer,intent(out) :: leaf
  integer, parameter :: DOWN=0, UP=1

  if ( pos < 0d0 ) then 
    leaf = DOWN
  else
    leaf = UP
  end if
end subroutine check_leaf

integer function calc_interval(logn, step)
  implicit none
  real*8, intent(in) :: logn
  integer, intent(in) :: step
  
  calc_interval = step
  if ( logn /= 1d0 ) then
    calc_interval = int(logn**dble(step))
  end if
end function calc_interval


subroutine warn_arg()
  implicit none
  print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  print *, '|                 !! ERROR !!                  |'
  print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  print *, ''
  print *, '! This program needs 4 arguments !'
  print *, ''
  print *, './tcf_transiton'
  print *, ''
  print *, 'a : atom1 = centering atom'
  print *, 'b : atom2 = distribution group'
  print *, 'c : z coordinate which splits region between 1 and 2'
  print *, 'd : z coordinate which splits region between 2 and 3'
  print *, ''
  print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  print *, '|               !! YARINAOSE !!                |'
  print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  stop

end subroutine warn_arg



end program tcf_for_stay_region
