program projectblue 
!! z_profile center = Zg (avarage position of p atom in each snapshot)
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
real*8, parameter :: bin = 0.1d0 
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
type(group_atom), allocatable :: grp_a(:)
integer, allocatable :: list_c(:), list_d(:), list_a(:)
! get_boxdim
real :: box_dim(3)
! allocatable
real*8, allocatable :: pos(:,:)
real*8, allocatable :: g(:)
integer, allocatable ::id_leaf(:,:)
integer*8, allocatable :: hb_count(:)
integer ::leaf
! others
integer :: i, j, k, l
integer :: nhist, ig, flag, id, id_min
integer :: id_1, id_2
integer :: step, num, npos_c, npos_d
real*8 :: z, r_max, dv, v(2)
real*8 :: posg, newpos(3)
real*8 :: zg(0:1), sum_box(2)
integer :: n_leaf(0:1)
integer :: id_h(2), ncite
character(50) :: outfile
character(50) :: argv
!======================================================================!
print *, '*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*'
print*, ''
print *, '                   == zp 0 = p =='
print *, ''

!MDinfo
call  mdinfo(trjfile,nstep,dt,1)
!MOLinfo
call  molinfo(mol, order)

!==================== get argumaents =====================
!Getargs
num = COMMAND_ARGUMENT_COUNT()
allocate(grp_a(num-2))
if ( num == 0 ) call warn_arg()
j=0
do i = 1, 2
  call getarg(i,argv)
  grp(i)%gname = trim(argv)
  call grpinfo(order, mol, grp(i))
end do
outfile='Data_zphb_near_'//trim(argv)//'.xvg'
if ( num >= 3 ) then
  do i = 3, num
    call getarg(i,argv)
    grp_a(i-2)%gname = argv
    call grpinfo(order, mol, grp_a(i-2))
    j = j + size(grp_a(i-2)%index_list)
  end do
  allocate(list_a(j))
  j = 0
  do i = 1, num - 2
    list_a(j+1:j+size(grp_a(i)%index_list)) = grp_a(i)%index_list(:)
    j = j + size(grp_a(i)%index_list)
  end do
end if

!Bondinfo
do i = 1, size(mol)
  if ( mol(i) % flag_com == 1 ) then
    call bondinfo(mol(i))
  end if
end do
print*, grp(2)%gname
call define_ncite(grp(2)%gname,mol,grp(2),ncite,id_h)

!==================== info from xtc =====================
call xtc % init(trjfile)
call xtc % read
call get_boxdim(xtc%box, box_dim)

r_max = box_dim(3)/2d0
nhist = ceiling(r_max/bin)
allocate(pos(3,xtc % natoms))
allocate(g(-nhist-100:nhist+100))
allocate(hb_count(-nhist-100:nhist+100))
g(:) = 0d0
hb_count(:) = 0
sum_box(:)=0
npos_c = size(grp(1) % index_list)
allocate(id_leaf(npos_c,0:1))
allocate(list_c(size(grp(1) % index_list)))
list_c = grp(1) % index_list
npos_d = size(grp(2) % index_list)
allocate(list_d(size(grp(2) % index_list)))
list_d = grp(2) % index_list

!Read trj
num = nstep
if(nstep == 0) num = 1
!========================================================
print*, 'reading trajctory and calculating length...'
do step = 1,nstep + 1
  if (mod(step*10,num) == 0) print *, int(dble(step*100)/dble(num)),'%'
  pos(:,:) = xtc % pos(:,:)
  call centering(mol, pos, posg,box_dim(3))
  if ( posg /= posg ) then
    posg=0
  end if
  zg = 0d0
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
!$omp parallel shared(pos,npos_d,box_dim,nhist,list_d,list_a,id_leaf,id_h,posg,n_leaf,ncite) &
!$omp & private(j,k,id_1,id_2,ig,newpos,leaf) &
!$omp & reduction(+:g, hb_count)
!$omp do  
  do i = 1, npos_d
    id_1 = list_d(i)
    call z_change(pos(1:3,id_1), box_dim(3), posg, newpos)
    call check_leaf(newpos(3), leaf)
    ig = get_bin(bin,pos,newpos(3),id_1,id_leaf,n_leaf(leaf),leaf)

    !if( dble(ig)*bin < box_dim(3)/4d0 + 0.35d0) then
      g(ig) = g(ig) + 1d0
  
      ! not water oxygen
      do j = 1, size(list_a)
        id_2 = list_a(j)
        k = check_hbond(pos,id_1,id_2,id_h,ncite,box_dim)
        hb_count(ig) = hb_count(ig) + k
      end do
      ! water
      do j = 1, npos_d
        id_2 = list_d(j)
        k = check_hbond(pos,id_1,id_2,id_h,ncite,box_dim)
        if ( k == 1 ) then
          hb_count(ig) = hb_count(ig) + 1
          call z_change(pos(1:3,id_2), box_dim(3), posg, newpos)
          call check_leaf(newpos(3),leaf)
          ig = get_bin(bin,pos,newpos(3),id_2,id_leaf,n_leaf(leaf),leaf)
          hb_count(ig) = hb_count(ig) + 1
        end if
      end do
  end do
!$omp end do
!$omp end parallel
  sum_box(1:2) = sum_box(1:2) + box_dim(1:2)
  call xtc % read
  call get_boxdim(xtc%box, box_dim)
end do
  print*, step,  hb_count(-nhist-100) ,ig
!stop
sum_box(1:2) = sum_box(1:2)/dble(nstep + 1)

call xtc % close
!========================================================
print*, 'outputting...'

print*, hb_count(-nhist-100) 
  print*, step,  hb_count(-nhist-100) ,ig
dv = sum_box(1)*sum_box(2)*bin
open(15,file=outfile)
do i = -ceiling(box_dim(3)/4d0/bin), ceiling(box_dim(3)/4d0/bin)
  if ( g(i) /= 0 ) then
    write(15,'(f15.8,f20.15)') dble(i)*bin, dble(hb_count(i))/(dble(g(i)))
  end if
end do
print*, dble(sum(hb_count))/dble(sum(g))

close(15)

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
  z = pos(3) !- box_dim(3)*nint(pos(3, posatm)/box_dim(3))
  if ( z < 0 ) z = z + box_dim

  if ( posg >= box_dim/2d0 ) then
    if ( abs(z) > box_dim) z = z - box_dim*nint(z/box_dim)
    z = z - posg
    if ( abs(z) > box_dim/2d0 ) z = z + box_dim
  else
    if ( abs(z) > box_dim) then
      z = z - box_dim*nint(z/box_dim)
    end if
    z = z - posg
    if ( abs(z) > box_dim/2d0 ) then 
      z = z - box_dim
    end if
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

integer function check_hbond(pos,id_1,id_2,id_h,ncite,box_dim)
  implicit none
  integer, parameter :: ON = 1, OFF = 0
  real*8, allocatable,intent(in) :: pos(:,:)
  integer, intent(in) :: id_1, id_2, id_h(2)
  integer, intent(in) :: ncite
  real, intent(in) :: box_dim(3)
  integer :: i_cite
  real*8 :: v_oo(3), v_oh(3)
  real*8 :: r_oo, r_oh
  real*8 ::beta

  check_hbond = 0
  v_oo(:) = pos(1:3,id_2) - pos(1:3,id_1)
  call length(box_dim, v_oo(1:3), r_oo)
  if ( r_oo > 0d0 .and. r_oo < 0.35d0) then
    do i_cite = 1, ncite
      v_oh = pos(1:3,id_1 + id_h(i_cite)) - pos(1:3,id_1)
      call length(box_dim, v_oh(1:3), r_oh)
      beta = calc_degree(v_oo, v_oh, r_oo, r_oh, ON)
      if ( beta < 30d0) then
        check_hbond = 1
      end if
    end do
  end if
end function check_hbond

integer function get_bin(bin,pos,newpos,id_1,id_leaf,n_leaf,leaf)
  implicit none
  real*8, intent(in) :: bin, newpos
  integer, intent(in) :: id_1, leaf, n_leaf
  real*8, allocatable, intent(in) :: pos(:,:)
  integer, allocatable, intent(in) :: id_leaf(:,:)
  real*8 :: z, v(2)
  integer :: j, id_min

  z = 9999999
  do j = 1, n_leaf
    v(:) = pos(1:2,id_1) - pos(1:2,id_leaf(j,leaf))
    if ( z > sqrt(sum(v**2)) ) then
      z = sqrt(sum(v**2))
      id_min = id_leaf(j,leaf)
    end if
  end do
  z = abs(newpos) - abs(pos(3,id_min))
  !write(*,'(2i10,3f10.5)') id_1, id_min, abs(newpos) , abs(pos(3,id_min)), z
  get_bin = ceiling(z/bin)
end function get_bin

subroutine define_ncite(argv,mol,grp,ncite,id_h)
  implicit none
  character(50), intent(in) :: argv
  type(MOL_atom), allocatable, intent(inout) :: mol(:)
  type(group_atom), intent(inout) :: grp
  integer, intent(out) :: ncite
  integer, intent(out) :: id_h(2)
  integer :: i, j, num
  integer, parameter :: OFF = 0, ON = 1

  if ( argv == 'ow') then
    call water_info(mol(grp%order))
    ncite = 2
    id_h(1) = 1
    id_h(2) = 2
  else
    !call bondinfo(mol(grp%order))
    do i = 1, grp%natoms
      do j = 1, mol(grp%order)%nbonds(grp%list(i))
        num = mol(grp%order)%list_bonds(j,grp%list(i))
        if (mol(grp%order)%flag_heavy(num) == OFF) then
          ncite = ncite + 1
          id_h(i) = num - grp%list(i)
        end if
      end do
    end do
  end if
end subroutine define_ncite

subroutine warn_arg()
  implicit none
  print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  print *, '|                 !! ERROR !!                  |'
  print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  print *, ''
  print *, '! This program needs more than 2 arguments !'
  print *, ''
  print *, './zphb_near a b ....'
  print *, ''
  print *, 'a    : atom1 = centering atom'
  print *, 'b    : atom2 = h-bond donor'
  print *, 'c... : atomN = h-bond acceptor'
  print *, ''
  print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  print *, '|                !! DAMEDESU !!                |'
  print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  stop

end subroutine warn_arg



end program projectblue
