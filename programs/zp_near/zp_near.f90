! Copyright (C) 2024 Kokoro Shikata
!
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

program polaris 
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
real*8, parameter :: bin = 0.025d0 
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
real*8, allocatable :: g(:)
integer, allocatable ::id_leaf(:,:)
integer ::leaf
! others
integer :: i, j, k, l
integer :: nhist, ig, flag, id, id_min
integer :: step, num, npos_c, npos_d
real*8 :: z, r_max, dv, v(2)
real*8 :: posg, newpos(3)
real*8 :: zg(0:1), sum_box(2)
integer :: n_leaf(0:1)
character(50) :: outfile
character(50) :: argv
!======================================================================!
print *, '*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*'
print*, ''
print *, '                   == zp 0 = p =='
print *, ''

!MDinfo
call  mdinfo(trjfile,nstep,dt,1)
call  molinfo(mol, order)

!MOLinfo
num = COMMAND_ARGUMENT_COUNT()
if ( num == 0 ) call warn_arg()
do i = 1, 2
  call getarg(i,argv)
  grp(i)%gname = trim(argv)
  call grpinfo(order, mol, grp(i))
end do
outfile='Data_zp_near_'//trim(argv)//'.xvg'

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
nhist = ceiling(r_max/bin)
allocate(pos(3,xtc % natoms))
allocate(g(-nhist:nhist))
g(:) = 0d0
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
!    zg(leaf) = zg(leaf) + newpos(3)
  end do
!  zg(:) = zg(:) / dble(n_leaf(:)) ! average
!$omp parallel shared(pos,npos_d,npos_c,box_dim,nhist,list_d,id_leaf,posg) &
!$omp & private(j,id,id_min,v,z,ig,newpos,leaf) &
!$omp & reduction(+:g)
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
    ig = ceiling(z/bin)
!    write(*,'(3f10.5,i5)') newpos(3), pos(3,id_min) ,newpos(3)- pos(3,id_min), ig
!    print* , pos(:,id)
!    print* , pos(:,id_min)
    
    if ( ig < nhist ) then
      g(ig) = g(ig) + 1d0
    end if
  end do
!$omp end do
!$omp end parallel

  sum_box(1:2) = sum_box(1:2) + box_dim(1:2)
  call xtc % read
  call get_boxdim(xtc%box, box_dim)
end do
sum_box(1:2) = sum_box(1:2)/dble(nstep + 1)

call xtc % close

print*, 'outputting...'

dv = sum_box(1)*sum_box(2)*bin
open(15,file=outfile)
do i = -nhist , nhist - 10
!  if ( g(i) /= 0 ) then
    write(15,'(f15.8,f20.15)') dble(i)*bin, g(i)/(2d0*dv*dble((nstep+1)))
!  end if
end do
print *, sum(g(:)/(dble(npos_d*(nstep+1))))
print *, sum(g)

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

subroutine warn_arg()
  implicit none
  print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  print *, '|                 !! ERROR !!                  |'
  print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  print *, ''
  print *, '! This program needs 1 arguments !'
  print *, ''
  print *, './zp_near'
  print *, ''
  print *, 'a : atom1 = centering atom'
  print *, 'b : atom2 = distribution group'
  print *, ''
  print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  print *, '|               !! YARINAOSE !!                |'
  print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  stop

end subroutine warn_arg



end program polaris
