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

program pz2
use info
use calc
use xdr, only: xtcfile
implicit none
type type_v
  integer*1, allocatable :: flag_z(:)
end type type_v

type(xtcfile) :: xtc

!======================================================================!
!============================== variable ==============================!
!======================================================================!
! parameter
integer, parameter :: ON = 1, UP = 1
integer, parameter :: OFF = 0, DOWN = 0
real*8, parameter :: dr = 0.002d0 !nm
real*8, parameter :: bin = 0.01d0 
real*8, parameter :: n_max = 3.0d0 !nm
real*8, parameter ::pi = dacos(-1d0)
! MDinfo
integer :: nstep
real*8 :: dt
character(50) :: trjfile
! MOLinfo
type(MOL_atom),allocatable :: mol(:)
character(5),allocatable :: order(:)
! GRPinfo
type(group_atom) :: grp
integer, allocatable :: list(:)
! get_boxdim
real :: box_dim(3)
type(type_v) :: types(2)
! allocatable
real*8, allocatable :: pos(:,:)
real*8, allocatable :: g(:,:)
integer, allocatable ::leaf(:)
real*8, allocatable ::z_list(:)
! others
integer :: i, j, k, l
integer :: nhist, ig, flag, id
integer :: step, num, npos
integer :: ini_pos(2), posatm
real*8 :: v(2), r, r_max, dv, rho
real*8 :: posg, newpos(3)
real*8 :: ave_z(0:1),s2(0:1)
integer :: n_leaf(0:1)
character(50) :: outfile
character(50) :: argv
!======================================================================!
print *, '*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*'
print*, ''
print *, '                      == p =='
print *, ''

!MDinfo
call  mdinfo(trjfile,nstep,dt,1)
call  molinfo(mol, order)

!MOLinfo
num = COMMAND_ARGUMENT_COUNT()
if ( num == 0 ) call warn_arg()
call getarg(1,argv)
grp%gname = trim(argv)
call grpinfo(order, mol, grp)

do i = 1, size(mol)
  if ( mol(i) % flag_com == 1 ) then
    call bondinfo(mol(i))
  end if
end do

call xtc % init(trjfile)
call xtc % read
call get_boxdim(xtc%box, box_dim)
r_max = box_dim(1)/2d0
nhist = ceiling(n_max/bin)
allocate(pos(3,xtc % natoms))
allocate(g(0:nhist,0:1))
g(:,:) = 0d0
npos = size(grp % index_list)
allocate(list(size(grp % index_list)))
allocate(leaf(npos))
allocate(z_list(npos))
list = grp % index_list

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
  ave_z = 0d0
  n_leaf(:) = 0
  s2 = 0d0

  do i = 1, npos
    id = list(i)
    call z_change(pos(1:3,id), box_dim(3), posg, newpos)
    if ( newpos(3) < 0d0 ) then 
      leaf(i) = DOWN
      n_leaf(DOWN) = n_leaf(DOWN) + 1
    else
      leaf(i) = UP
      n_leaf(UP) = n_leaf(UP) + 1
    end if
    ave_z(leaf(i)) = ave_z(leaf(i)) + newpos(3)
    z_list(i) = newpos(3)
  end do
  ave_z(:) = ave_z(:) / dble(n_leaf(:)) ! average
  do i = 1, npos
    ig = ceiling(abs(z_list(i)-ave_z(leaf(i)))/bin)
    if (ig ==0) print*, z_list(i)-ave_z(leaf(i))
    g(ig,leaf(i)) = g(ig,leaf(i)) + 1d0
  end do

  call xtc % read
  call get_boxdim(xtc%box, box_dim)
end do

call xtc % close

print*, 'outputting...'
outfile='Data_pz2.xvg'

open(15,file=outfile)
do i = 1, nhist
  !write(15,'(f8.3,f20.6)') (i -0.5)*dr, g(i)
  write(15,'(f15.8,f20.15)') dble(i)*bin, sum(g(i,:))/(bin*dble(npos*2*(nstep+1)))
!  if (i ==274)   write(*,'(f15.8,f10.5,f20.15)') dble(i)*bin,sum(g(i,:)), sum(g(i,:))/dble(2*(nstep+1))
!  if (i ==275)   write(*,'(f15.8,f10.5,f20.15)') dble(i)*bin,sum(g(i,:)), sum(g(i,:))/dble(2*(nstep+1))
!  if (i ==276)   write(*,'(f15.8,f10.5,f20.15)') dble(i)*bin,sum(g(i,:)), sum(g(i,:))/dble(2*(nstep+1))
  !write(15,'(f15.8,3f20.15)') dble(i)*bin, sum(g(i,:))/dble(2*(nstep+1)),g(i,:)/dble(nstep+1)
end do

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

subroutine warn_arg()
  implicit none
  print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  print *, '|                 !! ERROR !!                  |'
  print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  print *, ''
  print *, '! This program needs 2 arguments !'
  print *, ''
  print *, './pz2 a b'
  print *, ''
  print *, 'a : atom1 : center (e.g. phosphate)'
  print *, 'b : atom2 : distribution (e.g. ow)'
  print *, ''
  print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  print *, '|               !! YARINAOSE !!                |'
  print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  stop

end subroutine warn_arg



end program pz2
