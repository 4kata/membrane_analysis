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

program pmf
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
real*8, parameter :: dr = 0.002d0 !nm
real*8, parameter :: db = 1d0 !nm
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
real*8, allocatable :: g(:,:)
! others
integer :: i, j, k, l
integer :: nrhist, ig
integer :: nbhist, ib
integer :: step, num, npos1, npos2, id_1(0:1), id_2(0:1), id_h(2)
integer :: ini_pos(2), posatm, ncite
integer :: chain_id1, chain_id2, flag
real*8 :: v_oo(3), r_oo, v_oh(3), r_oh, beta
real*8 :: r_max, dv, rho, maxg
character(50) :: outfile,arg1, arg2
!======================================================================!
print *, '*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*'
print*, ''
print *, '                      == pmf =='
print *, ''

call  mdinfo(trjfile,nstep,dt,1)
call  molinfo(mol, order)
!num = iargc()
num=COMMAND_ARGUMENT_COUNT()
if ( num == 0 ) call warn_arg()
call getarg(1,arg1)
call getarg(2,arg2)
print*, 'donor side    : ', trim(arg1)
print*, 'accepcor side : ', trim(arg2)
print *, ''
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
r_max = box_dim(1)/2d0
nrhist = ceiling(r_max/dr)
nbhist = 180d0/db
allocate(pos(3,xtc % natoms))
allocate(g(0:nrhist,nbhist))
g(:,:) = 0d0
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

if (trim(grp(1)%resname) == trim(grp(2)%resname)) then
  flag = ON
else
  flag = OFF
end if


num = nstep
if(nstep == 0) num = 1
print*, 'reading trajctory and calculating length...'
do step = 0,nstep 
  if (mod(step*10,num) == 0) print *, int(dble(step*100)/dble(num)),'%'
  pos(:,:) = xtc % pos(:,:)
!$omp parallel shared(pos, npos1, npos2,ncite, box_dim, list1, list2, nrhist,id_h) &
!$omp & private(j,k,id_1, id_2, ig,ib, v_oo, r_oo,v_oh,r_oh,beta) &
!$omp & reduction(+:g)
!$omp do
  do i = 1, npos1
    id_1(0) = list1(i)
    id_1(1) = 1 +  int(dble(id_1(0))/dble(mol(grp(1)%order)%natoms)) 
    do j = 1,  npos2
      id_2(0) = list2(j)
      id_2(1) = 1 +  int(dble(id_2(0))/dble(mol(grp(2)%order)%natoms)) 
      v_oo(:) = pos(:,id_2(0)) - pos(:,id_1(0))
      call length(box_dim, v_oo, r_oo)
      do k = 1, ncite
        v_oh = pos(:,id_1(0) + id_h(k)) - pos(:,id_1(0))
        call length(box_dim, v_oh, r_oh)
        beta = calc_degree(v_oo,v_oh,r_oo,r_oh,ON)
        ig = ceiling(r_oo / dr)
        ib = ceiling(beta / db)
        if (ig <= nrhist .and. id_1(flag) /= id_2(flag)) then
          g(ig,ib) = g(ig,ib) + 1
        end if
      end do
    end do
  end do
!$omp end do
!$omp end parallel
  call xtc % read
  call get_boxdim(xtc%box, box_dim)
end do

call xtc % close
print *, ' normalizing ...'
! normaize of RDF
maxg = 0d0

rho = dble(npos1)/(box_dim(1)*box_dim(2)*box_dim(3))
do i = 0, nrhist
  do j = 1, nbhist
    ! Delta V nm
    dv = (2d0/3d0)*pi*(i**3 - (i-1)**3)*dr**3*(-cos(j*pi/180d0) +cos((j-1)*pi/180d0))
    g(i,j) = g(i,j)/(dble(nstep+1)*ncite*dble(npos2)*dv)
    if (g(i,j) > maxg) then
      maxg = g(i,j)
    end if
  end do
end do


g(:,:) = g(:,:)/maxg


print*, 'outputting...'
outfile='Data_pmf_'//trim(arg1)//'_'//trim(arg2)//'.xvg'
open(15,file=outfile)
do i = 0, nrhist
  beta = 0d0
  do j = 1, nbhist
    beta = beta + 1d0*db
    if (g(i,j) /= 0d0) then
      write(15,*) (i -0.5)*dr, beta, -log(g(i,j))
    else
      write(15,*) '?'
    end if
  end do
  write(15,*) ''
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
  print *, './pmf a b'
  print *, ''
  print *, 'a : atom1 (donar)'
  print *, 'b : atom2 (acceptor)'
  print *, ''
  print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  print *, '|               !! YARINAOSE  !!                |'
  print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  stop

end subroutine warn_arg



end program pmf
