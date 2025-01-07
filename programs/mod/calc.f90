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

module calc
use info
implicit none
contains

subroutine get_boxdim(inp_dim, box_dim)
  implicit none
  real, intent(in) :: inp_dim(3,3)
  real, intent(out) :: box_dim(3)

  box_dim(1) = inp_dim(1,1)
  box_dim(2) = inp_dim(2,2)
  box_dim(3) = inp_dim(3,3)

end subroutine get_boxdim

subroutine length(box_dim, v, r)
  implicit none
  real, intent(in) :: box_dim(3)
  real*8, intent(inout) :: v(3)
  real*8 :: r

  v(:) = v(:) - box_dim(:)*nint(v(:)/box_dim(:))
  r = sqrt(sum(v**2))

end subroutine length


real*8 function calc_degree(v1,v2,r1,r2,flag)
implicit none
  real*8, parameter :: pi = dacos(-1d0)
  integer, intent(in) :: flag
  real*8, intent(in) ::v1(3), v2(3)
  real*8, intent(in) ::r1, r2
  real*8 :: cosb,beta
  cosb = (sum(v1*v2))/(r1*r2)
  if ( cosb > 1d0 ) then
    cosb = 1d0
  elseif ( cosb < -1d0 ) then
    cosb = -1d0
  end if
  beta = dacos(cosb)
  if (flag == 1) then
    calc_degree = 180.0d0*beta/pi
  else
    calc_degree = beta
  end if
  
end function calc_degree

subroutine centering(mol, pos, posg, box_dim)
!subroutine centering(mol, pos, posg)
  implicit none
  type(MOL_atom), allocatable,intent(in) :: mol(:)
  real*8, allocatable, intent(in) :: pos(:,:)
  real*8, intent(out) :: posg
  real, intent(in) :: box_dim
  integer :: i, j, k
  integer :: posatm
  real*8 :: sum_mass, z

  posg = 0d0
  sum_mass = 0d0

  do i = 1, size(mol)
    if ( mol(i) % flag_com == 1 ) then
      posatm = 0
      do j = 1, mol(i) % order - 1
        posatm = posatm +  (mol(j) % natoms)*(mol(j) % nmol)
      end do
      do j = 1, mol(i) % nmol
        do k = 1, mol(i) % natoms
          posatm = posatm + 1
          z = pos(3, posatm) !- box_dim*nint(pos(3, posatm)/box_dim)
          if (abs(z) > box_dim) z = z - box_dim*nint(z/box_dim)
          posg = posg + z*(mol(i) % mass(k))
          !posg = posg + pos(3, posatm)*(mol(i) % mass(k))
          sum_mass = sum_mass + mol(i) % mass(k)
        end do
      end do
    end if
  end do

  posg = posg / sum_mass

end subroutine centering

integer function calc_logstep(nstep, logn)
  implicit none
  integer, intent(in) ::  nstep
  real*8, intent(in) :: logn
  integer :: logstep
  logstep = int(log(real(nstep))/log(logn))
  print*, 'logstep = ', logstep
  calc_logstep = logstep
end function calc_logstep

subroutine calc_com(mol, pos,box_dim)
  implicit none
  type(MOL_atom), intent(inout) :: mol
  real*8, intent(in) :: pos(:,:)
  real , intent(in) :: box_dim(:)
  real*8 :: first_pos(3), cood(3)
  integer :: j, k, l, ig, posatm
  ig  = 0
    mol%pos_com(:,:) = 0d0
  do j = 1, mol%nmol
    ig = mol%natoms*(j-1) + 1
    posatm = mol%list(ig)
    first_pos = pos(:,posatm)
    mol%pos_com(:,j) = mol%pos_com(:,j) + pos(:,posatm)*mol%mass(1)
    do k = 2, mol%natoms
      ig = mol%natoms*(j-1) + k
      posatm = mol%list(ig)
      cood(:) = pos(:,posatm)
      do l = 1, 3
        if ( abs(pos(l,posatm) - first_pos(l)) > (box_dim(l))/2d0 ) then
          if ( pos(l,posatm) > box_dim(l)/2d0 ) then
            cood(l) = pos(l,posatm) - box_dim(l)*nint(pos(l,posatm)/box_dim(l))
          else
            cood(l) = pos(l,posatm) + box_dim(l)*nint(pos(l,posatm)/box_dim(l))
          end if
        end if
      end do
      mol%pos_com(:,j) = mol%pos_com(:,j) + cood(:)*mol%mass(k)
    end do
    !print*, pos_com
    mol%pos_com(:,j) = mol%pos_com(:,j)/sum(mol%mass)
  end do
end subroutine calc_com


end module calc

