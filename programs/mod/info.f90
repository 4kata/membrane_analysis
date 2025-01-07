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


module info
implicit none

type MOL_atom
  character(5) :: resname
  character(5), allocatable :: atmname(:)
  integer :: natoms, nmol, order
  integer, allocatable :: list_bonds(:,:), nbonds(:)
  integer*1, allocatable :: flag_heavy(:)
  integer*1 :: flag_com
  real*8, allocatable :: mass(:), pos_com(:,:)
  integer, allocatable :: list(:)
  real*8, allocatable :: sigma(:)
end type MOL_atom
type group_atom
  character(50) :: gname
  character(5)  :: resname
  integer :: order, natoms
  integer, allocatable :: list(:), index_list(:)
end type group_atom

contains
subroutine mdinfo(trjfile,nstep,dt,flag)

integer, intent(in) :: flag
integer, intent(out) :: nstep
real*8, intent(out) :: dt
character(50), intent(out) :: trjfile
character(50) :: nojumpfile, wrapfile

namelist /system_info/nojumpfile, wrapfile, nstep, dt

  open (17,file='MDinfo', status = 'old')
  read(17,nml=system_info)
  close(17)
  if ( flag == 0) then
    trjfile = wrapfile
  else if (flag == 1) then
    trjfile = nojumpfile
  end if
  
  print *, ''
  print *, '= MDinfo ='
  write(*,'(a11, i8)') ' flames   : ', nstep
  write(*,'(a11, f9.3)') ' dt [ps] : ', dt
  write(*,'(a16, f10.2)') ' Total t [ns] : ', nstep*dt/1000
  
end subroutine mdinfo

subroutine molinfo(mol,order)
  implicit none
  type(MOL_atom), allocatable, intent(inout) :: mol(:)
  type(MOL_atom) ::  sol
  character(5), allocatable, intent(inout)  :: order(:)
  character(50) :: line, key, flag
  integer :: i
  
  sol % flag_com = 0
  open (17,file='MOLinfo', status = 'old')
  read(17,*) i
  allocate(mol(i))
  allocate(order(i))
  do 
    read(17,'(A50)',end = 999) line
    key = '['
    if (index(line,trim(key)) /= 0) then
      read(17,*) sol % order, sol % resname, sol % natoms, sol % nmol, flag
      if (index(flag,trim('COM')) /= 0) then
        sol % flag_com = 1
      end if
      mol(sol%order) = sol
      order(sol % order) = sol % resname
      print*,  trim(sol%resname),' : ', sol%nmol, 'molecules'
    end if
  end do
999  close(17)


end subroutine molinfo


subroutine grpinfo(order, mol, grp)
  implicit none
  character(5), intent(in) :: order(:)
  type(MOL_atom), intent(in) :: mol(:)
  type(group_atom), intent(inout) :: grp
  integer :: i,j, num
  character(50) :: line, key
  character(1000) :: aaa

  key = grp % gname
  print*, 'making group list...'
  open(17,file='IDinfo')
  do
    read(17,'(A50)',end=999) line
    if (index(line,trim(key)) /= 0) then
      read(17,*) grp%resname
      read(17,*) grp % natoms
      write(*,'(a14,a3,i3,a6) ') grp%gname, ' : ',grp% natoms,' atoms'
      allocate(grp % list(grp % natoms))
      read(17,*) grp % list(:)
      exit
    end if
  end do
999  close(17)

  do i = 1, size(order)
    if (trim(order(i)) == trim(grp%resname)) grp % order = i
  end do

  if ( grp % order == 0) then
    print*, '!!!ERROR!!!'
    print*, 'resname of IDinfo is wrong!!!'
    stop
  end if

  num = grp%natoms * mol(grp%order)%nmol
  allocate(grp%index_list(num))
  grp%index_list(:) = 0
  
  if ( grp % order > 1) then
    do i = 1, grp % order - 1
      grp%index_list(:) = grp%index_list(:) + mol(i)%natoms * mol(i) % nmol
    end do
  end if
  do i = 1, mol(grp%order)%nmol
    do j = 1, grp%natoms
      num = (i-1) * grp%natoms + j
      grp%index_list(num) = grp%index_list(num) + (i-1)*mol(grp%order)%natoms + grp%list(j)
    end do
  end do

end subroutine grpinfo

subroutine trgtinfo(resname, mol, trgt)
  implicit none
  character(5), intent(in) :: resname
  type(MOL_atom), intent(in) :: mol(:)
  type(MOL_atom), intent(out) :: trgt
  integer, parameter :: ON = 1, OFF = 0
  integer :: flag, i
  
  print*, 'detecting target...'
  flag = OFF
  do i = 1, size(mol)
    if (trim(mol(i)%resname) .eq. trim(resname)) then
      trgt = mol(i)
      flag = ON
      exit
    end if
  end do
  if ( flag == OFF) call warn_resname()
end subroutine trgtinfo

subroutine inipos(mol,trgt, ini_pos)
  implicit none
  type(MOL_atom), intent(in) :: trgt, mol(:)
  integer, intent(out) :: ini_pos
  integer :: i
  print*, 'searching 1st atom of the resname...'
  ini_pos = 0
  if (trgt%order /= 1) then
    do i = 1, trgt%order - 1
    ini_pos = ini_pos + mol(i)%natoms*mol(i)%nmol
    end do
  end if
end subroutine inipos

subroutine bondinfo(mol)
  implicit none
  type(MOL_atom), intent(inout) :: mol
  integer, parameter :: YES = 1
  integer, parameter :: NO = 0
  integer :: i, j
  integer :: ai, aj
  character(4) :: atmname
  character(5) :: resname
  character(50) :: itpfile, key, line
   
  print*, 'making bond list...'

  allocate(mol%flag_heavy(mol%natoms))
  allocate(mol%mass(mol%natoms))
  allocate(mol%nbonds(mol%natoms))
  allocate(mol%list_bonds(4,mol%natoms))
  mol%flag_heavy(:) = YES
  mol%nbonds(:) = 0
  mol%list_bonds(:,:) = 0
  resname = mol%resname
  print *, 'open itp file ...'
  itpfile = 'toppar/'//trim(resname)//'.itp'
  open(17,file=itpfile,status='old')
  key = 'atoms'
  do
    read(17,'(A)') line
    if (index(line,trim(key))/=0) then
      read(17,*) 
      exit
    end if
  end do
  print *, 'defining atom types...'
  key = 'H'
  do i = 1, mol%natoms
    read(17,'(13x,a4,47x,f7.4)') atmname, mol%mass(i)
    do j = 1, len(atmname) 
      if (index(trim(atmname(j:j)),trim(key)) /= 0) then
        mol%flag_heavy(i) = NO
        exit
      else if (index(trim(atmname(j:j)),'') == 0) then
        exit
      end if
    end do
  end do
  key = 'bonds'
  print *, 'defining bond pair...'
  do
    read(17,'(A)')  line
    if (index(line,trim(key)) /= 0) then
      read(17,'(A)') line 
      exit
    end if
  end do
  key = 'pairs'
  do
    read(17,*) ai, aj, i
    mol%nbonds(ai) = mol%nbonds(ai) + 1
    mol%nbonds(aj) = mol%nbonds(aj) + 1
    mol%list_bonds(mol%nbonds(ai),ai) = aj
    mol%list_bonds(mol%nbonds(aj),aj) = ai
    read(17,*)
    read(17,'(A)') line
    if ( ai == 3 .or. aj ==3) then
      print *,  mol%list_bonds(mol%nbonds(ai),ai), mol%list_bonds(mol%nbonds(aj),aj)
    end if
  
    if (index(line,trim(key)) /= 0) then
      exit
    else
      do j =1,2
        backspace 17
      end do
    end if
  end do

999 close(17)
end subroutine bondinfo 

subroutine water_info(mol)
  implicit none
  type(MOL_atom), intent(inout) :: mol
  integer, parameter :: ON = 1
  integer, parameter :: OFF = 0
  
  allocate (mol%flag_heavy(3))
  allocate (mol%nbonds(3))
  allocate (mol%list_bonds(3,4))
  mol%list_bonds(:,:) = 0
  print*, size(mol%flag_heavy)
  mol%flag_heavy(1) = ON
  mol%flag_heavy(2) = OFF
  mol%flag_heavy(3) = OFF
  mol%nbonds(1) = 2
  mol%nbonds(2) = 1
  mol%nbonds(3) = 1
  mol%list_bonds(1,1)  = 2
  mol%list_bonds(2,1)  = 3
  mol%list_bonds(1,2)  = 1
  mol%list_bonds(1,3)  = 1
end subroutine water_info


integer function calc_posatm(ini_pos, i, natoms, id)
  implicit none
  integer, intent(in) :: ini_pos, i, natoms, id
  calc_posatm = ini_pos + (i-1)*natoms + id
end function calc_posatm



subroutine warn_resname()
  implicit none
  print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  print *, '|                 !! ERROR !!                  |'
  print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  print *, ''
  print *, '! Maybe wrong resname... Check again!'
  print *, ''
  print *, './msd a'
  print *, ''
  print *, 'a : resname'
  print *, ''
  print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  print *, '|               !! YARINAOSE !!                |'
  print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  stop

end subroutine warn_resname
end module info









