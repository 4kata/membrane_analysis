program z_profile
use info
use calc
use xdr, only: xtcfile
implicit none
type(xtcfile) :: xtc

!======================================================================!
!============================== variable ==============================!
!======================================================================!
! parameter
real, parameter :: dz = 0.01e0 !nm
! MDinfo
integer :: nstep
real*8 :: dt
character(50) :: trjfile
! MOLinfo
type(MOL_atom),allocatable :: mol(:)
character(5), allocatable :: order(:)
! GRPinfo
type(group_atom) :: grp
integer, allocatable :: list(:)
! get_boxdim
real :: box_dim(3)
! allocatable
real*8, allocatable :: pos(:,:)
real*8, allocatable :: g(:)
! others
integer :: i, j, k, l
integer :: step, num
integer :: nhist
integer :: posatm, ig
real :: r_max, dumy
real*8 :: z
real*8 :: sum_box(3), dv
real*8 :: posg
character(50) :: argv, outfile

!======================================================================!
print *, '*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*'
print*, ''
print *, '                      == z_profile =='
print *, ''

call  mdinfo(trjfile,nstep,dt,1)
call  molinfo(mol, order)

call getarg(1,argv)
grp % gname = argv
call grpinfo(order, mol, grp)
allocate(list(size(grp%index_list)))
list = grp%index_list

call xtc % init(trjfile)
call xtc % read
call get_boxdim(xtc % box, box_dim)
do i = 1, size(mol)
  if ( mol(i) % flag_com == 1 ) then
    call bondinfo(mol(i))
  end if
end do
r_max = box_dim(3)/2d0
nhist = ceiling(r_max/dz)
allocate(g(0:nhist-1))
allocate(pos(3,xtc % natoms))
g(:) = 0d0
sum_box(:) = 0d0

num = nstep
if(nstep == 0) num = 1

do step = 1, nstep + 1
  if (mod(step*10,num) == 0) print *, int(dble(step*100)/dble(num)),'%'
  pos(:,:) = xtc % pos(:,:)
!  call centering(mol, pos, posg)
  call centering(mol, pos, posg,box_dim(3))
!$omp parallel shared(pos,grp,mol,box_dim,nhist,list,posg) &
!$omp & private(z, posatm, ig) &
!$omp & reduction(+:g)
!$omp do
  do i = 1, grp%natoms * mol(grp%order)%nmol
    posatm = list(i)
    !z = abs(pos(3,posatm))
    if ( posg >= box_dim(3)/2d0 ) then
      z = pos(3,posatm) !- box_dim(3)*nint(pos(3, posatm)/box_dim(3))
      if ( abs(z) > box_dim(3)) z = z - box_dim(3)*nint(z/box_dim(3))
      z = z - posg 
      if ( abs(z) > box_dim(3)/2d0 ) z = z + box_dim(3)
    else
      z = pos(3,posatm) !- box_dim(3)*nint(pos(3, posatm)/box_dim(3))
      if ( abs(z) > box_dim(3)) z = z - box_dim(3)*nint(z/box_dim(3))
      z = z - posg 
      if ( abs(z) > box_dim(3)/2d0 ) z = z - box_dim(3)
    end if
    z = abs(z)
      

    !ig = nint(z*10d0/dz)
    !ig = int(real(ig)/10d0)
      !!ig = ceiling(real(ig)/10d0)
    ig = int(z/dz)
    !if ( ig == 1 .or. ig == 0) print*, z/dz , ig
    if (ig < nhist) then
      g(ig) = g(ig) + 1d0
    end if
  end do
!$omp end do
!$omp end parallel
  sum_box(:) = sum_box(:) + box_dim(:)
  call xtc % read 
  call get_boxdim(xtc % box, box_dim)
end do
sum_box(:) = sum_box(:)/dble(nstep + 1)
dv = sum_box(1)*sum_box(2)*dble(dz)
print*,sum_box(:),dv
!g(:) = g(:)/(dv*dble(((nstep + 1) * (grp%natoms * mol(grp%order)%nmol))))
g(:) = g(:)/(2d0*dv*dble((nstep + 1)))

outfile='Data_zp_'//trim(argv)//'.xvg'
open(15,file=outfile)
do i = -(nhist-1), nhist -1
  write(15,'(f7.3,f20.6)') i*dz, g(abs(i))
!  print*,  i,(i -0.5)*dz, g(i)
end do



end program z_profile
