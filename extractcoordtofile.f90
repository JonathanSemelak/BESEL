program extractcoordtofile
use readandget
use netcdf
implicit none
character(len=500) :: infile, inname
character(len=500) :: chi, iname
integer :: i, j, k, ati, nsteps, spatial, natoms, nrestr
integer, allocatable, dimension (:) :: mask
real(4), allocatable, dimension (:) :: coordx,coordy,coordz
double precision, allocatable, dimension(:,:,:) :: coordall

call readinputextrator(nrestr,mask,infile,i)

if (i .le. 9) write(chi,'(I1)') i
if (i .gt. 9 .and. i .le. 99) write(chi,'(I2)') i
iname = trim(infile) // "_"
iname = trim(iname) // trim(chi)
iname = trim(iname) // ".nc"
chi = "feneb.traj."//chi
ati=mask(i)

call getdims(iname,nsteps,spatial,natoms)
if (allocated(coordx)) deallocate(coordx)
if (allocated(coordy)) deallocate(coordy)
if (allocated(coordz)) deallocate(coordz)
if (allocated(coordall)) deallocate(coordall)
allocate(coordx(nsteps),coordy(nsteps),coordz(nsteps))
allocate(coordall(3,nrestr,nsteps))

call getcoordfromnetcdf(iname,nsteps,natoms,spatial,coordx,coordy,coordz,coordall,nrestr,mask)

open (unit=2000, file=chi)
write(2000,*) nsteps
do j=1,nrestr
  do k=1,nsteps
    write(2000,*) coordall(1,j,k),coordall(2,j,k),coordall(3,j,k)
  end do
end do
close(unit=2000)
end program extractcoordtofile
