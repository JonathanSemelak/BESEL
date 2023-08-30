program integrator
use netcdf
use readandget
implicit none
double precision, dimension (:,:,:), allocatable :: rav,fav
double precision, dimension (:,:), allocatable :: profile
integer :: i, j, k, trash, nlines, io, atnum, atnum_old, nrestr, nrep

open(unit=25, file="feneb.gradients")
nlines=0
atnum=0
atnum_old=0
nrestr=0
nrep=0
do
  read(25,*,iostat=io) atnum
  if (atnum .eq. atnum_old+1) then
    nrestr=nrestr+1
    atnum_old=atnum
  end if
  if (io/=0) exit
  nlines = nlines + 1
end do
nrep=nlines/nrestr
close(25)

allocate (rav(3,nrestr,nrep),fav(3,nrestr,nrep),profile(2,nrep))

open(unit=25, file="feneb.gradients")
do i=1,nrep
  do j=1,nrestr
    read(25,*) trash, rav(1:3,j,i), fav(1:3,j,i)
  end do
  read(25,*)
end do
close(25)

call getprofile(rav,fav,nrep,nrestr,profile)

end program integrator
