subroutine readinput(nrep,infile,reffile,outfile,mask,nrestr)

implicit none
character(len=50) :: infile, reffile, outfile
integer :: nrestr, nrep, i
logical ::  per, vel
double precision :: kref, kspring, steep_size, ftol
integer, allocatable, dimension (:) :: mask
double precision, allocatable, dimension(:,:,:) :: rav, fav, tang

open (unit=1000, file="feneb.in", status='old', action='read') !read align.in file
read(1000,*) infile
read(1000,*) reffile
read(1000,*) outfile
read(1000,*) per, vel
read(1000,*) nrep
read(1000,*) nrestr
if (nrep .eq. 1) read(1000,*) kref
if (nrep .gt. 1) read(1000,*) kref, kspring
if (nrep .gt. 1) allocate(tang(3,nrestr,nrep))
allocate(mask(nrestr),rav(3,nrestr,nrep),fav(3,nrestr,nrep))
read(1000,*) steep_size
read(1000,*) ftol
read(1000,*) (mask(i),i=1,nrestr)
close (unit=1000)

end subroutine readinput
