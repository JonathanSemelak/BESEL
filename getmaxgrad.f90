program getmaxgrad
implicit none
double precision, dimension (:,:,:), allocatable :: rav,fav
double precision, dimension (:,:), allocatable :: gradnorm
double precision, dimension (:), allocatable :: maxgradnorms
double precision :: maxgradnorm
integer :: i, j, k, trash, nlines, io, atnum, atnum_old, nrestr, nrep
integer :: maxgradreplica

open(unit=25, file="Pos_forces.dat")
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

allocate (rav(3,nrestr,nrep),fav(3,nrestr,nrep))

open(unit=25, file="Pos_forces.dat")
do i=1,nrep
  do j=1,nrestr
    read(25,*) trash, rav(1:3,j,i), fav(1:3,j,i)
  end do
  read(25,*)
end do
close(25)

allocate (gradnorm(nrestr,nrep),maxgradnorms(nrep))

do i=1,nrep
  maxgradnorms(i)=gradnorm(1,i)
  do j=2,nrestr
    if (gradnorm(j,i) .gt. maxgradnorms(i)) maxgradnorms(i)=gradnorm(j,i)
  end do
end do

maxgradnorm=maxgradnorms(i)
maxgradreplica=1
do i=1,nrep
    if (maxgradnorms(i) .gt. maxgradnorms(i)) then
      maxgradnorm=maxgradnorms(i)
      maxgradreplica=i
    end if
end do

write(*,*) maxgradnorm

end program getmaxgrad
