subroutine getprofile(rav,fav,nrep,nrestr,profile)
implicit none
double precision, dimension(3,nrestr,nrep), intent(in) :: fav, rav
double precision, dimension(2,nrestr,nrep) :: profileall, proftest
double precision, dimension(2,nrep), intent(out) :: profile
integer, intent(in) :: nrestr, nrep
integer :: i,j

profile=0.d0
profileall=0.d0
do i=1,nrestr
  do j=1,nrep-1
    profileall(1,i,j)=j
    if (j.eq.1) then
      profileall(2,i,j)=0.d0
    else
      profileall(2,i,j)=- (((fav(1,i,j)+fav(1,i,j+1))*(rav(1,i,j+1)-rav(1,i,j))/2.d0)  &
                          + ((fav(2,i,j)+fav(2,i,j+1))*(rav(2,i,j+1)-rav(2,i,j))/2.d0)  &
                          + ((fav(3,i,j)+fav(3,i,j+1))*(rav(3,i,j+1)-rav(3,i,j))/2.d0))
    end if
  end do
end do

do j=1,nrep
  do i=1,nrestr
    profile(1,j)=profileall(1,i,j)
    profile(2,j)=profile(2,j)+profileall(2,i,j)
  end do
end do

do j=2,nrep
  profile(2,j)=profile(2,j-1)+profile(2,j)
end do

! profile(2,1:nrep-1)=profile(2,1:nrep-1)-profile(2,1)

open(unit=1659, file="profile.dat", position='append')
do j=1,nrep
  write(1659,'(2x, f5.1, 2x, f20.10)') profile(1,j), profile(2,j)
end do
close(1659)

end subroutine

subroutine getbarrier(profile, nrep, barrier, minpoint, maxpoint)
implicit none

double precision, dimension(2,nrep-1), intent(in) :: profile
integer, intent(in) :: nrep
double precision :: minim, barrier, minpoint, maxpoint
integer :: maxindex,j

minpoint=profile(1,1)
maxpoint=profile(1,1)
minim=profile(2,1)
barrier=profile(2,1)
maxindex=1

do j=2,nrep-1
  if (profile(2,j) .gt. barrier) then
    barrier=profile(2,j)
    maxpoint=profile(1,j)
    maxindex=j
  end if
end do


do j=2,maxindex
  if (profile(2,j) .lt. minim) then
    minim=profile(2,j)
    minpoint=profile(1,j)
  end if
end do

end subroutine

subroutine getrmsd(fav, kref, nrep, nrestr,rmsd)
implicit none
double precision, dimension(3,nrestr,nrep), intent(in) :: fav
double precision, dimension(nrep), intent(out) :: rmsd
double precision,intent(in) :: kref
integer, intent(in) :: nrestr, nrep
integer :: i,j,n

rmsd=0.d0
do n=1,nrep
  do i=1,nrestr
    do j=1,3
      rmsd(n)=rmsd(n)+(fav(j,i,n)/kref)**2
    end do
  end do
end do

rmsd=dsqrt(rmsd/dble(nrestr))

end subroutine

subroutine getselfdist(rav, rrefall, nrep, nrestr, selfdist)
implicit none
double precision, dimension(3,nrestr,nrep), intent(in) :: rav, rrefall
integer, intent(in) :: nrestr, nrep
double precision, dimension(2,nrestr,nrep-1), intent(out) :: selfdist
integer :: i,j,n

selfdist=0.d0
do n=1,nrep-1
  do i=1,nrestr
    do j=1,3
      selfdist(1,i,n)=selfdist(1,i,n)+(rrefall(j,i,n)-rrefall(j,i,n+1))**2
      selfdist(2,i,n)=selfdist(2,i,n)+(rav(j,i,n)-rav(j,i,n+1))**2
    end do
    selfdist(1,i,n)=dsqrt(selfdist(1,i,n))
    selfdist(2,i,n)=dsqrt(selfdist(2,i,n))
  end do
end do

end subroutine

subroutine getdistrightminusleft(rav, nrep, nrestr, equispaced)
implicit none
double precision, dimension(3,nrestr,nrep), intent(in) :: rav
double precision :: distright, distleft
integer, intent(in) :: nrestr, nrep
logical, intent(inout) :: equispaced
integer :: i,j

equispaced=.True.

do i=2,nrep-1
  do j=1,nrestr
      distright=(rav(1,j,i+1)-rav(1,j,i))**2+(rav(2,j,i+1)-rav(2,j,i))**2+(rav(3,j,i+1)-rav(3,j,i))**2
      distleft=(rav(1,j,i)-rav(1,j,i-1))**2+(rav(2,j,i)-rav(2,j,i-1))**2+(rav(3,j,i)-rav(3,j,i-1))**2
  end do
  distright=sqrt(distright)
  distleft=sqrt(distleft)
  equispaced=(equispaced .and. (abs(distright-distleft) .lt. 0.0001d0))
end do

end subroutine getdistrightminusleft

subroutine geterror(rav,fav,nrep,nrestr,profile)
implicit none
double precision, dimension(3,nrestr,nrep), intent(in) :: fav, rav
double precision, dimension(2,nrestr,nrep) :: profileall, proftest
double precision, dimension(2,nrep), intent(out) :: profile
integer, intent(in) :: nrestr, nrep
integer :: i,j

profile=0.d0
profileall=0.d0
do i=1,nrestr
  do j=1,nrep-1
    profileall(1,i,j)=j
    if (j.eq.1) then
      profileall(2,i,j)=0.d0
    else
      profileall(2,i,j)=- (((fav(1,i,j)+fav(1,i,j+1))*(rav(1,i,j+1)-rav(1,i,j))/2.d0)  &
                          + ((fav(2,i,j)+fav(2,i,j+1))*(rav(2,i,j+1)-rav(2,i,j))/2.d0)  &
                          + ((fav(3,i,j)+fav(3,i,j+1))*(rav(3,i,j+1)-rav(3,i,j))/2.d0))
    end if
  end do
end do

do j=1,nrep
  do i=1,nrestr
    profile(1,j)=profileall(1,i,j)
    profile(2,j)=profile(2,j)+profileall(2,i,j)
  end do
end do

do j=2,nrep
  profile(2,j)=profile(2,j-1)+profile(2,j)
end do

! profile(2,1:nrep-1)=profile(2,1:nrep-1)-profile(2,1)

open(unit=1660, file="error.dat", position='append')
do j=1,nrep
  write(1660,'(2x, f5.1, 2x, f20.10)') profile(1,j), profile(2,j)
end do
close(1660)

end subroutine geterror
