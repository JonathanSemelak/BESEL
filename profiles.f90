subroutine getprofile(rav,fav,nrep,nrestr,profile)
implicit none
double precision, dimension(3,nrestr,nrep), intent(in) :: fav, rav
double precision, dimension(2,nrestr,nrep-1) :: profileall
double precision, dimension(2,nrep-1), intent(out) :: profile
integer, intent(in) :: nrestr, nrep
integer :: i,j

profile=0.d0
profileall=0.d0

do i=1,nrestr
  do j=1,nrep-1
    profileall(1,i,j)=(dble(j)+dble(j+1))/2.d0
    profileall(2,i,j)=- (profileall(2,i,j) + ((fav(1,i,j)+fav(1,i,j+1))*(rav(1,i,j+1)-rav(1,i,j))/2.d0)  &
                                           + ((fav(2,i,j)+fav(2,i,j+1))*(rav(2,i,j+1)-rav(2,i,j))/2.d0)  &
                                           + ((fav(3,i,j)+fav(3,i,j+1))*(rav(3,i,j+1)-rav(3,i,j))/2.d0))
  end do
end do

do j=1,nrep-1
  do i=1,nrestr
    profile(1,j)=profileall(1,i,j)
    profile(2,j)=profile(2,j)+profileall(2,i,j)
  end do
end do

do j=2,nrep-1
  profile(2,j)=profile(2,j-1)+profile(2,j)
end do

profile(2,1:nrep-1)=profile(2,1:nrep-1)-profile(2,1)

open(unit=1659, file="profile.dat", position='append')
do j=1,nrep-1
  write(1659,'(2x, f4.1, 2x, f20.10)') profile(1,j), profile(2,j)
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
