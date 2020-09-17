subroutine getmaxforce(nrestr,nrep,rep,fav,maxforce,ftol,relaxd)
  implicit none
  double precision, dimension(3,nrestr,nrep), intent(in) :: fav
  integer, intent(in) :: nrep, rep, nrestr
  double precision, intent(in) :: ftol
  double precision, intent(out) :: maxforce
  logical, intent(out) :: relaxd
  double precision :: fmax2
  integer :: i,j

  relaxd=.FALSE.
  maxforce=0.d0
  do i=1,nrestr
     do j=1,3
       fmax2=fav(j,i,rep)**2
       if (fmax2 .gt. maxforce) maxforce=fmax2
    end do
  end do
  maxforce=dsqrt(maxforce)
  if(maxforce .le. ftol) relaxd=.TRUE.

  write(9999,*) "maxforce: ", maxforce

end subroutine getmaxforce

subroutine steep(rav,fav,nrep,rep,steep_size,maxforce,nrestr)
implicit none
double precision, dimension(3,nrestr,nrep), intent(inout) :: rav
double precision, dimension(3,nrestr,nrep), intent(in) :: fav
integer, intent(in) :: nrep, rep, nrestr
double precision, intent(in) :: steep_size
double precision, intent(inout) :: maxforce
double precision :: step
integer :: i,j

if (maxforce .lt. 1d-30) then
  step=0.d0
else
  step=steep_size/maxforce
end if
do i=1,nrestr
    do j=1,3
      rav(j,i,rep)=rav(j,i,rep)+step*fav(j,i,rep)
    end do
  end do

end subroutine steep
