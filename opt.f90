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

end subroutine getmaxforce

subroutine steep(rav,fav,nrep,rep,steep_size,maxforce,nrestr,lastmforce,stepl)
implicit none
double precision, dimension(3,nrestr,nrep), intent(inout) :: rav
double precision, dimension(3,nrestr,nrep), intent(in) :: fav
double precision, intent(out) :: stepl
double precision :: lastmforce, steep_size, step
integer, intent(in) :: nrep, rep, nrestr
double precision, intent(inout) :: maxforce
integer :: i,j

stepl=steep_size

if (maxforce .gt. lastmforce) stepl=stepl*0.85d0
if (stepl .lt. 1d-5) then

  write(9999,*) "-----------------------------------------------------"
  write(9999,*) "Warning: max precision reached on atomic displacement"
  write(9999,*) "step length will be set to zero"
  write(9999,*) "-----------------------------------------------------"

  stepl=0.d0
end if
!write (*,*) rep, nrep, nrep-1, (rep .eq. nrep .or. rep .eq. nrep-1)
!if (rep .eq. nrep .or. rep .eq. nrep-1) write(9999,*) "Step length: ", step

if (maxforce .lt. 1d-30) then
  stepl=0.d0
else

step=stepl/maxforce
do i=1,nrestr
    do j=1,3
      rav(j,i,rep)=rav(j,i,rep)+step*fav(j,i,rep)
    end do
  end do
end if

end subroutine steep
