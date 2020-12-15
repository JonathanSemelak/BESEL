subroutine getmaxforce(nrestr,nrep,rep,fav,maxforce,ftol,relaxd,maxforceat,rms)
  implicit none
  double precision, dimension(3,nrestr,nrep), intent(in) :: fav
  integer, intent(in) :: nrep, rep, nrestr
  double precision, intent(in) :: ftol
  double precision, intent(inout) :: maxforce, rms
  logical, intent(out) :: relaxd
  double precision :: fmax2
  integer :: i,j,maxforceat

  relaxd=.FALSE.
  maxforce=0.d0
  maxforceat=0
  do i=1,nrestr
    fmax2=0.d0
    do j=1,3
      fmax2=fmax2+fav(j,i,rep)**2
    end do
    rms=rms+fmax2
    if (fmax2 .gt. maxforce) maxforce=fmax2
    if (fmax2 .gt. maxforce) maxforceat=i
  end do
  maxforce=dsqrt(maxforce)
  if(maxforce .le. ftol) relaxd=.TRUE.

end subroutine getmaxforce

subroutine steep(rav,fav,nrep,rep,steep_size,maxforce,nrestr,lastmforce,stepl,deltaA,dontg)
implicit none
double precision, dimension(3,nrestr,nrep), intent(inout) :: rav
double precision, dimension(3,nrestr,nrep) :: rnew
double precision, dimension(3,nrestr,nrep), intent(in) :: fav, dontg
double precision, intent(out) :: stepl
double precision :: lastmforce, steep_size, step, deltaA, n1, n2, n3
integer, intent(in) :: nrep, rep, nrestr
double precision, intent(inout) :: maxforce
integer :: i,j,auxunit
logical :: moved


stepl=steep_size

if (maxforce .lt. 1d-30) stepl=0.d0
step=stepl/maxforce

  do i=1,nrestr
    do j=1,3
      rav(j,i,rep)=rav(j,i,rep)+step*fav(j,i,rep)
    end do
  end do

  if (stepl .lt. 1d-10) then
    moved=.true.
    stepl=0.d0
    write(*,*) "Max precision reached"
  end if

end subroutine steep
