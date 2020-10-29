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
    fmax2=0.d0
    do j=1,3
      fmax2=fmax2+fav(j,i,rep)**2
    end do
    if (fmax2 .gt. maxforce) maxforce=fmax2
  end do
  maxforce=dsqrt(maxforce)
  if(maxforce .le. ftol) relaxd=.TRUE.

end subroutine getmaxforce

subroutine steep(rav,fav,nrep,rep,steep_size,maxforce,nrestr,lastmforce,stepl,deltaA)
implicit none
double precision, dimension(3,nrestr,nrep), intent(inout) :: rav
double precision, dimension(3,nrestr,nrep) :: rnew
double precision, dimension(3,nrestr,nrep), intent(in) :: fav
double precision, intent(out) :: stepl
double precision :: lastmforce, steep_size, step, deltaA, n1, n2
integer, intent(in) :: nrep, rep, nrestr
double precision, intent(inout) :: maxforce
integer :: i,j,auxunit
logical :: moved
stepl=steep_size
!moved=.false.
!if (maxforce .gt. lastmforce) stepl=stepl*0.85d0
!if (maxforce .lt. lastmforce) stepl=stepl*1.1d0
!if (nrep .gt. 1 .and. maxforce .gt. lastmforce) stepl=stepl*0.85d0
!if (stepl .lt. 1d-5 .or. maxforce .lt. 1d-30) stepl=0.d0
if (maxforce .lt. 1d-30) stepl=0.d0
step=stepl/maxforce
write(1810,*) rep, stepl, maxforce, step

!   do while (.not. moved)
!   deltaA=0.d0
!   step=stepl/maxforce

  do i=1,nrestr
    n1=dble(fav(1,i,rep)**2+fav(2,i,rep)**2+fav(3,i,rep)**2)
    n2=dble(rav(1,i,rep)**2+rav(2,i,rep)**2+rav(3,i,rep)**2)
    auxunit=3000+i
    write(auxunit,*) rep, n2, n1*step
    do j=1,3
      rav(j,i,rep)=rav(j,i,rep)+step*fav(j,i,rep)
    end do
  end do

  ! do i=1,nrestr
  !   do j=1,3
  !     deltaA=deltaA-fav(j,i,rep)*(rnew(j,i,rep)-rav(j,i,rep))
  !   end do
  ! end do
  !write(*,*) "deltaA", deltaA
  ! if (deltaA .lt. 0.d0) then
  !   rav=rnew
  !   moved=.true.
  ! else
  !   stepl=stepl*0.85
  ! end if

  if (stepl .lt. 1d-10) then
    moved=.true.
    stepl=0.d0
    write(*,*) "Max precision reached"
  end if
! end do

! do i=1,nrestr
!   do j=1,3
!     rav(j,i,rep)=rav(j,i,rep)+step*fav(j,i,rep)
!   end do
! end do


end subroutine steep
