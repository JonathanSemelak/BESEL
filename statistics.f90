subroutine getsstatistics(coord,nsteps,skip,nevalfluc,dt,Z,H0,minsegmentlenght,goodrav)
implicit none
real (kind=4), DIMENSION(nsteps) :: coord
integer :: skip, nsteps, nskip, nmax, nevalfluc, nsegment, segmentlenght, minsegmentlenght
double precision, allocatable, dimension(:) :: segmentedcoord, segmenteddev, coordav
double precision, dimension(38) :: test
double precision ::  dt, S, m, Z, xi, xj, U, V, devav, goodrav
integer :: i,j,k,start,end,count,N
logical :: toosmall
logical, intent(inout) :: H0
! nmax=0
H0=.False.
toosmall=.False.
nskip=skip
! nskip=0
do while((.not.H0) .and. (.not.toosmall))
  ! write(*,*) (.not.H0), (.not.toosmall), ((.not.H0) .and. (.not.toosmall))
  if (allocated(coordav)) deallocate(coordav)
  ! write(*,*) nsteps, nskip, nsteps-nskip
  allocate(coordav(nsteps-nskip))
  coordav=0.d0
  do k=nskip+1,nsteps
    coordav(k-nskip)=(coordav(k-nskip-1)*(dble(k-1-nskip))+coord(k))/dble(k-nskip)
    write(8888,*) nskip, k-nskip, coord(k), coordav(k-nskip)
  end do
  ! write(8888,*)


if ((nsteps-nskip).le.nevalfluc) toosmall=.True.
if ((nsteps-nskip).le.nevalfluc) write(*,*) "Not enough data 1"

nmax=0
do i=2, nevalfluc-1
  ! write(*,*)  coord(i-1),coord(i-1),coord(i)
 if ((coordav(i) .gt. coordav(i-1)) .and. (coordav(i) .gt. coordav(i+1))) nmax=nmax+1

end do

if (nmax.eq.0) then
  segmentlenght=minsegmentlenght
else
  segmentlenght=nevalfluc/nmax
end if
write(9999,*) "Fluctuation time is:", dble(segmentlenght*dt)
! write(*,*) nevalfluc, nmax, segmentlenght

if (segmentlenght.lt.minsegmentlenght) segmentlenght=minsegmentlenght
nsegment=(nsteps-nskip)/segmentlenght

if (nsegment.gt.100) then
  nsegment=100
  segmentlenght=(nsteps-nskip)/nsegment
endif
write(9999,*) "Number of segments:", nsegment
write(9999,*) "Number of points in each segment:", segmentlenght
if (nsegment.lt.24) toosmall=.True.
if (nsegment.lt.24) write(*,*) "Not enough data 2"

if (allocated(segmentedcoord)) deallocate(segmentedcoord)
if (allocated(segmenteddev)) deallocate(segmenteddev)
allocate(segmentedcoord(nsegment),segmenteddev(nsegment))


!computes mean
do i=1,nsegment
   segmentedcoord(i)=0.d0
   start=(i-1)*segmentlenght+1
   end=start+segmentlenght-1
   do j=start,end
     count=1
     segmentedcoord(i)=(segmentedcoord(i)*(dble(count-1))+coordav(j))/dble(count)
     count=count+1
   end do
end do

!computes standard deviation
do i=1,nsegment
   segmenteddev(i)=0.d0
   start=(i-1)*segmentlenght+1
   end=start+segmentlenght-1
   do j=start,end
     segmenteddev(i)=segmenteddev(i)+(coordav(j)-segmentedcoord(i))**2
   end do
   segmenteddev(i)=dsqrt(segmenteddev(i)/(segmentlenght-1))
end do


!do Mann-Kendall test for trend
S=0.d0
N=nsegment

do i=1,N-1
   do j=i+1,N
     xi=segmentedcoord(i)
     xj=segmentedcoord(j)
     if ((xj-xi).gt.0.d0) S = S+1.d0
     if ((xj-xi).lt.0.d0) S = S-1.d0
     if ((xj-xi).eq.0.d0) S = S+0.d0
   end do
end do

if (S.gt.0.d0) m=-1.d0
if (S.lt.0.d0) m=1.d0
if (S.eq.0.d0) m=0.d0

V=dsqrt((N*(N-1))/18.d0)
V=V*dsqrt(dble(2*N+5))
Z=(S-m)/V
U=1.96d0
H0=.False.
! write(*,*) nskip,abs(Z)
if (abs(Z).lt.U) then
   H0=.True.
   write(*,*) "Gotta drop ", nskip, " steps"
else
  nskip=nskip+1000
endif
! write(*,*) H0, toosmall
! write(*,*) "ASD1"
end do
! write(*,*) "ASD2"

goodrav=0.d0
! write(*,*) "ASD3"

do i=1,nsegment
  goodrav=goodrav+segmentedcoord(i)
end do
! write(*,*) "ASD4"

goodrav=goodrav/dble(nsegment)

! write(*,*) "ASD5"

devav=0.d0
do j=1,nsteps-nskip
   devav=devav+(coordav(j)-goodrav)**2
end do
devav=dsqrt(devav/(nsteps-nskip-1))
! write(*,*) "ASD6"


write(6666,*) goodrav
write(7777,*) devav
! write(*,*) "ASD7"

end subroutine getsstatistics
