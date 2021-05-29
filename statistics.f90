subroutine getsstatistics(coord,nsteps,skip,nevalfluc,dt,Z,H0,minsegmentlenght,goodrav,gooddevav)
implicit none
real (kind=4), DIMENSION(nsteps) :: coord
integer :: skip, nsteps, nskip, nmax, nevalfluc, nsegment, segmentlenght, minsegmentlenght
double precision, allocatable, dimension(:) :: segmentedcoord, segmenteddev, coordav
double precision, dimension(38) :: test
double precision ::  dt, S, m, Z, xi, xj, U, V, gooddevav, goodrav
integer :: i,j,k,start,end,count,N
logical :: toosmall, nsegmentok
logical, intent(inout) :: H0
! nmax=0
H0=.False.
toosmall=.False.
nskip=skip

do while((.not.H0) .and. (.not.toosmall))
  if (allocated(coordav)) deallocate(coordav)
  allocate(coordav(nsteps-nskip))
  coordav=0.d0
  do k=nskip+1,nsteps
    coordav(k-nskip)=(coordav(k-nskip-1)*(dble(k-1-nskip))+coord(k))/dble(k-nskip)
  end do


  if ((nsteps-nskip).le.nevalfluc) toosmall=.True.
  if ((nsteps-nskip).le.nevalfluc) write(*,*) "Not enough data 1"

  nmax=0
  do i=2, nevalfluc-1
    if ((coordav(i) .gt. coordav(i-1)) .and. (coordav(i) .gt. coordav(i+1))) nmax=nmax+1
  end do

  if (nmax.eq.0) then
    segmentlenght=minsegmentlenght
  else
    segmentlenght=nevalfluc/nmax
  end if
  write(9999,*) "Fluctuation time is:", dble(segmentlenght*dt)

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

  nsegmentok=.False.
  N=24
  do while ((.not. nsegmentok) .and. (N.le.nsegment))
    S=0.d0
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
    if (abs(Z).lt.U) then
      H0=.True.
      nsegmentok=.True.
      ! write(*,*) N,nsegmentok
      ! write(*,*) "Gotta use:",N*segmentlenght, " steps (",N, "segments) Drop ", nskip, " steps"
    ! else
    !   nskip=nskip+1000
    endif
    N=N+1
  end do
  if (H0) then
    write(*,*) "Steps used for analysis:",N*segmentlenght
    write(*,*) "Segments:", N
    write(*,*) "Droped:", nskip
    write(*,*) "Total:", nskip+N*segmentlenght
  else
    nskip=nskip+1000
  endif
end do
N=N-1
do i=1,N
    write(8888,*) i, segmentedcoord(i)
end do
write(8888,*)

goodrav=0.d0
do i=1,N
  goodrav=goodrav+segmentedcoord(i)
end do
goodrav=goodrav/dble(nsegment)

gooddevav=0.d0
do j=1,N
   gooddevav=gooddevav+(segmentedcoord(j)-goodrav)**2
end do
gooddevav=dsqrt(gooddevav/(N-1))
write(7777,*) gooddevav


end subroutine getsstatistics
