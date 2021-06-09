!Original subroutines by N. Foglia 03/2018
!Oriented to free energy optimizations by J. Semelak 09/2020

subroutine gettang(rav,tang,nrestr,nrep)
implicit none
double precision, dimension(3,nrestr,nrep), intent(in) :: rav
double precision, dimension(3,nrestr,nrep), intent(out) :: tang
integer, intent(in) :: nrestr, nrep
double precision, dimension(nrestr,nrep) :: norm
integer :: i,j

tang=0.d0
norm=0.d0
do i=2,nrep-1
  do j=1,nrestr
    tang(1:3,j,i) = rav(1:3,j,i+1) - rav(1:3,j,i-1)
    norm(j,i) = tang(1,j,i)**2 + tang(2,j,i)**2 + tang(3,j,i)**2
    norm(j,i) = dsqrt(norm(j,i))
    if (norm(j,i) .lt. 1d-5) then
      tang(1:3,j,i) = 0.d0
    else
      tang(1:3,j,i) = tang(1:3,j,i)/norm(j,i)
    endif
  end do
end do

end subroutine gettang


subroutine getnebforce(rav,devav,fav,tang,nrestr,nrep,kspring,maxforceband,ftol,&
                       relaxd,ftrue,ftang,fperp,fspring,wrmforce,dontg)
implicit none
double precision, dimension(3,nrestr,nrep), intent(inout) :: fav
double precision, dimension(3,nrestr,nrep), intent(in) :: rav, devav, tang
double precision, intent(out) :: maxforceband
integer :: maxforcerep,maxforceat,maxstdrep,maxstdat
integer, intent(in) :: nrestr, nrep
double precision, intent(in) :: kspring, ftol
double precision, dimension(3,nrestr,nrep) :: fspring, ftang, fperp, ftrue, dontg
double precision, dimension(nrestr,nrep) :: fproj
double precision :: distright, distleft, maxforce, rms, maxforceband2,n1,n2,n3,n4,n5
double precision :: maxstdband,maxstd
integer :: i,j,auxunit
logical :: relaxdrep,relaxd,wrmforce

	fspring=0.d0
  fproj=0.d0
  ftang=0.d0
  ftrue=0.d0
  maxforceband=0.d0
  maxstdband=0.d0
  relaxd=.FALSE.
  do i=2,nrep-1
	  do j=1, nrestr
        !Saves true force
        ftrue(1:3,j,i)=fav(1:3,j,i)
        !Computes spring force
	      distright=(rav(1,j,i+1)-rav(1,j,i))**2+(rav(2,j,i+1)-rav(2,j,i))**2+(rav(3,j,i+1)-rav(3,j,i))**2
        distleft=(rav(1,j,i)-rav(1,j,i-1))**2+(rav(2,j,i)-rav(2,j,i-1))**2+(rav(3,j,i)-rav(3,j,i-1))**2
        distright=sqrt(distright)
	      distleft=sqrt(distleft)
	      fspring(1:3,j,i)=kspring*(distright-distleft)*tang(1:3,j,i)

        !Computes force component on tangent direction
        fproj(j,i)=fav(1,j,i)*tang(1,j,i)+fav(2,j,i)*tang(2,j,i)+fav(3,j,i)*tang(3,j,i)
        !Computes neb force
        ftang(1:3,j,i)=fproj(j,i)*tang(1:3,j,i)
        fperp(1:3,j,i)=fav(1:3,j,i)-ftang(1:3,j,i)
        fav(1:3,j,i)=fperp(1:3,j,i)+fspring(1:3,j,i)

	  end do
    if (wrmforce) then
      call getmaxforce(nrestr,nrep,i,fperp,maxforce,ftol,relaxdrep,maxforceat,rms)
      if (maxforce .gt. maxforceband) maxforcerep=i
      if (maxforce .gt. maxforceband) maxforceband=maxforce
      write(9999,*) "Replica: ", i, "Max force: ", maxforce, "Converged: ", relaxdrep
      call getmaxstd(nrestr,nrep,i,fperp,devav,maxstd,maxstdat)
      if (maxstd .gt. maxstdband) maxstdrep=i
      if (maxstd .gt. maxstdband) maxstdband=maxstd
    else
      call getmaxforce(nrestr,nrep,i,fspring,maxforce,ftol,relaxdrep,maxforceat,rms)
      if (maxforce .gt. maxforceband) maxforcerep=i
      if (maxforce .gt. maxforceband) maxforceband=maxforce
    end if
  end do
  if (wrmforce) then
    write(9999,*) "-----------------------------------------------------------------"
    write(9999,*) "Band max force: ", maxforceband, "on replica: ", maxforcerep
    write(9999,*) "-----------------------------------------------------------------"
    write(9999,*) "Band max STD: ", maxstdband, "on replica: ", maxstdrep
    write(9999,*) "-----------------------------------------------------------------"
    if(maxforceband .le. ftol) relaxd=.TRUE.
    if (relaxd) write(9999,*) "System converged: T"
    if (.not. relaxd) write(9999,*) "System converged: F"
  endif

!Test ----------------
  maxforceband2=0.d0
  do i=2,nrep-1
    if (wrmforce) then
      call getmaxforce(nrestr,nrep,i,fav,maxforce,ftol,relaxdrep,maxforceat,rms)
      if (maxforce .gt. maxforceband2) maxforcerep=i
      if (maxforce .gt. maxforceband2) maxforceband2=maxforce
    end if
  end do
  if (wrmforce) then
    write(9999,*) "-----------------------------------------------------------------"
    write(9999,*) "Band max fneb: ", maxforceband2, "on replica: ", maxforcerep
    write(9999,*) "-----------------------------------------------------------------"
  endif

  maxforceband2=0.d0
  do i=2,nrep-1
    if (wrmforce) then
      call getmaxforce(nrestr,nrep,i,fspring,maxforce,ftol,relaxdrep,maxforceat,rms)
      if (maxforce .gt. maxforceband2) maxforcerep=i
      if (maxforce .gt. maxforceband2) maxforceband2=maxforce
    end if
  end do

  if (wrmforce) then
    write(9999,*) "Band max fspring0: ", maxforceband2, "on replica: ", maxforcerep
    write(9999,*) "-----------------------------------------------------------------"
  endif
!---------------------

end subroutine getnebforce
