program freeenergysegments
use readandget
use netcdf
implicit none
character(len=50) :: infile, reffile, iname, rname, hname, sname, chrep, chat
logical :: netcdf,per,velin,startfound,endfound,normalize
integer :: i, j, k, at, nsteps, spatial, natoms, nrestr, rep, auxunit, auxunit2, bins, start, end, cutlenght
integer, allocatable, dimension (:) :: mask
real(4), allocatable, dimension (:) :: coordx, coordy, coordz
double precision, allocatable, dimension (:) :: hist_x, hist_y, segment_x, segment_y, singlecoord
double precision, allocatable, dimension (:) :: hist_x_cut, hist_y_cut, segment_x_cut, segment_y_cut
double precision, allocatable, dimension (:,:) :: rref
double precision, allocatable, dimension(:,:,:) :: coordall
double precision :: max, min, kref, temp, truncar, R, ini, fin
double precision, dimension(6) :: boxinfo

R=0.001987d0 !kcal/K mol
call readinputsegmentsfreenergy(nrestr,mask,infile,reffile,rep,netcdf,bins,kref,temp,per,velin,boxinfo,truncar,normalize)

allocate(hist_x(bins),hist_y(bins))
if (rep .le. 9) write(chrep,'(I1)') rep
if (rep .gt. 9 .and. rep .le. 99) write(chrep,'(I2)') rep
rname = trim(reffile) // "_"
rname = trim(rname) // trim(chrep)
rname = trim(rname) // ".rst7"

if (netcdf) then
  iname = trim(infile) // "_"
  iname = trim(iname) // trim(chrep)
  iname = trim(iname) // ".nc"
  call getdims(iname,nsteps,spatial,natoms)
  if (allocated(coordx)) deallocate(coordx)
  if (allocated(coordy)) deallocate(coordy)
  if (allocated(coordz)) deallocate(coordz)
  if (allocated(coordall)) deallocate(coordall)
  allocate(coordx(nsteps),coordy(nsteps),coordz(nsteps),singlecoord(nsteps))
  allocate(coordall(3,nrestr,nsteps),segment_x(nsteps),segment_y(nsteps),rref(3,natoms))
  call getcoordfromnetcdf(iname,nsteps,natoms,spatial,coordx,coordy,coordz,coordall,nrestr,mask)
else
  call getcoordfromfenebtraj(nsteps,coordall,nrestr,rep)
  allocate(singlecoord(nsteps),segment_x(nsteps),segment_y(nsteps))
endif

segment_x=0.d0
segment_y=0.d0
call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin)

do j=1,nrestr
  at=mask(j)
  do i=1,3
    do k=1,nsteps
       singlecoord(k)=coordall(i,j,k)
    end do
    min = MINVAL(singlecoord)
    max = MAXVAL(singlecoord)
    call histogram( singlecoord, nsteps, min, max, bins, hist_x, hist_y )

    k=1
    startfound=.false.
    endfound=.false.
    start=1
    end=bins
    do while (k .le. bins .and. (.not. startfound .or. .not. endfound))
      if (hist_y(k) .ge. truncar .and. .not. startfound) then
         start = k
         startfound = .true.
      endif
      if (hist_y(bins-k+1) .ge. truncar .and. .not. endfound) then
         end = bins-k+1
         endfound = .true.
      end if
    k = k + 1
    end do
    if (.not. startfound .or. .not. endfound) write(*,*) "WARNING, maybe histograms are looking weird"

    do k=start,end
      segment_x(k)=hist_x(k)
      segment_y(k)=-R*temp*LOG(hist_y(k))-(kref/2.d0)*(hist_x(k)-rref(i,at))**2
    end do

    cutlenght=end-start+1
    if (allocated(hist_x_cut)) deallocate(hist_x_cut)
    if (allocated(hist_y_cut)) deallocate(hist_y_cut)
    if (allocated(segment_x_cut)) deallocate(segment_x_cut)
    if (allocated(segment_y_cut)) deallocate(segment_y_cut)
    allocate(hist_x_cut(cutlenght), hist_y_cut(cutlenght), segment_x_cut(cutlenght), segment_y_cut(cutlenght))

    do k=1,cutlenght
      hist_x_cut(k)=hist_x(start+k-1)
      hist_y_cut(k)=hist_y(start+k-1)
      segment_x_cut(k)=segment_x(start+k-1)
      segment_y_cut(k)=segment_y(start+k-1)
    end do
    min=MINVAL(segment_y_cut)
    do k=1,cutlenght
      segment_y_cut(k)=segment_y_cut(k)-min
    end do

    if(normalize) then
      ini=segment_x_cut(1)
      fin=segment_x_cut(cutlenght)
      do k=1,cutlenght
        segment_x_cut(k)=(segment_x_cut(k)-ini)/(fin-ini)
      end do
      ! ini=MINVAL(segment_y_cut)
      ! fin=MAXVAL(segment_y_cut)
      ! do k=1,cutlenght
      !   segment_y_cut(k)=(segment_y_cut(k)-ini)/(fin-ini)
      ! end do
    end if



    if (i.eq.1) hname = "h.x.rep."//chrep
    if (i.eq.2) hname = "h.y.rep."//chrep
    if (i.eq.3) hname = "h.z.rep."//chrep
    if (i.eq.1) sname = "s.x.rep."//chrep
    if (i.eq.2) sname = "s.y.rep."//chrep
    if (i.eq.3) sname = "s.z.rep."//chrep
    if (at .le. 9) write(chat,'(I1)') at
    if (at .gt. 9 .and. at .le. 99) write(chat,'(I2)') at
    hname = trim(hname)//".at."
    hname = trim(hname)//trim(chat)
    auxunit = 1000000*i+10000*j+rep
    sname = trim(sname)//".at."
    sname = trim(sname)//trim(chat)
    auxunit2 = 99000000*i+990000*j+rep
    if (auxunit.eq.auxunit2) write(*,*) "WARNING, auxunit = auxunit2"
    open(auxunit, file = hname, status = 'unknown')
    open(auxunit2, file = sname, status = 'unknown')
    do k=1,cutlenght
        write(auxunit,*) hist_x_cut(k), hist_y_cut(k)
        write(auxunit2,*) segment_x_cut(k), segment_y_cut(k)
    end do
    close(auxunit)
    close(auxunit2)

  end do
end do

end program freeenergysegments



SUBROUTINE histogram( A, L, x_start, x_end, N, hist_x, hist_y )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Histogram takes an array A, of length L,
! starting and ending x values, x_start, x_end,
! and the number of bins, N, and outputs
! 2 arrays, hist_x and hist_y which contain the histogram
! data from A.
!
! extracted from https://searchcode.com/codesearch/raw/96378353/
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	IMPLICIT NONE

	!! Inputs
	INTEGER										:: L
	INTEGER										:: N
	REAL(KIND=8)							:: x_start
	REAL(KIND=8)							:: x_end
	REAL(KIND=8),DIMENSION(L)	:: A

	!! Outputs
	REAL(KIND=8),DIMENSION(N)	:: hist_x
	REAL(KIND=8),DIMENSION(N)	:: hist_y

	!! Internal
	INTEGER										:: i
	INTEGER										:: j
	INTEGER										:: tot
	REAL(KIND=8)							:: dx
	REAL(KIND=8)							:: left
	REAL(KIND=8)							:: right
  REAL(KIND=8)							:: normalization_factor

	!! Bin Width
	dx = (x_end - x_start)/DBLE(N)

	!! Loop over all bins
	DO i=1,N
		!! Set up hist_x, evenly spaced bin array with values of the bin centers
		hist_x(i) = x_start + (i-1)*dx + dx/2.0D0

		!! bin left value
		left      = hist_x(i) - dx/2.0D0
		!! bin right value
		right     = hist_x(i) + dx/2.0D0
		!! initialize counter for bin to 0
		tot       = 0
		!! Loop over entire array A
		DO j=1,L
			!! Increment counter if value is within bin left and right
			IF ( A(j) .GT. left .AND. A(j) .LT. right ) tot = tot+1
		END DO
		!! Set hist_y to bin counter
		hist_y(i) = DBLE(tot)
	END DO

! normalization (added by J. Semelak 2021)
  normalization_factor=0.d0
  DO i = 1, N
    normalization_factor = normalization_factor + hist_y(i)*dx
  END DO
  DO i = 1, N
    hist_y(i)=hist_y(i)/normalization_factor
  END DO

END SUBROUTINE histogram
