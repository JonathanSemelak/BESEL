program freeenergysegments
use readandget
use netcdf
implicit none
character(len=50) :: infile, chrep, iname, hname, chat
logical :: netcdf
integer :: i, j, k, at, nsteps, spatial, natoms, nrestr, rep, auxunit, bins
integer, allocatable, dimension (:) :: mask
real(4), allocatable, dimension (:) :: coordx, coordy, coordz
double precision, allocatable, dimension (:) :: hist_x, hist_y, singlecoord
double precision, allocatable, dimension(:,:,:) :: coordall
double precision :: max, min, kref, temp

call readinputsegmentsfreenergy(nrestr,mask,infile,rep,netcdf,bins,kref,temp)

allocate(hist_x(bins),hist_y(bins))
if (rep .le. 9) write(chrep,'(I1)') rep
if (rep .gt. 9 .and. rep .le. 99) write(chrep,'(I2)') rep
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
  allocate(coordall(3,nrestr,nsteps))
  call getcoordfromnetcdf(iname,nsteps,natoms,spatial,coordx,coordy,coordz,coordall,nrestr,mask)
else
  call getcoordfromfenebtraj(nsteps,coordall,nrestr,rep)
  allocate(singlecoord(nsteps))
endif

do j=1,nrestr
  do i=1,3
    do k=1,nsteps
       singlecoord(k)=coordall(i,j,k)
    end do
    min = MINVAL(singlecoord)
    max = MAXVAL(singlecoord)
    call histogram( singlecoord, nsteps, min, max, bins, hist_x, hist_y )

    if (i.eq.1) hname = "h.x.rep."//chrep
    if (i.eq.2) hname = "h.y.rep."//chrep
    if (i.eq.3) hname = "h.z.rep."//chrep
    at=mask(j)
    if (at .le. 9) write(chat,'(I1)') at
    if (at .gt. 9 .and. at .le. 99) write(chat,'(I2)') at
    hname = trim(hname)//".at."
    hname = trim(hname)//trim(chat)
    auxunit = 1000000*i+10000*j+rep
    open(auxunit, file = hname, status = 'unknown')
    do k=1,bins
      write(auxunit,*) hist_x(k), hist_y(k)
    end do
    close(auxunit)
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
