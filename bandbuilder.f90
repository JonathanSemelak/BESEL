	program bandbuilder
!
!This program generates an initial set of replicas for NEB calculations in Amber .rst7 format,
!using only reactives and products restarts (and TS if possible)
!Original subroutine by N. Foglia 03/2018
!Adapted to amber restars by J. Semelak 09/2020
use readandget
use netcdf
implicit none
character(len=50) :: rcfile, pcfile, tsfile, prefix, chi, oname
integer :: nrestr, nrep, i, j, k, middlepoint, natoms
logical ::  usets, per, velin, velout, test
double precision, dimension(3) :: BAND_slope
double precision, dimension(3) :: BAND_const
double precision, dimension(6):: boxinfo
integer, allocatable, dimension (:) :: mask
double precision, allocatable, dimension(:,:) :: rclas, selfdist, rref
double precision, allocatable, dimension(:,:,:) :: rav


!reads imputfile
call readinputbuilder(rcfile, pcfile, tsfile, prefix, nrestr, nrep, usets, per, velin, velout, rav, mask, test)
if (.not. test) then

if (usets) then
	!reads middlepoint coordinates
	if (allocated(rclas)) deallocate(rclas)
	call getrefcoord(tsfile,nrestr,mask,natoms,rclas,boxinfo,per,velin)
	middlepoint=nrep+1
	middlepoint=middlepoint/2
	do i=1,nrestr
		rav(1:3,i,middlepoint)=rclas(1:3,mask(i))
	end do
else
	middlepoint=nrep
end if

!reads product complex coordinates
if (allocated(rclas)) deallocate(rclas)
call getrefcoord(pcfile,nrestr,mask,natoms,rclas,boxinfo,per,velin)

do i=1,nrestr
	rav(1:3,i,nrep)=rclas(1:3,mask(i))
end do

!reads reactant complex coordinates
if (allocated(rclas)) deallocate(rclas)
call getrefcoord(rcfile,nrestr,mask,natoms,rclas,boxinfo,per,velin)

do i=1,nrestr
	rav(1:3,i,1)=rclas(1:3,mask(i))
end do


!generate sinitial middleimages
do i=1,nrestr
  BAND_slope(1:3)= rav(1:3,i,middlepoint)-rav(1:3,i,1)
  BAND_slope=BAND_slope/(dble(middlepoint) - 1.d0)
  BAND_const=rav(1:3,i,1)-BAND_slope(1:3)
  do k=1, middlepoint
    rav(1:3,i,k)=BAND_slope(1:3)*dble(k) + BAND_const(1:3)
  end do

!usign TS state case
  if (middlepoint .ne. nrep) then
    BAND_slope(1:3)= rav(1:3,i,nrep) - rav(1:3,i,middlepoint)
    BAND_slope=BAND_slope/(dble(nrep) - dble(middlepoint))
    BAND_const=rav(1:3,i,middlepoint) - dble(middlepoint)*BAND_slope(1:3)
    do k=middlepoint, nrep
	    rav(1:3,i,k)=BAND_slope(1:3)*dble(k) + BAND_const(1:3)
	  end do
	end if
end do

!write output files

do i=1,nrep
			if (i .le. 9) write(chi,'(I1)') i
			if (i .gt. 9 .and. i .le. 99) write(chi,'(I2)') i
	    if (i .gt. 99 .and. i .le. 999) write(chi,'(I3)') i
			oname = trim(prefix) // "_"
  		oname = trim(oname) // trim(chi)
			oname = trim(oname) // ".rst7"
	    call writenewcoord(oname,rclas,boxinfo,natoms,nrestr,mask,per,velout,rav,nrep,i)
end do

allocate(selfdist(nrestr,nrep-1))
call selfdistiniband(rav, nrep, nrestr, selfdist)

open(unit=1111, file="selfdist.dat")
do i=1,nrestr
	do j=1,nrep-1
		write(1111,'(2x, I6,2x, f20.10)') j, selfdist(i,j)
	end do
	write(1111,*)
end do
!-------------TESTTTTTTTTTTTTTTTT


else
	allocate(rref(3,natoms))
	do i=1,nrep
		if (i .le. 9) write(chi,'(I1)') i
		if (i .gt. 9 .and. i .le. 99) write(chi,'(I2)') i
		if (i .gt. 99 .and. i .le. 999) write(chi,'(I3)') i
		oname = trim(prefix) // "_"
		oname = trim(oname) // trim(chi)
		oname = trim(oname) // ".rst7"
    ! call getfilenames(i,chi,infile,reffile,outfile,iname,rname,oname) !rname = NAME_r_i.rst7 ; i=replica
	  call getrefcoord(oname,nrestr,mask,natoms,rref,boxinfo,per,velin)
	  call getcoordextrema(rref,natoms,rav,nrestr,nrep,i,mask)
  end do
	allocate(selfdist(nrestr,nrep-1))
	call selfdistiniband(rav, nrep, nrestr, selfdist)

	open(unit=1111, file="selfdist.dat")
	do i=1,nrestr
		do j=1,nrep-1
			write(1111,'(2x, I6,2x, f20.10)') j, selfdist(i,j)
		end do
		write(1111,*)
	end do


end if!END IF TEST
end program bandbuilder


subroutine selfdistiniband(rav, nrep, nrestr, selfdist)

implicit none
double precision, dimension(3,nrestr,nrep), intent(in) :: rav
integer, intent(in) :: nrestr, nrep
double precision, dimension(nrestr,nrep-1), intent(out) :: selfdist
integer :: i,j,n

selfdist=0.d0
do n=1,nrep-1
  do i=1,nrestr
    do j=1,3
      selfdist(i,n)=selfdist(i,n)+(rav(j,i,n)-rav(j,i,n+1))**2
    end do
    selfdist(i,n)=dsqrt(selfdist(i,n))
  end do
end do
end subroutine
