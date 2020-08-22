program feneb
implicit none
character(len=50) :: infile,  chi, fname2
character(len=51) ::fname
integer :: nsteps, spatial, natoms, nrestr, nrep
integer :: i, j
integer, allocatable, dimension (:) :: mask
real(4) :: coordinate
real(4), allocatable, dimension (:) :: coordinate_evol
integer, dimension (3) :: point
double precision :: kref
double precision, allocatable, dimension(:,:,:) :: rav, fav
!------------ Read input

  open (unit=55, file="aka.in", status='old', action='read') !read align.in file
  read(55,*) infile
  read(55,*) nrep
  read(55,*) nrestr
  read(55,*) kref
  allocate(mask(nrestr),rav(3,nrestr,nrep),fav(3,nrestr,nrep))
  read(55,*) (mask(i),i=1,nrestr)
  close (unit=55)
!------------

  open(unit=66, file="Pos_forces.dat")
  do i=1,nrep
    if (i .le. 9) write(chi,'(I1)') i
    if (i .gt. 9 .and. i .le. 99) write(chi,'(I2)') i
    fname = trim(infile) // "_"
    fname = trim(fname) // trim(chi)
    fname = trim(fname) // ".nc"
    call getdims(fname,nsteps,spatial,natoms)
    allocate(coordinate_evol(nsteps))
    call getandwritecoord(fname,nsteps,natoms,spatial,point,coordinate_evol,&
    nrestr,mask,kref,rav,fav,nrep,i)
    do j=1,nrestr
      write(66,'(2x, I6,2x, 6(f20.10,2x))') j, rav(1:3,j,i), fav(1:3,j,i)
    end do
    write(66,'(2x, I6,2x, 6(f20.10,2x))')
  end do
  close(66)


contains

SUBROUTINE getdims(infile,nsteps,spatial,natoms)
USE netcdf
IMPLICIT NONE
INTEGER(KIND=4), INTENT(OUT) :: nsteps,spatial,natoms
INTEGER(KIND=4) :: ncid
CHARACTER(LEN=50), INTENT(IN) :: infile
CHARACTER(LEN=50) :: xname, yname, zname

CALL check(nf90_open(infile, nf90_nowrite, ncid))
CALL check(nf90_inquire_dimension(ncid,1,xname,nsteps))
CALL check(nf90_inquire_dimension(ncid,2,yname,spatial))
CALL check(nf90_inquire_dimension(ncid,3,zname,natoms))
CALL check(nf90_close(ncid))
END SUBROUTINE getdims

SUBROUTINE getandwritecoord(infile,nsteps,natoms,spatial,point, &
                            coordinate_evol,nrestr,mask,kref,rav,fav,nrep,rep)
USE netcdf
IMPLICIT NONE
REAL (KIND=4), DIMENSION(3) :: coordinate_evol, ref, av
INTEGER(KIND=4), INTENT(IN) :: natoms, nsteps, spatial, nrestr
INTEGER(KIND=4) :: ncid, xtype, ndims, varid
CHARACTER(LEN=50), INTENT(IN) :: infile
CHARACTER(LEN=50) :: xname, vname, chi, chrep
double precision :: kref
double precision, dimension(3,nrestr,nrep) :: rav,fav
integer, dimension(3) :: point,endp
integer, dimension(nrestr) :: mask
integer :: i,j,k,ati,atf,nrep,rep


CALL check(nf90_open(infile, nf90_nowrite, ncid))
do i = 1,nrestr !natoms
  if (i .le. 9) write(chi,'(I1)') i
  if (i .gt. 9 .and. i .le. 99) write(chi,'(I2)') i
  if (rep .le. 9) write(chrep,'(I1)') rep
  if (rep .gt. 9 .and. rep .le. 99) write(chrep,'(I2)') rep
  chi = "at."//chi
  chi = trim(chi)//".rep."
  chi = trim(chi)//trim(chrep)
  open(i, file = chi, status = 'unknown')
  ati=mask(i)
  atf=mask(1)
  point = (/ 1,ati,1 /)
  endp = (/ 3,atf,1 /)
  CALL check(nf90_get_var(ncid,3,ref,start = point,count = endp))
  av=0.d0
  do j = 1,nsteps
      point = (/ 1,ati,j /)
      endp = (/ 3,atf,1 /)
      CALL check(nf90_get_var(ncid,3,coordinate_evol,start = point,count = endp))
      av(1:3)=(av(1:3)*(dble(j)-1)+coordinate_evol(1:3))/dble(j)
      write(i,*) j, -kref*(av(1:3)-ref(1:3))
  enddo
  rav(1:3,i,rep)=av(1:3)
  fav(1:3,i,rep)=-kref*(av(1:3)-ref(1:3))
  close(i)
enddo

CALL check(nf90_close(ncid))
END SUBROUTINE getandwritecoord

SUBROUTINE check(istatus)
USE netcdf
IMPLICIT NONE
INTEGER, INTENT (IN) :: istatus
IF (istatus /= nf90_noerr) THEN
write(*,*) TRIM(ADJUSTL(nf90_strerror(istatus)))
END IF
END SUBROUTINE check


end program feneb
