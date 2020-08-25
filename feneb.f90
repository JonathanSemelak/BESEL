program feneb
implicit none
character(len=50) :: infile, reffile, ofile, chi, fname, rname, oname
integer :: nsteps, spatial, natoms, nrestr, nrep
integer :: i, j
integer, allocatable, dimension (:) :: mask
real(4) :: coordinate
real(4), allocatable, dimension (:) :: coordinate_evol
integer, dimension (3) :: point
double precision :: kref, steep_size, ftol, maxforce
double precision, dimension(6) :: boxinfo
double precision, allocatable, dimension(:,:) :: rref
double precision, allocatable, dimension(:,:,:) :: rav, fav
logical ::  per, vel, relaxd
!------------ Read input

  open (unit=1000, file="feneb.in", status='old', action='read') !read align.in file
  read(1000,*) infile
  read(1000,*) reffile
  read(1000,*) ofile
  read(1000,*) per, vel
  read(1000,*) nrep
  read(1000,*) nrestr
  read(1000,*) kref
  read(1000,*) steep_size
  read(1000,*) ftol
  allocate(mask(nrestr),rav(3,nrestr,nrep),fav(3,nrestr,nrep))
  read(1000,*) (mask(i),i=1,nrestr)
  close (unit=1000)
!------------

  open(unit=1001, file="Pos_forces.dat")
  do i=1,nrep
    if (i .le. 9) write(chi,'(I1)') i
    if (i .gt. 9 .and. i .le. 99) write(chi,'(I2)') i
    fname = trim(infile) // "_"
    fname = trim(fname) // trim(chi)
    fname = trim(fname) // ".nc"
    rname = trim(reffile) // "_"
    rname = trim(rname) // trim(chi)
    rname = trim(rname) // ".rst7"
    oname = trim(ofile) // "_"
    oname = trim(oname) // trim(chi)
    oname = trim(oname) // ".rst7"

    call getdims(fname,nsteps,spatial,natoms)
    if (allocated(coordinate_evol)) deallocate(coordinate_evol)
    if (allocated(rref)) deallocate(rref)
    allocate(coordinate_evol(nsteps),rref(3,natoms))
    call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,vel)

    call getandwritecoord(fname,nsteps,natoms,spatial,coordinate_evol,&
    nrestr,mask,kref,rav,fav,nrep,i)

    do j=1,nrestr
      write(1001,'(2x, I6,2x, 6(f20.10,2x))') j, rav(1:3,j,i), fav(1:3,j,i)
    end do
    write(1001,'(2x, I6,2x, 6(f20.10,2x))')

    call max_force(nrestr,nrep,i,fav,maxforce,ftol,relaxd)

    if (.not. relaxd) then
      call steep(rav,fav,nrep,i,steep_size,maxforce)
      call writenewcoord(oname,rref,boxinfo,natoms,nrestr,mask,per,rav,i)
    else
      write(*,*) "Convergence criteria of ", ftol, " (kcal/mol A) achieved"
    endif
  end do

  close(1001)


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

SUBROUTINE getandwritecoord(infile,nsteps,natoms,spatial, &
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


subroutine getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,vel)

implicit none
integer :: nrestr, natoms
character(len=50), intent(in) :: rname
double precision, dimension(3,natoms), intent(inout) :: rref
double precision, dimension(6) :: boxinfo
integer, dimension(nrestr),intent(in) :: mask
integer :: i
logical ::  per, vel

i=1
open (unit=1002, file=rname, status='old', action='read') !read ref file
read(1002,*)
read(1002,*)
do while (i .le. natoms/2)
  read(1002,'(6(f12.7))') rref(1,2*i-1), rref(2,2*i-1), rref(3,2*i-1), &
                          rref(1,2*i), rref(2,2*i), rref(3,2*i)
  i = i + 1
enddo
if (mod(natoms,2) .ne. 0) read(1002,'(6(f12.7))') rref(2*i,1:3)
if (vel) then
  i=1
  do while (i .le. natoms/2)
    read(1002,*)
    i = i + 1
  enddo
  if (mod(natoms,2) .ne. 0) read(1002,*)
endif
if (per) read(1002,'(6(f12.7))') boxinfo(1:6)
close (unit=1002)
end subroutine getrefcoord

subroutine writenewcoord(oname,rref,boxinfo,natoms,nrestr,mask,per,rav,rep)

implicit none
character(len=50), intent(in) :: oname
integer, intent(in) :: natoms, nrestr, rep
double precision, dimension(3,nrestr,nrep) :: rav
double precision, dimension(3,natoms), intent(in) :: rref
double precision, dimension(3,natoms) :: rout
double precision, dimension(6), intent(in) :: boxinfo
integer, dimension(nrestr), intent(in) :: mask
logical, intent(in) :: per
integer :: i, j, at, auxunit

rout=rref
auxunit=2000+rep

do i=1,nrestr
  do j=1,3
    at=mask(i)
    rout(j,at)=rav(j,i,rep)
  end do
end do

open (unit=auxunit, file=oname)
write(auxunit,*) "FENEB restart, replica: ", rep
write(auxunit,*) natoms
i=1
do while (i .le. natoms/2)
  write(auxunit,'(6(f12.7))') rout(1,2*i-1), rout(2,2*i-1), rout(3,2*i-1), &
                          rout(1,2*i), rout(2,2*i), rout(3,2*i)
  i = i + 1
enddo
if (mod(natoms,2) .ne. 0) write(auxunit,'(6(f12.7))') rout(2*i,1:3)
i=1
do while (i .le. natoms/2)
  write(auxunit,'(6(f12.7))') 0.d0,0.d0,0.d0,0.d0,0.d0,0.d0
  i = i + 1
enddo
  if (mod(natoms,2) .ne. 0) write(auxunit,'(6(f12.7))') 0.d0, 0.d0, 0.d0
  if (per) write(auxunit,'(6(f12.7))') boxinfo(1:6)
close (unit=auxunit)

end subroutine writenewcoord


subroutine max_force(nrestr,nrep,rep,fav,maxforce,ftol,relaxd)
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

end subroutine max_force

subroutine steep(rav,fav,nrep,rep,steep_size,maxforce)
implicit none
double precision, dimension(3,nrestr,nrep), intent(inout) :: rav
double precision, dimension(3,nrestr,nrep), intent(in) :: fav
integer, intent(in) :: nrep, rep
double precision, intent(in) :: steep_size
double precision, intent(inout) :: maxforce
double precision :: step
integer :: i,j


step=steep_size/maxforce
  write(*,*) "maxforce: ", maxforce, "step size"
  do i=1,nrestr
    do j=1,3
      rav(j,i,rep)=rav(j,i,rep)+step*fav(j,i,rep)
    end do
  end do

end subroutine steep

SUBROUTINE check(istatus)
USE netcdf
IMPLICIT NONE
INTEGER, INTENT (IN) :: istatus
IF (istatus /= nf90_noerr) THEN
write(*,*) TRIM(ADJUSTL(nf90_strerror(istatus)))
END IF
END SUBROUTINE check







end program feneb
