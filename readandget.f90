module readandget
implicit none
contains
subroutine readinput(nrep,infile,reffile,outfile,mask,nrestr, &
           rav,fav,tang,kref,kspring,steep_size,ftol,per,vel)
implicit none
character(len=50) :: infile, reffile, outfile
integer :: nrestr, nrep, i
logical ::  per, vel
double precision :: kref, kspring, steep_size, ftol
integer, allocatable, dimension (:), intent(inout) :: mask
double precision, allocatable, dimension(:,:,:), intent(inout) :: rav, fav, tang

open (unit=1000, file="feneb.in", status='old', action='read') !read align.in
read(1000,*) infile
read(1000,*) reffile
read(1000,*) outfile
read(1000,*) per, vel
read(1000,*) nrep
read(1000,*) nrestr
if (nrep .eq. 1) read(1000,*) kref
if (nrep .gt. 1) read(1000,*) kref, kspring
if (nrep .gt. 1) allocate(tang(3,nrestr,nrep))
allocate(mask(nrestr),rav(3,nrestr,nrep),fav(3,nrestr,nrep))
read(1000,*) steep_size
read(1000,*) ftol
read(1000,*) (mask(i),i=1,nrestr)
close (unit=1000)

end subroutine readinput


subroutine getfilenames(rep,chrep,infile,reffile,outfile,iname,rname,oname)

implicit none
integer, intent(in) :: rep
character(len=50), intent(out) :: infile, reffile, outfile, chrep, iname, rname, oname

  if (rep .le. 9) write(chrep,'(I1)') rep
  if (rep .gt. 9 .and. rep .le. 99) write(chrep,'(I2)') rep
  if (rep .gt. 99 .and. rep .le. 999) write(chrep,'(I3)') rep

  iname = trim(infile) // "_"
  iname = trim(iname) // trim(chrep)
  iname = trim(iname) // ".nc"
  rname = trim(reffile) // "_"
  rname = trim(rname) // trim(chrep)
  rname = trim(rname) // ".rst7"
  oname = trim(outfile) // "_"
  oname = trim(oname) // trim(chrep)
  oname = trim(oname) // ".rst7"

end subroutine getfilenames

subroutine getdims(iname,nsteps,spatial,natoms)
use netcdf
implicit none
integer(kind=4), intent(out) :: nsteps,spatial,natoms
integer(kind=4) :: ncid
character(len=50), intent(in) :: iname
character(len=50) :: xname, yname, zname

call check(nf90_open(iname, nf90_nowrite, ncid))
call check(nf90_inquire_dimension(ncid,1,xname,nsteps))
call check(nf90_inquire_dimension(ncid,2,yname,spatial))
call check(nf90_inquire_dimension(ncid,3,zname,natoms))
call check(nf90_close(ncid))
end subroutine getdims

subroutine getavcoordanforces(iname,nsteps,natoms,spatial, &
                            coordx,coordy,coordz,nrestr,mask, &
                            kref,rav,fav,nrep,rep,rref)
use netcdf
implicit none
real (kind=4), DIMENSION(nsteps) :: coordx, coordy, coordz, ref, av
integer(kind=4), intent(in) :: natoms, nsteps, spatial, nrestr
integer(kind=4) :: ncid, xtype, ndims, varid
character(len=50), intent(in) :: iname
character(len=50) :: xname, vname, chi, chrep
double precision :: kref
double precision, dimension(3,nrestr,nrep) :: rav,fav
double precision, dimension(3,natoms), intent(inout) :: rref
integer, dimension(3) :: point,endp
integer, dimension(nrestr) :: mask
integer :: i,j,k,ati,atf,nrep,rep,auxunit

call check(nf90_open(iname, nf90_nowrite, ncid))
do i = 1,nrestr !natoms
  if (i .le. 9) write(chi,'(I1)') i
  if (i .gt. 9 .and. i .le. 99) write(chi,'(I2)') i
  if (rep .le. 9) write(chrep,'(I1)') rep
  if (rep .gt. 9 .and. rep .le. 99) write(chrep,'(I2)') rep
  chi = "at."//chi
  chi = trim(chi)//".rep."
  chi = trim(chi)//trim(chrep)
  auxunit = 1000*i+rep
  open(auxunit, file = chi, status = 'unknown')
  ati=mask(i)
  av=0.d0
  point = (/ 1,ati,1 /)
  endp = (/ 1,1,nsteps /)
  call check(nf90_get_var(ncid,3,coordx,start = point,count = endp))
  point = (/ 2,ati,1 /)
  endp = (/ 1,1,nsteps /)
  call check(nf90_get_var(ncid,3,coordy,start = point,count = endp))
  point = (/ 3,ati,1 /)
  endp = (/ 1,1,nsteps /)
  call check(nf90_get_var(ncid,3,coordz,start = point,count = endp))

  do k=1,nsteps
    av(1)=(av(1)*(dble(k-1))+coordx(k))/dble(k)
    av(2)=(av(2)*(dble(k-1))+coordy(k))/dble(k)
    av(3)=(av(3)*(dble(k-1))+coordz(k))/dble(k)
  write(auxunit,*) k, kref*(av(1:3)-rref(1:3,ati))
  end do
  rav(1:3,i,rep)=av(1:3)
  fav(1:3,i,rep)=kref*(av(1:3)-rref(1:3,ati))
  close(auxunit)
enddo

call check(nf90_close(ncid))
end subroutine getavcoordanforces


subroutine getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,vel)

implicit none
integer :: nrestr
integer, intent(out) :: natoms
character(len=50), intent(in) :: rname
double precision, dimension(:,:), allocatable :: rref
double precision, dimension(6), intent(out) :: boxinfo
integer, dimension(nrestr),intent(in) :: mask
integer :: i
logical ::  per, vel
if (allocated(rref)) deallocate(rref)
i=1
open (unit=1002, file=rname, status='old', action='read') !read ref file
read(1002,*)
read(1002,*) natoms
allocate(rref(3,natoms))
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

subroutine getcoordextrema(rref,natoms,rav,nrestr,nrep,rep,mask)
implicit none
double precision, dimension(3,nrestr,nrep), intent(out) :: rav
double precision, dimension(3,natoms), intent(in) :: rref
integer, intent(in) :: natoms, nrestr, nrep, rep
integer :: at,i,j
integer, dimension(nrestr) :: mask

  do i=1,nrestr
    do j=1,3
      at=mask(i)
      rav(j,i,rep) = rref(j,at)
    end do
  end do

end subroutine getcoordextrema

subroutine check(istatus)
use netcdf
implicit none
integer, intent (in) :: istatus
if (istatus /= nf90_noerr) then
write(9999,*) trim(adjustl(nf90_strerror(istatus)))
end if
end subroutine check

end module readandget
