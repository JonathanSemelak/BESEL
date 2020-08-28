program feneb
implicit none
character(len=50) :: infile, reffile, outfile, chi, iname, rname, oname
integer :: nsteps, spatial, natoms, nrestr, nrep
integer :: i, j
integer, allocatable, dimension (:) :: mask
real(4) :: coordinate
real(4), allocatable, dimension (:) :: coordx,coordy,coordz
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
  read(1000,*) outfile
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

 open(unit=9999, file="feneb.out") !Opten file for feneb output

!------------ Main loop

  do i=1,nrep

    call getfilenames(i,chi,infile,reffile,outfile,iname,rname,oname)

    call getdims(iname,nsteps,spatial,natoms)

    if (allocated(coordx)) deallocate(coordx)
    if (allocated(coordy)) deallocate(coordy)
    if (allocated(coordz)) deallocate(coordz)
    if (allocated(rref)) deallocate(rref)
    allocate(coordx(nsteps),coordy(nsteps),coordz(nsteps),rref(3,natoms))

    call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,vel)

    call getavcoordanforces(iname,nsteps,natoms,spatial,coordx,coordy, coordz,&
                    nrestr,mask,kref,rav,fav,nrep,i,rref)

    call writeposforces(rav,fav,nrestr,i)

    call getmaxforce(nrestr,nrep,i,fav,maxforce,ftol,relaxd)


  end do

!----------- End main loop
  if (nrep .gt. 1) then

  if (.not. relaxd) then
    call steep(rav,fav,nrep,i,steep_size,maxforce)
    call writenewcoord(oname,rref,boxinfo,natoms,nrestr,mask,per,rav,i)
  else
    write(999,*) "Convergence criteria of ", ftol, " (kcal/mol A) achieved"
  endif
   end if
  close(unit=9999)

contains

subroutine getfilenames(rep,chrep,infile,reffile,outfile,iname,rname,oname)

implicit none
integer, intent(in) :: rep
character(len=50), intent(out) :: infile, reffile, outfile, chrep, iname, rname, oname

  if (rep .le. 9) write(chi,'(I1)') rep
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

subroutine writeposforces(rav,fav,nrestr,rep)

implicit none
integer, intent(in) :: nrestr, rep
double precision, dimension(3,nrestr,nrep), intent(in) :: rav, fav

open(unit=1001, file="Pos_forces.dat")
do i=1,nrestr
write(1001,'(2x, I6,2x, 6(f20.10,2x))') rep, rav(1:3,i,rep), fav(1:3,i,rep)
end do
write(1001,'(2x, I6,2x, 6(f20.10,2x))')
close(1001)

end subroutine writeposforces

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
integer :: i,j,k,ati,atf,nrep,rep

call check(nf90_open(iname, nf90_nowrite, ncid))
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
  av=0.d0
  point = (/ 1,ati,1 /)
  endp = (/ 1,1,nsteps /)
  call check(nf90_get_var(ncid,3,coordx,start = point,count = endp))
  point = (/ 2,ati,1 /)
  endp = (/ 1,1,nsteps-1 /)
  call check(nf90_get_var(ncid,3,coordy,start = point,count = endp))
  point = (/ 3,ati,1 /)
  endp = (/ 1,1,nsteps-1 /)
  call check(nf90_get_var(ncid,3,coordz,start = point,count = endp))

  do k=1,nsteps
    av(1)=(av(1)*(dble(k-1))+coordx(k))/dble(k)
    av(2)=(av(2)*(dble(k-1))+coordy(k))/dble(k)
    av(3)=(av(3)*(dble(k-1))+coordz(k))/dble(k)
  write(i,*) k, kref*(av(1:3)-rref(1:3,ati))
  end do
  rav(1:3,i,rep)=av(1:3)
  fav(1:3,i,rep)=kref*(av(1:3)-rref(1:3,ati))
  close(i)
enddo

call check(nf90_close(ncid))
end subroutine getavcoordanforces


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
write(auxunit,'(I8)') natoms
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
     do j=1,3
       fmax2=fav(j,i,rep)**2
       if (fmax2 .gt. maxforce) maxforce=fmax2
    end do
  end do
  maxforce=dsqrt(maxforce)
  if(maxforce .le. ftol) relaxd=.TRUE.

  write(9999,*) "maxforce: ", maxforce


end subroutine getmaxforce

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
  do i=1,nrestr
    do j=1,3
      rav(j,i,rep)=rav(j,i,rep)+step*fav(j,i,rep)
    end do
  end do

end subroutine steep

subroutine check(istatus)
use netcdf
implicit none
integer, intent (in) :: istatus
if (istatus /= nf90_noerr) then
write(9999,*) trim(adjustl(nf90_strerror(istatus)))
end if
end subroutine check







end program feneb
