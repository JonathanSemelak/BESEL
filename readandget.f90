module readandget
implicit none
contains
subroutine readinput(nrep,infile,reffile,outfile,topfile,mask,nrestr, &
           rav,ravout,devav,fav,ftrue,ftang,fperp,fspring,tang,kref,kspring,steep_size,steep_spring,ftol,per, &
           velin,velout,wgrad,wtemp,dt,wtempstart,wtempend,wtempfrec,mass,rrefall,nscycle,dontg,ravprevsetp, &
           rextrema, skip, dostat, minsegmentlenght, nevalfluc,rfromtraj,usensteps,nstepsexternal,smartstep, &
           typicalneb, tangoption, optoption, FIRE_dt_max, tangrecalc, maxdist, stopifconverged)
implicit none
character(len=50) :: infile, reffile, outfile, line, exp, keyword, topfile, prefix
integer :: nrestr, nrep, i, ierr, nscycle,skip, wtempfrec, wtempstart, wtempend
integer :: minsegmentlenght, nevalfluc, nstepsexternal, tangoption, optoption
logical :: per, velin, velout, wgrad, rextrema, wtemp, dostat, rfromtraj, usensteps, smartstep, typicalneb, tangrecalc
logical :: stopifconverged
double precision :: kref, kspring, steep_size, steep_spring, ftol, dt, FIRE_dt_max, maxdist
integer, allocatable, dimension (:), intent(inout) :: mask
double precision, allocatable, dimension(:,:,:), intent(inout) :: rav, fav, tang, ftang, ftrue,fperp, rrefall, ravprevsetp
double precision, allocatable, dimension(:,:,:), intent(inout) :: fspring, dontg, devav, ravout
double precision, allocatable, dimension(:), intent(inout) :: mass
integer :: start_range, end_range, i_mask, mask_index, comment_loc
character(len=20) :: substr
character(len=20), dimension(2) :: range

! set some default variables
 nscycle=50000
 rextrema=.False.
 steep_spring=0.001d0
 steep_size=0.01d0
 smartstep=.True.
 skip=0
 typicalneb=.False.
 wtemp=.False.
 dostat=.False.
 rfromtraj=.False.
 usensteps=.False.
 nstepsexternal=5000
 wtempfrec=1
 dt=0.001
 minsegmentlenght=100
 nevalfluc=1000
 tangoption=1
 optoption=0
 FIRE_dt_max=1.5d0
 tangrecalc=.True.
 maxdist=0.001d0
 stopifconverged=.False.
 wgrad=.False.
 per=.True.
 velin=.False.
 velout=.False.
 kspring=500.d0
 ftol=2.25d0

open (unit=1000, file='feneb.in', status='old', action='read') !read feneb.in
do
   read (1000,"(a)",iostat=ierr) line ! read line into character variable
   if (ierr /= 0) exit
   ! Trim line at the first "!" character to remove comments
   comment_loc = index(line, "!")
   if (comment_loc /= 0) then
       line = line(1:comment_loc-1)
   end if
   read (line,*,iostat=ierr) keyword ! read first keyword of line
   if (ierr /= 0) cycle
   if (keyword == 'prefix') read(line,*) exp,prefix
   if (keyword == 'per') read(line,*) exp, per
   if (keyword == 'velin') read(line,*) exp, velin
   if (keyword == 'velout') read(line,*) exp, velout
   if (keyword == 'nrep') read(line,*) exp, nrep
   if (keyword == 'nrestr') read(line,*) exp, nrestr
   if (keyword == 'kref') read(line,*) exp, kref
   if (keyword == 'kspring') read(line,*) exp, kspring
   if (keyword == 'steepsize') read(line,*) exp, steep_size
   if (keyword == 'steepspring') read(line,*) exp, steep_spring
   if (keyword == 'ftol') read(line,*) exp, ftol
   if (keyword == 'wgrad') read(line,*) exp, wgrad
   if (keyword == 'nscycle') read(line,*) exp, nscycle
   if (keyword == 'rextrema') read(line,*) exp, rextrema
   if (keyword == 'skip') read(line,*) exp, skip
   if (keyword == 'dt') read(line,*) exp, dt
   if (keyword == 'wtemp') read(line,*) exp, wtemp, wtempstart, wtempend, wtempfrec
   if (keyword == 'dostat') read(line,*) exp, dostat, nevalfluc, minsegmentlenght
   if (keyword == 'rfromtraj') read(line,*) exp, rfromtraj
   if (keyword == 'usensteps') read(line,*) exp, usensteps, nstepsexternal
   if (keyword == 'smartstep') read(line,*) exp, smartstep
   if (keyword == 'typicalneb') read(line,*) exp, typicalneb
   if (keyword == 'tangoption') read(line,*) exp, tangoption
   if (keyword == 'optoption') read(line,*) exp, optoption
   if (keyword == 'dtmax') read(line,*) exp, FIRE_dt_max
   if (keyword == 'tangrecalc') read(line,*) exp, tangrecalc
   if (keyword == 'maxdist') read(line,*) exp, maxdist
   if (keyword == 'stopifconverged') read(line,*) exp, stopifconverged
end do
close (unit=1000)
if (nrep .gt. 1) allocate(tang(3,nrestr,nrep),ftang(3,nrestr,nrep),ftrue(3,nrestr,nrep),&
                          fperp(3,nrestr,nrep),fspring(3,nrestr,nrep),dontg(3,nrestr,nrep))
allocate(mass(nrestr),mask(nrestr),rav(3,nrestr,nrep),fav(3,nrestr,nrep),rrefall(3,nrestr,nrep),&
        ravprevsetp(3,nrestr,nrep),devav(3,nrestr,nrep),ravout(3,nrestr,nrep))

open (unit=1000, file="feneb.in", status='old', action='read') !read feneb.in now that mask is allocated
mask_index = 1
do
    read (1000,"(a)",iostat=ierr) line ! read line into character variable
    if (ierr /= 0) exit
           read (line,*) keyword ! read first keyword of line
    if (keyword == 'mask') then
      line=line(6:) ! remove the first five characters which are 'mask'
      substr = ''
    do i = 1, len_trim(line)+1
        if (line(i:i) == ' ' .or. i==len_trim(line)+1) then
                 ! Here you have completed a "word" and can handle it as a single number
                 ! or the start of a range
          if (index(substr, '-') > 0) then
            call get_range(substr, range) ! parse range of numbers
            read(range(1), *) start_range
            read(range(2), *) end_range
            do i_mask = start_range, end_range
                if (mask_index > nrestr) then
                  print *, "Error: Too many indices for 'mask'. Increase 'nrestr'."
                  stop
                end if
                mask(mask_index) = i_mask
                mask_index = mask_index + 1
            end do
          else
            read(substr, *) i_mask
            if (mask_index > nrestr) then
                print *, "Error: Too many indices for 'mask'. Increase 'nrestr'."
                stop
            end if
            mask(mask_index) = i_mask
            mask_index = mask_index + 1
          end if
          substr = ''  ! Reset substr for the next "word"
      else
          substr = trim(substr) // line(i:i)
      end if
    end do
    end if
    if (keyword == 'all') then
      do i=1,nrestr
        mask(i)=i
      end do
    end if
end do
close (unit=1000)

! Finally, defines files names that use to be part of the input file
! In the future this should be cleanner...

infile = trim(prefix) // "_f"
reffile = trim(prefix) // "_r"
outfile = trim(prefix) // "_o"
topfile = trim(prefix) // ".prmtop"

end subroutine readinput

subroutine get_range(input_string, output_array)
   implicit none
   character(len=*), intent(in) :: input_string
   character(len=20), dimension(2), intent(out) :: output_array
   integer :: pos

   pos = index(input_string, '-')
   output_array(1) = input_string(:pos-1)
   output_array(2) = input_string(pos+1:)
end subroutine get_range

subroutine readtop(topfile,natoms,mask,mass,nrestr)
implicit none
character(len=50) :: topfile
double precision, dimension(natoms) :: massall
double precision, dimension(nrestr) :: mass
integer, dimension(nrestr) :: mask
integer :: natoms,nrestr, i, at, ierr
character(len=10) :: buffer
character(len=4)  :: ismass

open (unit=80008, file=topfile, status='old', action='read')
do
  read(80008,*,iostat=ierr) buffer, ismass
  if (ierr /= 0) exit
    if (ismass == "MASS") exit
enddo
read (80008,*)
i=1
do while (i .le. (natoms/5)*5)
  read(80008,*) massall(i), massall(i+1), massall(i+2), massall(i+3), massall(i+4)
  i = i + 5
enddo
  if (mod(natoms,5) .ne. 0) read(80008,*) massall(i:natoms)
close(80008)

do i=1,nrestr
  at=mask(i)
  mass(i)=massall(at)
end do

end subroutine readtop


subroutine readinputbuilder(rcfile, pcfile, tsfile, prefix, nrestr, nrep, usets, per, velin, velout,&
  rav, mask, iddp, nmax, onlytest, wselfdist)
implicit none
character(len=200) :: rcfile, pcfile, tsfile, prefix, exp, keyword, line, all
integer :: nrestr, nrep, i, ierr, nmax
logical ::  usets, per, velin, velout, onlytest, iddp, wselfdist
integer, allocatable, dimension (:), intent(inout) :: mask
double precision, allocatable, dimension(:,:,:), intent(inout) :: rav
integer :: start_range, end_range, i_mask, mask_index
character(len=20) :: substr
character(len=20), dimension(2) :: range
onlytest = .false.
usets = .false.
iddp = .false.
nmax = 2500
tsfile = 'None'
iddp = .False.
nmax = 500
wselfdist = .False.
per=.True.
velin=.False.
velout=.False.
open (unit=1000, file="bandbuilder.in", status='old', action='read') !read align.in
do
   read (1000,"(a)",iostat=ierr) line ! read line into character variable
   if (ierr /= 0) exit
   read (line,*) keyword ! read first keyword of line
   if (keyword == 'prefix') read(line,*) exp,prefix
   if (keyword == 'rcfile') read(line,*) exp,rcfile
   if (keyword == 'pcfile') read(line,*) exp,pcfile
   if (keyword == 'usets') read(line,*) exp,usets
   if (keyword == 'tsfile') read(line,*) exp,tsfile
   if (keyword == 'nrep') read(line,*) exp, nrep
   if (keyword == 'nrestr') read(line,*) exp, nrestr
   if (keyword == 'per') read(line,*) exp, per
   if (keyword == 'velin') read(line,*) exp, velin
   if (keyword == 'velout') read(line,*) exp, velout
   if (keyword == 'onlytest') read(line,*) exp, onlytest
   if (keyword == 'nmax') read(line,*) exp, nmax
   if (keyword == 'iddp') read(line,*) exp, iddp
   if (keyword == 'wselfdist') read(line,*) exp, wselfdist
end do

close (unit=1000)

allocate(mask(nrestr),rav(3,nrestr,nrep))
open (unit=1000, file="bandbuilder.in", status='old', action='read') !read feneb.in now that mask is allocated
mask_index = 1
do
    read (1000,"(a)",iostat=ierr) line ! read line into character variable
    if (ierr /= 0) exit
           read (line,*) keyword ! read first keyword of line
    if (keyword == 'mask') then
      line=line(6:) ! remove the first five characters which are 'mask'
      substr = ''
    do i = 1, len_trim(line)+1
        if (line(i:i) == ' ' .or. i==len_trim(line)+1) then
                 ! Here you have completed a "word" and can handle it as a single number
                 ! or the start of a range
          if (index(substr, '-') > 0) then
            call get_range(substr, range) ! parse range of numbers
            read(range(1), *) start_range
            read(range(2), *) end_range
            do i_mask = start_range, end_range
                if (mask_index > nrestr) then
                  print *, "Error: Too many indices for 'mask'. Increase 'nrestr'."
                  stop
                end if
                mask(mask_index) = i_mask
                mask_index = mask_index + 1
            end do
          else
            read(substr, *) i_mask
            if (mask_index > nrestr) then
                print *, "Error: Too many indices for 'mask'. Increase 'nrestr'."
                stop
            end if
            mask(mask_index) = i_mask
            mask_index = mask_index + 1
          end if
          substr = ''  ! Reset substr for the next "word"
      else
          substr = trim(substr) // line(i:i)
      end if
    end do
    end if
    if (keyword == 'all') then
      do i=1,nrestr
        mask(i)=i
      end do
    end if
end do
close (unit=1000)

end subroutine readinputbuilder

subroutine readinputextrator(nrestr,mask,infile,i)

implicit none
character(len=50) :: infile, exp, keyword, line
integer  :: nrestr, i, j, ierr
integer, allocatable, dimension (:) :: mask

open (unit=1000, file="extractor.in", status='old', action='read') !read align.in
do
   read (1000,"(a)",iostat=ierr) line ! read line into character variable
   if (ierr /= 0) exit
   read (line,*) keyword ! read first keyword of line
   if (keyword == 'infile') read(line,*) exp,infile
   if (keyword == 'nrestr') read(line,*) exp, nrestr
   if (keyword == 'rep') read(line,*) exp, i
end do
close (unit=1000)
allocate(mask(nrestr))
open (unit=1000, file="extractor.in", status='old', action='read') !read align.in
do
   read (1000,"(a)",iostat=ierr) line ! read line into character variable
   if (ierr /= 0) exit
   read (line,*) keyword ! read first keyword of line
   if (keyword == 'mask') read(line,*) exp, mask(1:nrestr)
   if (keyword == 'all') then
     do j=1,nrestr
       mask(j)=j
     end do
   end if
end do
close (unit=1000)

end subroutine readinputextrator

subroutine readinputsegmentsfreenergy(nrestr,mask,infile,reffile,i,netcdf,bins,kref,temp,per,velin,boxinfo,truncar,normalize)
implicit none
character(len=50) :: infile, reffile, exp, keyword, line
integer  :: nrestr, i, j, ierr, bins
logical :: netcdf, per, velin, normalize
integer, allocatable, dimension (:) :: mask
double precision :: kref, temp, truncar
double precision, dimension(6) :: boxinfo

netcdf=.false.
temp=300.d0
kref=500.d0
bins=50
truncar=0
open (unit=1000, file="segments.in", status='old', action='read') !read align.in
do
   read (1000,"(a)",iostat=ierr) line ! read line into character variable
   if (ierr /= 0) exit
   read (line,*) keyword ! read first keyword of line
   if (keyword == 'infile') read(line,*) exp, infile
   if (keyword == 'reffile') read(line,*) exp, reffile
   if (keyword == 'nrestr') read(line,*) exp, nrestr
   if (keyword == 'rep') read(line,*) exp, i
   if (keyword == 'netcdf') read(line,*) exp, netcdf
   if (keyword == 'temp') read(line,*) exp, temp
   if (keyword == 'kref') read(line,*) exp, kref
   if (keyword == 'bins') read(line,*) exp, bins
   if (keyword == 'per') read(line,*) exp, per
   if (keyword == 'velin') read(line,*) exp, velin
   if (keyword == 'truncar') read(line,*) exp, truncar
   if (keyword == 'normalize') read(line,*) exp, normalize
end do
close (unit=1000)
allocate(mask(nrestr))
open (unit=1000, file="segments.in", status='old', action='read') !read align.in
do
   read (1000,"(a)",iostat=ierr) line ! read line into character variable
   if (ierr /= 0) exit
   read (line,*) keyword ! read first keyword of line
   if (keyword == 'mask') read(line,*) exp, mask(1:nrestr)
   if (keyword == 'all') then
     do j=1,nrestr
       mask(j)=j
     end do
   end if
end do
close (unit=1000)

end subroutine readinputsegmentsfreenergy

subroutine readhistory(iteration,FIRE_Ndescend,FIRE_dt,FIRE_alpha,FIRE_vel,nrep,nrestr)
implicit none
integer, intent(in) :: nrep, nrestr
integer, intent(out) :: FIRE_Ndescend, iteration
double precision, dimension(nrep), intent(out):: FIRE_dt, FIRE_alpha
double precision, dimension(3,nrestr,nrep), intent(out) :: FIRE_vel
integer :: i ,j

  write(9999,*) "Reading feneb.history"
  open (unit=2200011, file="feneb.history", action='read')
  read(2200011,*) iteration,FIRE_Ndescend
  do i=1,nrep
    read(2200011,*) FIRE_dt(i),FIRE_alpha(i)
  end do
  do i=1,nrestr
    do j=1,nrep
      read(2200011,*) FIRE_vel(1,i,j),FIRE_vel(2,i,j),FIRE_vel(3,i,j)
    end do
  end do
  close(unit=2200011)

end subroutine readhistory

subroutine getfilenames(rep,chrep,infile,reffile,outfile,iname,rname,oname,avname)

implicit none
integer, intent(in) :: rep
character(len=50), intent(in) :: infile, reffile, outfile
character(len=50), intent(out) :: chrep, iname, rname, oname, avname

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
avname = trim(outfile) // "_av_"
avname = trim(avname) // trim(chrep)
avname = trim(avname) // ".rst7"
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

subroutine getcoordfromnetcdf(iname,nsteps,natoms,spatial,coordx,coordy,coordz,coordall,nrestr,mask)
use netcdf
implicit none
real (kind=4), DIMENSION(nsteps) :: coordx, coordy, coordz
integer(kind=4), intent(in) :: natoms, nsteps, spatial, nrestr
integer(kind=4) :: ncid
character(len=50), intent(in) :: iname
integer, dimension(3) :: point,endp
integer, dimension(nrestr) :: mask
double precision, dimension(3,nrestr,nsteps) :: coordall
integer :: i,k,ati

call check(nf90_open(iname, nf90_nowrite, ncid))
do i = 1,nrestr
  ati=mask(i)
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
    coordall(1,i,k)=coordx(k)
    coordall(2,i,k)=coordy(k)
    coordall(3,i,k)=coordz(k)
  end do
enddo
call check(nf90_close(ncid))
end subroutine getcoordfromnetcdf

subroutine getcoordfromfenebtraj(nsteps,coordall,nrestr,i)
implicit none
character(len=50) :: chi
integer :: i, j, k, nsteps, nrestr
double precision, allocatable, dimension(:,:,:), intent(inout) :: coordall

if (i .le. 9) write(chi,'(I1)') i
if (i .gt. 9 .and. i .le. 99) write(chi,'(I2)') i
chi = "feneb.traj."//chi
open (unit=22000, file=chi, action='read')
read(22000,*) nsteps
allocate(coordall(3,nrestr,nsteps))
do j=1,nrestr
  do k=1,nsteps
    read(22000,*) coordall(1,j,k),coordall(2,j,k),coordall(3,j,k)
  end do
end do
close(unit=22000)
end subroutine getcoordfromfenebtraj


subroutine getravfav(coordall,nsteps,natoms,nrestr,mask,kref,rav,devav,fav,nrep,rep,rref,wgrad, &
                            skip,wtemp,dt,mass,tempfilesize,temp)
implicit none
integer, intent(in) :: nrestr,nrep,rep,nsteps,natoms
integer, dimension(nrestr), intent(in) :: mask
double precision :: kref, dt, vat, ekin
double precision, dimension(3) :: av
double precision, dimension(3,nrestr,nrep) :: rav,fav,devav
double precision, dimension(3,nrestr,nsteps) :: coordall
double precision, dimension(3,natoms), intent(in) :: rref
double precision, dimension(tempfilesize,nrep), intent(inout) :: temp
double precision, dimension(nrestr), intent(in) :: mass
character(len=50) :: chi, chrep
integer :: i,j,k,ati,auxunit,skip,tempfilesize,nsteps2
logical :: wgrad,wtemp

ekin=0.d0
temp(1:tempfilesize,rep)=0.d0

do i=1,nrestr
  if (i .le. 9) write(chi,'(I1)') i
  if (i .gt. 9 .and. i .le. 99) write(chi,'(I2)') i
  if (rep .le. 9) write(chrep,'(I1)') rep
  if (rep .gt. 9 .and. rep .le. 99) write(chrep,'(I2)') rep
  chi = "at."//chi
  chi = trim(chi)//".rep."
  chi = trim(chi)//trim(chrep)
  ati=mask(i)
  av=0
  auxunit = 1000*i+rep
  if (wgrad) open(auxunit, file = chi, status = 'unknown')
  do k=skip+1,nsteps
    do j=1,3
      av(j)=(av(j)*(dble(k-1-skip))+coordall(j,i,k))/dble(k-skip)
    end do
    if (wgrad) write(auxunit,*) k-skip, av(1:3)
  end do

  rav(1:3,i,rep)=av(1:3)
  fav(1:3,i,rep)=kref*(av(1:3)-rref(1:3,ati))

  if (wgrad) close(auxunit)
  av=0
  do k=skip+1,nsteps
    do j=1,3
     av(j)=(av(j)*(dble(k-1-skip))+coordall(j,i,k))/dble(k-skip)
     devav(j,i,rep)=devav(j,i,rep)+(av(j)-rav(j,i,rep))**2
    end do
  end do
  do j=1,3
    devav(j,i,rep)=dsqrt(devav(j,i,rep)/(nsteps-skip-1))
  end do


!For comparision with .rst7 file, remember that
!AMBER VEL UNITS are Angstroms per 1/20.455 ps
  if (wtemp) then
    do k=1,nsteps-1
       vat=((coordall(1,i,k+1)-coordall(1,i,k))/dt)**2 + &
           ((coordall(2,i,k+1)-coordall(2,i,k))/dt)**2 + &
           ((coordall(3,i,k+1)-coordall(3,i,k))/dt)**2
       temp(k,rep)=temp(k,rep)+ &
       (mass(i)*vat*10d0/(8.314472d0*3.d0*nrestr))
    end do
  end if

end do
end subroutine getravfav

subroutine getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin)

implicit none
integer :: nrestr
integer, intent(out) :: natoms
character(len=50), intent(in) :: rname
double precision, dimension(:,:), allocatable :: rref
double precision, dimension(6), intent(out) :: boxinfo
integer, dimension(nrestr),intent(in) :: mask
integer :: i
logical ::  per, velin

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
if (mod(natoms,2) .ne. 0) read(1002,'(3(f12.7))') rref(1:3,2*i-1)
if (velin) then
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

subroutine getposforcesextrema(rav,fav,nrestr,nrep)
implicit none
double precision, dimension(3,nrestr,nrep), intent(inout) :: rav, fav
integer, intent(in) :: nrestr,nrep
integer :: i, exp

write(9999,*) "Reading extrema coordinates and forces from feneb.extrema1 and feneb.extrema2 files"

open(unit=210421, file="feneb.extrema1", status='old', action='read')
do i=1,nrestr
  read(210421,'(2x, I6,2x, 6(f20.10,2x))') exp, rav(1:3,exp,1), fav(1:3,exp,1)
end do
close(210421)
open(unit=210422, file="feneb.extrema2", status='old', action='read')
do i=1,nrestr
  read(210422,'(2x, I6,2x, 6(f20.10,2x))') exp, rav(1:3,exp,nrep), fav(1:3,exp,nrep)
end do
close(210422)
end subroutine getposforcesextrema


subroutine check(istatus)
use netcdf
implicit none
integer, intent (in) :: istatus
if (istatus /= nf90_noerr) then
write(9999,*) trim(adjustl(nf90_strerror(istatus)))
end if
end subroutine check

end module readandget
