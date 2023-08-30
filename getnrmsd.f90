program getnrmsd
  use netcdf
  use readandget
  implicit none
  double precision, dimension (:,:,:), allocatable :: rav
  double precision, dimension (:,:), allocatable :: rref
  double precision, dimension (:), allocatable :: nrmsd
  double precision, dimension(6) :: boxinfo
  integer, allocatable, dimension (:) :: mask
  double precision :: dihedro2
  integer :: i,j,k, rep, nrep, trash, at1, at2, at3, at4, ierr
  integer :: nlines, io, atnum, atnum_old, nrestr, natoms
  character(len=50) :: exp, line, keyword, chrep, rname,reffile
  logical :: readfromrst7, velin, per



  readfromrst7 = .False.
  open (unit=1000, file='getnrmsd.in', status='old', action='read') !read align.in
  do
     read (1000,"(a)",iostat=ierr) line ! read line into character variable
     if (ierr /= 0) exit
     read (line,*) keyword ! read first keyword of line
     if (keyword == 'readfromrst7') read(line,*) exp,readfromrst7
  end do
  close (unit=1000)
  if (readfromrst7) then
    open (unit=1000, file='getnrmsd.in', status='old', action='read') !read align.in
    do
       read (1000,"(a)",iostat=ierr) line ! read line into character variable
       if (ierr /= 0) exit
       read (line,*) keyword ! read first keyword of line
       if (keyword == 'nrestr') read(line,*) exp,nrestr
       if (keyword == 'nrep') read(line,*) exp,nrep
       if (keyword == 'reffile') read(line,*) exp,reffile
       if (keyword == 'velin') read(line,*) exp,velin
       if (keyword == 'per') read(line,*) exp,per
    end do
    close (unit=1000)
    allocate(mask(nrestr),rav(3,nrestr,nrep), nrmsd(nrep))
    open (unit=1000, file='getnrmsd.in', status='old', action='read') !read align.in
    do
       read (1000,"(a)",iostat=ierr) line ! read line into character variable
       if (ierr /= 0) exit
       read (line,*) keyword ! read first keyword of line
       if (keyword == 'mask') read(line,*) exp,mask(1:nrestr)
       if (keyword == 'all') then
         do i=1,nrestr
           mask(i)=i
         end do
      end if
    end do
    close (unit=1000)
  endif
  !determine nrep and nrestr
  if (.not. readfromrst7) then
    open(unit=25, file="feneb.gradients")
    nlines=0
    atnum=0
    atnum_old=0
    nrestr=0
    nrep=0
    DO
      READ(25,*,iostat=io) atnum
      if (atnum .eq. atnum_old+1) then
        nrestr=nrestr+1
        atnum_old=atnum
      end if
      IF (io/=0) EXIT
      nlines = nlines + 1
    END DO
    nrep=nlines/nrestr
    close(25)

    allocate (rav(3,nrestr,nrep))
    allocate (nrmsd(nrep))

    open(unit=25, file="feneb.gradients")
    do i=1,nrep
    do j=1,nrestr
      read(25,*) trash, rav(1:3,j,i)
    end do
      read(25,*)
    end do
    close(25)
  else
    do rep=1,nrep
      if (rep .le. 9) write(chrep,'(I1)') rep
      if (rep .gt. 9 .and. rep .le. 99) write(chrep,'(I2)') rep
      if (rep .gt. 99 .and. rep .le. 999) write(chrep,'(I3)') rep
      rname = trim(reffile) // "_"
      rname = trim(rname) // trim(chrep)
      rname = trim(rname) // ".rst7"
      call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin)
      call getcoordextrema(rref,natoms,rav,nrestr,nrep,rep,mask)
    end do
  nrep=nrep
  end if

!coordinates in rav
  nrmsd=0.d0
  do i=2,nrep
    do j=1,nrestr
       nrmsd(i)=(rav(1,j,i)-rav(1,j,1))**2 + &
                (rav(2,j,i)-rav(2,j,1))**2 + &
                (rav(3,j,i)-rav(3,j,1))**2
    end do
    nrmsd(i)=dsqrt(nrmsd(i))
  end do
  nrmsd=nrmsd/nrmsd(nrep)

  open(unit=1001, file="nrms.dat", position='append')
  do i=1,nrep
  write(1001,'(2x,2(f20.10,2x))') nrmsd(i)
  end do
  close(1001)

end program getnrmsd
