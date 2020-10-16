program getdihe
  implicit none
  double precision, dimension (:,:,:), allocatable :: rav
  double precision, dimension (:), allocatable :: phi, psi
  integer, dimension (4) :: mask1, mask2
  double precision :: dihedro2
  integer :: i,j,k, nrestr, npict, trash, at1, at2, at3, at4, ierr
  integer :: nlines, io, atnum, atnum_old
  character(len=50) :: exp, line, keyword


  open (unit=1000, file='getdihe.in', status='old', action='read') !read align.in
  do
     read (1000,"(a)",iostat=ierr) line ! read line into character variable
     if (ierr /= 0) exit
     read (line,*) keyword ! read first keyword of line
     if (keyword == 'phi') read(line,*) exp,mask1(1:4)
     if (keyword == 'psi') read(line,*) exp,mask2(1:4)
  end do
  close (unit=1000)

  !determine npict and nrestr
  open(unit=25, file="Pos_forces.dat")
  nlines=0
  atnum=0
  atnum_old=0
  nrestr=0
  npict=0
  DO
    READ(25,*,iostat=io) atnum
    if (atnum .eq. atnum_old+1) then
      nrestr=nrestr+1
      atnum_old=atnum
    end if
    IF (io/=0) EXIT
    nlines = nlines + 1
  END DO
  npict=nlines/nrestr
  close(25)

  allocate (rav(3,nrestr,npict), phi(npict), psi(npict))

  open(unit=25, file="Pos_forces.dat")
  do i=1,npict
  do j=1,nrestr
    read(25,*) trash, rav(1:3,j,i)
  end do
    read(25,*)
  end do
  close(25)

  do i=1,npict
    at1=mask1(1)
    at2=mask1(2)
    at3=mask1(3)
    at4=mask1(4)

    phi(i)=dihedro2(rav(1,at1,i),rav(2,at1,i),rav(3,at1,i), &
                rav(1,at2,i),rav(2,at2,i),rav(3,at2,i), &
                rav(1,at3,i),rav(2,at3,i),rav(3,at3,i), &
                rav(1,at4,i),rav(2,at4,i),rav(3,at4,i))

    at1=mask2(1)
    at2=mask2(2)
    at3=mask2(3)
    at4=mask2(4)

    psi(i)=dihedro2(rav(1,at1,i),rav(2,at1,i),rav(3,at1,i), &
                rav(1,at2,i),rav(2,at2,i),rav(3,at2,i), &
                rav(1,at3,i),rav(2,at3,i),rav(3,at3,i), &
                rav(1,at4,i),rav(2,at4,i),rav(3,at4,i))
  end do

  open(unit=1001, file="phipsi.dat", position='append')
  do i=1,npict
  write(1001,'(2x,2(f20.10,2x))') phi(i), psi(i)
  end do
  close(1001)

end program getdihe

 function dihedro2(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4)
        implicit none
        double precision dihedro2,l1,l2,pi
        double precision x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
        double precision l,A,B,C,D,m,n,nx,ny,nz,mx,my,mz,scalar

! c       vectors n and m
! c       m = vectorial product 21 and 23
! c       n = vectorial product 23 and 34

        pi=DACOS(-1.d0)
        mx = (y1-y2)*(z3-z2) - (z1-z2)*(y3-y2)
        my = -((x1-x2)*(z3-z2) -(z1-z2)*(x3-x2))
        mz = (x1-x2)*(y3-y2)-(y1-y2)*(x3-x2)

        nx = (y2-y3)*(z4-z3) - (z2-z3)*(y4-y3)
        ny = -((x2-x3)*(z4-z3) -(z2-z3)*(x4-x3))
        nz = (x2-x3)*(y4-y3)-(y2-y3)*(x4-x3)

! c       scalar product n*m

        scalar = mx*nx + my*ny + mz*nz
        m = mx**2d0 + my**2d0 + mz**2d0
        m = m**(0.5d0)
        n = nx**2d0 + ny**2d0 + nz**2d0
        n = n**(0.5d0)

        if (m*n.eq.0.d0) STOP "problem in bias dihedral definition"
        dihedro2 = ACOS(scalar/(m*n))
        dihedro2 = dihedro2*180d0/pi

! c       which cuadrant?
! c       plane generation with points 1 2 3

        A = (y2-y1)*(z3-z1)-(z2-z1)*(y3-y1)
        B = (z2-z1)*(x3-x1)-(x2-x1)*(z3-z1)
        C = (x2-x1)*(y3-y1)-(y2-y1)*(x3-x1)
        D = -x1*A-y1*B-z1*C

! c       distance(l) at4 to ABCD=0 plane

        l1 = A*x4+B*y4+C*z4+D
        l2 = (A**2d0+B**2d0+C**2d0)
        l2 = l2**0.5d0
        l= l1/l2

! c       if l>0 -> dihe<0 , else >0

        if(l.lt.0) then
                ! dihedro2 = 360d0-dihedro2
                dihedro2 = -dihedro2
        endif

        end function dihedro2
