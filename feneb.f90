program feneb
use netcdf
use readandget
implicit none
character(len=50) :: infile, reffile, outfile, chi, iname, rname, oname
integer :: nsteps, spatial, natoms, nrestr, nrep
integer :: i, j
integer, allocatable, dimension (:) :: mask
real(4) :: coordinate
real(4), allocatable, dimension (:) :: coordx,coordy,coordz
integer, dimension (3) :: point
double precision :: kref, steep_size, ftol, maxforce, kspring, maxforceband, lastmforce
double precision :: stepl
double precision, dimension(6) :: boxinfo
double precision, allocatable, dimension(:,:) :: rref
double precision, allocatable, dimension(:,:,:) :: rav, fav, tang, ftang
logical ::  per, velin, velout, relaxd, converged

!------------ Read input
    call readinput(nrep,infile,reffile,outfile,mask,nrestr,lastmforce, &
                 rav,fav,ftang,tang,kref,kspring,steep_size,ftol,per,velin,velout)
!------------


 open(unit=9999, file="feneb.out") !Opten file for feneb output
!------------ Main loop
  if (nrep .eq. 1) then !FE opt only
    write(9999,*) "---------------------------------------------------"
    write(9999,*) "Performing FE full optmization for a single replica"
    write(9999,*) "---------------------------------------------------"

    call getfilenames(nrep,chi,infile,reffile,outfile,iname,rname,oname)
    call getdims(iname,nsteps,spatial,natoms)

    if (allocated(coordx)) deallocate(coordx)
    if (allocated(coordy)) deallocate(coordy)
    if (allocated(coordz)) deallocate(coordz)
    if (allocated(rref)) deallocate(rref)

    allocate(coordx(nsteps),coordy(nsteps),coordz(nsteps),rref(3,natoms))

    call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin)
    call getavcoordanforces(iname,nsteps,natoms,spatial,coordx,coordy,coordz,&
                        nrestr,mask,kref,rav,fav,nrep,nrep,rref)
    call writeposforces(rav,fav,nrestr,nrep)
    call getmaxforce(nrestr,nrep,nrep,fav,maxforce,ftol,relaxd)
    write(9999,*) "Max force: ", maxforce
    if (.not. relaxd) then
       call steep(rav,fav,nrep,nrep,steep_size,maxforce,nrestr,lastmforce,stepl)
       call writenewcoord(oname,rref,boxinfo,natoms,nrestr,mask,per,velout,rav,nrep,nrep)
    else
       write(9999,*) "Convergence criteria of ", ftol, " (kcal/mol A) achieved"
    endif

  elseif (nrep .gt. 1) then !NEB on FE surface

    write(9999,*) "---------------------------------------------------"
    write(9999,*) "       Performing NEB on the FE surface"
    write(9999,*) "---------------------------------------------------"
!------------ Get coordinates for previously optimized extrema
!------------ And set forces to zero
    !Reactants
    call getfilenames(1,chi,infile,reffile,outfile,iname,rname,oname)
    call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin)
    call getcoordextrema(rref,natoms,rav,nrestr,nrep,1,mask)
    !Products
    call getfilenames(nrep,chi,infile,reffile,outfile,iname,rname,oname)
    call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin)
    call getcoordextrema(rref,natoms,rav,nrestr,nrep,nrep,mask)

    !Forces set to zero
    fav=0.d0

!------------ Band loop

    do i=2,nrep-1
      call getfilenames(i,chi,infile,reffile,outfile,iname,rname,oname)
      call getdims(iname,nsteps,spatial,natoms)

      if (allocated(coordx)) deallocate(coordx)
      if (allocated(coordy)) deallocate(coordy)
      if (allocated(coordz)) deallocate(coordz)
      if (allocated(rref)) deallocate(rref)
      allocate(coordx(nsteps),coordy(nsteps),coordz(nsteps),rref(3,natoms))

      call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin)

      call getavcoordanforces(iname,nsteps,natoms,spatial,coordx,coordy, coordz,&
                    nrestr,mask,kref,rav,fav,nrep,i,rref)
    end do

!----------- Compute tangent and nebforce

    call gettang(rav,tang,nrestr,nrep)
    call getnebforce(rav,fav,tang,nrestr,nrep,kspring,maxforceband,ftol,converged,ftang)

!----------- Write mean pos and forces
    do i=1,nrep
      call writeposforces(rav,ftang,nrestr,i,nrep)
    end do


!----------- moves the band
    ! converged = .TRUE.
    if (.not. converged) then
      do i=1,nrep
        call getmaxforce(nrestr,nrep,i,fav,maxforce,ftol,relaxd)
        ! write(9999,*) "Replica: ", i, "Converged: " relaxd
        if (.not. relaxd) call steep(rav,fav,nrep,i,steep_size,maxforceband,nrestr,lastmforce,stepl)
        call getfilenames(i,chi,infile,reffile,outfile,iname,rname,oname)
        call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin)
        call writenewcoord(oname,rref,boxinfo,natoms,nrestr,mask,per,velout,rav,nrep,i)
      ! converged = (converged .and. relaxd)
      end do
      write(9999,'(1x,a,f8.6)') "Step length: ", stepl
    end if
  end if


  close(unit=9999)

end program feneb
