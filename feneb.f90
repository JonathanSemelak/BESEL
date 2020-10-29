 program feneb
use netcdf
use readandget
implicit none
character(len=50) :: infile, reffile, outfile, chi, iname, rname, oname
integer :: nsteps, spatial, natoms, nrestr, nrep, nscycle
integer :: i, j, k
integer, allocatable, dimension (:) :: mask
real(4) :: coordinate
real(4), allocatable, dimension (:) :: coordx,coordy,coordz
integer, dimension (3) :: point
double precision :: kref, steep_size, ftol, maxforce, kspring, maxforceband, lastmforce
double precision :: stepl, deltaA
double precision, dimension(6) :: boxinfo
double precision, allocatable, dimension(:,:) :: rref
double precision, allocatable, dimension(:,:,:) :: rav, fav, tang, ftang, ftrue, fperp, rrefall
double precision, allocatable, dimension(:,:,:) :: fspring
logical ::  per, velin, velout, relaxd, converged, wgrad

!------------ Read input
    call readinput(nrep,infile,reffile,outfile,mask,nrestr,lastmforce, &
                 rav,fav,ftrue,ftang,fperp,fspring,tang,kref,kspring,steep_size, &
                 ftol,per,velin,velout,wgrad,rrefall,nscycle)
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
                        nrestr,mask,kref,rav,fav,nrep,nrep,rref,wgrad)

    call writeposforces(rav,fav,nrestr,nrep,nrep)

    call getmaxforce(nrestr,nrep,nrep,fav,maxforce,ftol,relaxd)

    write(9999,*) "Max force: ", maxforce

    ! stepl=0.05d0
    ! DeltaA=-0.15d0
    ! write(9999,'(1x,a,f20.16)') "Step length: ", stepl, "DeltaA: ", deltaA
    ! write(9999,*) "System converged: F"
    ! call writenewcoord(oname,rref,boxinfo,natoms,nrestr,mask,per,velout,rav,nrep,nrep)


    if (.not. relaxd) then
       call steep(rav,fav,nrep,nrep,steep_size,maxforce,nrestr,lastmforce,stepl,deltaA)
       write(9999,'(1x,a,f20.16)') "Step length: ", stepl, "DeltaA: ", deltaA

       if (stepl .lt. 1d-10) then
         write(9999,*) "-----------------------------------------------------"
         write(9999,*) "Warning: max precision reached on atomic displacement"
         write(9999,*) "step length has been set to zero"
         write(9999,*) "-----------------------------------------------------"
       end if
       call writenewcoord(oname,rref,boxinfo,natoms,nrestr,mask,per,velout,rav,nrep,nrep)
       write(9999,*) "System converged: F"
    else
       write(9999,*) "System converged: T"
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

    ! write(*,*) "rc"
    ! do i=1,nrestr
    ! write(*,'(2x, I6,2x, 6(f20.10,2x))') i, rav(1:3,i,1)
    ! end do

    !Products
    call getfilenames(nrep,chi,infile,reffile,outfile,iname,rname,oname)
    call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin)
    call getcoordextrema(rref,natoms,rav,nrestr,nrep,nrep,mask)

    ! write(*,*) "rc"
    ! do i=1,nrestr
    ! write(*,'(2x, I6,2x, 6(f20.10,2x))') i, rav(1:3,i,nrep)
    ! end do
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
                    nrestr,mask,kref,rav,fav,nrep,i,rref,wgrad)
    end do


!----------- Computes
    do i=1,nrep
      call getfilenames(i,chi,infile,reffile,outfile,iname,rname,oname)
      call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin)
      call getcoordextrema(rref,natoms,rrefall,nrestr,nrep,i,mask)
    end do

!----------- Computes tangent and nebforce
    !call gettang(rrefall,tang,nrestr,nrep)

    call gettang(rav,tang,nrestr,nrep)

    call getnebforce(rav,fav,tang,nrestr,nrep,kspring,maxforceband,ftol,converged,&
                    ftrue,ftang,fperp,fspring,.true.)
! fav ---> fneb
!----------- Write mean pos and forces
    do i=1,nrep
!      call writeposforces(rav,fav,nrestr,i,nrep)
      call writeposforces(rav,ftang,nrestr,i,nrep)
    end do


!----------- moves the band
    ! converged = .TRUE.
    if (.not. converged) then

      ! if (nscycle .eq. 0) then

        do i=1,nrep
          call getmaxforce(nrestr,nrep,i,fperp,maxforce,ftol,relaxd)
          ! write(9999,*) "Replica: ", i, "Converged: " relaxd
          if (.not. relaxd) call steep(rav,fperp,nrep,i,steep_size,maxforceband,nrestr,lastmforce,stepl,deltaA)
          ! call getfilenames(i,chi,infile,reffile,outfile,iname,rname,oname)
          ! call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin)
          ! call writenewcoord(oname,rref,boxinfo,natoms,nrestr,mask,per,velout,rav,nrep,i)
      ! converged = (converged .and. relaxd)
        end do
        write(9999,'(1x,a,f8.6)') "Step length: ", stepl
        if (stepl .lt. 1d-5) then
          write(9999,*) "-----------------------------------------------------"
          write(9999,*) "Warning: max precision reached on atomic displacement"
          write(9999,*) "step length has been set to zero"
          write(9999,*) "-----------------------------------------------------"
        end if

      if (nscycle .gt. 0) then

        write(9999,*) "-----------------------------------------------------"
        write(9999,*) "Performing extra optimization steps using fspring    "
        write(9999,*) "get a better distribution of replicas.               "
        write(9999,'(1x,a,I4)') "Extra optmization movements: ", nscycle
        write(9999,*) "-----------------------------------------------------"
      end if
      !   do i=1,nrep
      !     call getmaxforce(nrestr,nrep,i,fperp,maxforce,ftol,relaxd)
      !     ! write(9999,*) "Replica: ", i, "Converged: " relaxd
      !     if (.not. relaxd) call steep(rav,fperp,nrep,i,steep_size,maxforceband,nrestr,lastmforce,stepl,deltaA)
      !     call getfilenames(i,chi,infile,reffile,outfile,iname,rname,oname)
      !     call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin)
      !     call writenewcoord(oname,rref,boxinfo,natoms,nrestr,mask,per,velout,rav,nrep,i)
      ! ! converged = (converged .and. relaxd)
      !   end do

        do k=1,nscycle
          !Computes new tangent
          !call gettang(rav,tang,nrestr,nrep)
          !Computes spring force and others
          call getnebforce(rav,fav,tang,nrestr,nrep,kspring,maxforceband,ftol,converged,&
                          ftrue,ftang,fperp,fspring,.false.)
          !Moves band using spring force only
          do i=2,nrep-1
            !call getmaxforce(nrestr,nrep,i,fperp,maxforce,ftol,relaxd)
            !write(*,*) i, relaxd
            call steep(rav,fspring,nrep,i,steep_size,maxforceband,nrestr,lastmforce,stepl,deltaA)
          end do
        end do

        do i=1,nrep
          call getfilenames(i,chi,infile,reffile,outfile,iname,rname,oname)
          call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin)
          call writenewcoord(oname,rref,boxinfo,natoms,nrestr,mask,per,velout,rav,nrep,i)
      ! converged = (converged .and. relaxd)
        end do

      ! end if
    end if

  end if


  close(unit=9999)

end program feneb
