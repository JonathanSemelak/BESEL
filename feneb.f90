 program feneb
use netcdf
use readandget
implicit none
character(len=50) :: infile, reffile, outfile, chi, iname, rname, oname
integer :: nsteps, spatial, natoms, nrestr, nrep, nscycle,maxforceat
integer :: i, j, k
integer, allocatable, dimension (:) :: mask
real(4) :: coordinate
real(4), allocatable, dimension (:) :: coordx,coordy,coordz
integer, dimension (3) :: point
double precision :: kref, steep_size, ftol, maxforce, kspring, maxforceband, lastmforce, maxforcebandprevsetp, steep_spring
double precision :: stepl, deltaA, rmsfneb, minpoint, maxpoint, barrier
double precision, dimension(6) :: boxinfo
double precision, allocatable, dimension(:) :: rmsd
double precision, allocatable, dimension(:,:) :: rref, profile
double precision, allocatable, dimension(:,:,:) :: rav, fav, tang, ftang, ftrue, fperp, rrefall, ravprevsetp
double precision, allocatable, dimension(:,:,:) :: fspring, dontg
logical ::  per, velin, velout, relaxd, converged, wgrad, moved, maxpreached

!------------ Read input
    call readinput(nrep,infile,reffile,outfile,mask,nrestr,lastmforce, &
                 rav,fav,ftrue,ftang,fperp,fspring,tang,kref,kspring,steep_size, &
                 ftol,per,velin,velout,wgrad,rrefall,nscycle,dontg,ravprevsetp)
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
                        nrestr,mask,kref,rav,fav,nrep,nrep,rref,wgrad,dontg)

    call writeposforces(rav,fav,nrestr,nrep,nrep)

    call getmaxforce(nrestr,nrep,nrep,fav,maxforce,ftol,relaxd,maxforceat,rmsfneb)


    write(9999,*) "Max force: ", maxforce

    if (.not. relaxd) then
       call steep(rav,fav,nrep,nrep,steep_size,maxforce,nrestr,lastmforce,stepl,deltaA,dontg)
       if (stepl .lt. 1d-10) then
         write(9999,*) "-----------------------------------------------------"
         write(9999,*) "Warning: max precision reached on atomic displacement"
         write(9999,*) "step length has been set to zero"
         write(9999,*) "-----------------------------------------------------"
       end if

       call getfilenames(nrep,chi,infile,infile,outfile,iname,rname,oname) !toma ultima foto p/ siguiente paso
       call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin)

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



    !Products
    call getfilenames(nrep,chi,infile,reffile,outfile,iname,rname,oname)
    call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin)
    call getcoordextrema(rref,natoms,rav,nrestr,nrep,nrep,mask)

    !Forces set to zero
    fav=0.d0

!------------ Band loop
    do i=1,nrep

      call getfilenames(i,chi,infile,reffile,outfile,iname,rname,oname)
      call getdims(iname,nsteps,spatial,natoms)

      if (allocated(coordx)) deallocate(coordx)
      if (allocated(coordy)) deallocate(coordy)
      if (allocated(coordz)) deallocate(coordz)
      if (allocated(rref)) deallocate(rref)
      allocate(coordx(nsteps),coordy(nsteps),coordz(nsteps),rref(3,natoms))

      call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin)

      call getavcoordanforces(iname,nsteps,natoms,spatial,coordx,coordy, coordz,&
                    nrestr,mask,kref,rav,fav,nrep,i,rref,wgrad,dontg)
    end do


!----------- Puts reference values in a single array (rrefall). Currently not used.
    do i=1,nrep
      call getfilenames(i,chi,infile,reffile,outfile,iname,rname,oname)
      call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin)
      call getcoordextrema(rref,natoms,rrefall,nrestr,nrep,i,mask)
    end do


!----------- Write mean pos and forces
    do i=1,nrep
      call writeposforces(rav,fav,nrestr,i,nrep)
    end do

!----------- Write RMSD

    allocate(rmsd(nrep))
    call getrmsd(fav, kref, nrep, nrestr,rmsd)
    open(unit=40000, file="rmsd.dat")
      write(40000,*) rmsd(1:nrep)
    close(40000)

!----------- Compute the free energy profile by umbrella integration
    allocate(profile(2,nrep-1))
    call getprofile(rav,fav,nrep,nrestr,profile)
    call getbarrier(profile, nrep, barrier, minpoint, maxpoint)

!----------- Computes tangent and nebforce
    call gettang(rav,tang,nrestr,nrep)

    call getnebforce(rav,fav,tang,nrestr,nrep,kspring,maxforceband,ftol,relaxd,&
                    ftrue,ftang,fperp,fspring,.true.,dontg)
! fav ---> fneb

!----------- moves the band
    if (.not. converged) then
      ! rav=rrefall !TEST USAR RESTR POS

       do i=2,nrep-1
          if (.not. relaxd) call steep(rav,fperp,nrep,i,steep_size,maxforceband,nrestr,lastmforce,stepl,deltaA,dontg)
        end do

        write(9999,'(1x,a,f8.6)') "Step length: ", stepl
        if (stepl .lt. 1d-5) then
          write(9999,*) "-----------------------------------------------------"
          write(9999,*) "Warning: max precision reached on atomic displacement"
          write(9999,*) "step length has been set to zero"
          write(9999,*) "-----------------------------------------------------"
        end if

        rmsfneb=0.d0
        do i=1,nrep
          call getmaxforce(nrestr,nrep,i,fav,maxforce,ftol,relaxd,maxforceat,rmsfneb)
        end do
        rmsfneb=dsqrt(rmsfneb/dble(nrep*nrestr))

        write(9999,'(1x,a,f8.6)') "rmsfneb(FNEB): ", rmsfneb/nrep

        if (nscycle .eq. 1) write(9999,*) "WARNING: Using only fperp to move the band!"
        if (nscycle .gt. 1) then

          write(9999,*) "-----------------------------------------------------"
          write(9999,*) "Performing extra optimization steps using fspring    "
          write(9999,*) "to get a better distribution of replicas.            "
          write(9999,'(1x,a,I4)') "Extra optmization movements: ", nscycle
          write(9999,*) "-----------------------------------------------------"
        end if

        steep_spring=steep_size
        maxpreached=.False.
        do k=1,nscycle
          if (.not. maxpreached) then
          !Computes spring force and others
          call getnebforce(rav,fav,tang,nrestr,nrep,kspring,maxforceband,ftol,converged,&
                          ftrue,ftang,fperp,fspring,.false.,dontg)

          !como wrmforce es false, acá usa fspring para determinar maxforceband
          !Moves band using spring force only
          write(9999,'(A6,2x,I3,2x,A15,2x,f8.6,2x,A20,2x,f8.6)') "Step: ", k, "Step length", &
          steep_spring, "Band max fspringN: ", maxforceband
          ! write(*,*) k,maxforceband
          dontg=0.d0
          ravprevsetp=rav
          maxforcebandprevsetp=maxforceband
          moved=.False.
          do while (.not. moved)
            do i=2,nrep-1
              call steep(rav,fspring,nrep,i,steep_spring,maxforcebandprevsetp,nrestr,lastmforce,stepl,deltaA,dontg)
            end do
            call getnebforce(rav,fav,tang,nrestr,nrep,kspring,maxforceband,ftol,converged,&
                            ftrue,ftang,fperp,fspring,.false.,dontg)
            if ((maxforceband .lt. maxforcebandprevsetp) .or.  (stepl .lt. 1d-10)) then
               moved=.True.
            else
               steep_spring=steep_spring*0.5d0
               rav=ravprevsetp
            end if
          end do
          if (stepl .lt. 1d-10) then
            write(9999,*) "Max precision reached"
            maxpreached=.True.
            ! write(9999,*) "Band max fspringLast: ", maxforceband
          end if
          end if
          if (k .eq. nscycle) then
            write(9999,*) "-----------------------------------------------------"
            write(9999,*) "Band max fspringLast: ", maxforceband
            write(9999,*) "-----------------------------------------------------"
          end if
        end do


    !------------ Get coordinates for previously optimized extrema
        !Reactants
        call getfilenames(1,chi,infile,reffile,outfile,iname,rname,oname)
        call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin)
        call getcoordextrema(rref,natoms,rav,nrestr,nrep,1,mask)

        !Products
        call getfilenames(nrep,chi,infile,reffile,outfile,iname,rname,oname)
        call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin)
        call getcoordextrema(rref,natoms,rav,nrestr,nrep,nrep,mask)


        do i=1,nrep
          ! call getfilenames(i,chi,infile,reffile,outfile,iname,rname,oname)
          call getfilenames(i,chi,infile,infile,outfile,iname,rname,oname) !toma ultima foto p/ siguiente paso
          call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin)
          call writenewcoord(oname,rref,boxinfo,natoms,nrestr,mask,per,velout,rav,nrep,i)
        end do
    end if !converged
  write(9999,*) "-----------------------------------------------------"
  write(9999,*) "                   Extrema info"
  write(9999,*) "Barrier: ", barrier
  write(9999,*) "Minimum point: ", minpoint
  write(9999,*) "Maximum point: ", maxpoint
  write(9999,*) "-----------------------------------------------------"

  end if !nrep gt 1

  close(unit=9999)

end program feneb
