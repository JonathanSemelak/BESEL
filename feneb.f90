 program feneb
use netcdf
use readandget
implicit none
character(len=50) :: infile, reffile, outfile,topfile, chi, iname, rname, oname, tempname
integer :: nsteps, spatial, natoms, nrestr, nrep, nscycle,maxforceat, atj
integer :: i, j, k, n, start, end, skip, wtempstart, wtempend, wtempfrec, tempfilesize, minsegmentlenght, nevalfluc
integer, allocatable, dimension (:) :: mask
real(4) :: coordinate
real(4), allocatable, dimension (:) :: coordx,coordy,coordz, coordstat
integer, dimension (3) :: point
double precision :: kref, steep_size, ftol, maxforce, kspring, maxforceband, lastmforce, maxforcebandprevsetp, steep_spring
double precision :: stepl, deltaA, rmsfneb, minpoint, maxpoint, barrier, dt, Z, goodrav
double precision, dimension(6) :: boxinfo
double precision, allocatable, dimension(:) :: rmsd, mass
double precision, allocatable, dimension(:,:) :: rref, profile, temp
double precision, allocatable, dimension(:,:,:) :: rav, fav, tang, ftang, ftrue, fperp, rrefall, ravprevsetp
double precision, allocatable, dimension(:,:,:) :: fspring, dontg, selfdist,coordall
logical ::  per, velin, velout, relaxd, converged, wgrad, wtemp, moved, maxpreached, equispaced, rextrema, test
logical ::  dostat, H0, H0T



!------------ Read input
    call readinput(nrep,infile,reffile,outfile,topfile,mask,nrestr,lastmforce, &
                 rav,fav,ftrue,ftang,fperp,fspring,tang,kref,kspring,steep_size,steep_spring, &
                 ftol,per,velin,velout,wgrad,wtemp,dt,wtempstart,wtempend,wtempfrec,mass,rrefall, &
                 nscycle,dontg,ravprevsetp, &
                 rextrema, skip,dostat, minsegmentlenght,nevalfluc)

!------------
 test=.False.

 open(unit=9999, file="feneb.out") !Opten file for feneb output
!------------ Main loop
  if (nrep .eq. 1) then !FE opt only
    write(9999,*) "---------------------------------------------------"
    write(9999,*) "Performing FE full optmization for a single replica"
    write(9999,*) "---------------------------------------------------"

    call getfilenames(nrep,chi,infile,reffile,outfile,iname,rname,oname)
    call getdims(iname,nsteps,spatial,natoms)
    ! nsteps=1000 ! TEST
    tempfilesize=(nsteps-1)
    allocate(temp(tempfilesize,nrep))
    if(wtemp) call readtop(topfile,natoms,mask,mass,nrestr)
    call readtop(topfile,natoms,mask,mass,nrestr)
    if (allocated(coordx)) deallocate(coordx)
    if (allocated(coordy)) deallocate(coordy)
    if (allocated(coordz)) deallocate(coordz)
    if (allocated(rref)) deallocate(rref)

    allocate(coordx(nsteps),coordy(nsteps),coordz(nsteps),rref(3,natoms))
    allocate(coordall(3,nrestr,nsteps),coordstat(nsteps))
    call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin)

    call getavcoordanforces(iname,nsteps,natoms,spatial,coordx,coordy,coordz,&
                        coordall,nrestr,mask,kref,rav,fav,nrep,nrep,rref,wgrad,dontg,&
                        skip,wtemp,dt,mass,tempfilesize,temp)

    if (dostat) then
    write(9999,*) "---------------------------------------------------"
    write(9999,*) "Statistic stuff"
    write(9999,*) "---------------------------------------------------"
    H0T=.True.
    do j=1,nrestr
      atj=mask(j)
      write(9999,*) "Atom:",j
      write(9999,*) "MK test H0"
      coordstat(1:nsteps)=coordall(1,j,1:nsteps)
      call getsstatistics(coordstat,nsteps,skip,nevalfluc,dt,Z,H0,minsegmentlenght,goodrav)
      write(*,*) "ASD2"
      write(9999,*) "coord x:",H0
      rav(1,j,nrep)=goodrav
      fav(1,j,nrep)=kref*(rav(1,j,nrep)-rref(1,atj))
      H0T=(H0T.and.H0)

      coordstat(1:nsteps)=coordall(2,j,1:nsteps)

      call getsstatistics(coordstat,nsteps,skip,nevalfluc,dt,Z,H0,minsegmentlenght,goodrav)
      write(*,*) "ASD2"
      write(9999,*) "coord y:",H0
      rav(2,j,nrep)=goodrav
      fav(2,j,nrep)=kref*(rav(2,j,nrep)-rref(2,atj))
      H0T=(H0T.and.H0)

      coordstat(1:nsteps)=coordall(3,j,1:nsteps)
      call getsstatistics(coordstat,nsteps,skip,nevalfluc,dt,Z,H0,minsegmentlenght,goodrav)
      write(*,*) "ASD2"
      write(9999,*) "coord z:",H0
      rav(3,j,nrep)=goodrav
      fav(3,j,nrep)=kref*(rav(3,j,nrep)-rref(3,atj))
      H0T=(H0T.and.H0)
    end do
    write(*,*) "Trend free:", H0T
    write(9999,*) "---------------------------------------------------"

    end if
    if (wtemp) then
      open(unit=2203280, file="temperature.dat")
      do j=wtempstart, wtempend
        if (mod(j,wtempfrec) .eq. 0) write(2203280,*) j, temp(j,1)
      end do
      close(2203280)
    end if

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
       call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,.True.)
       call writenewcoord(oname,rref,boxinfo,natoms,nrestr,mask,per,velout,rav,nrep,nrep,test)
       write(9999,*) "System converged: F"
    else
       write(9999,*) "System converged: T"
       write(9999,*) "Convergence criteria of ", ftol, " (kcal/mol A) achieved"
    endif

  elseif (nrep .gt. 1) then !NEB on FE surface

    write(9999,*) "---------------------------------------------------"
    write(9999,*) "       Performing NEB on the FE surface"
    write(9999,*) "---------------------------------------------------"

!------------ Set forces to zero

    fav=0.d0

!------------ Band loop
    if (rextrema) then
      start=2
      end=nrep-1
    else
      start=1
      end=nrep
    endif

    do i=start,end
      call getfilenames(i,chi,infile,reffile,outfile,iname,rname,oname)
      call getdims(iname,nsteps,spatial,natoms)
      tempfilesize=(nsteps-1)
      if(i .eq. start) allocate(temp(tempfilesize,nrep))
      if(i .eq. start .and. wtemp) call readtop(topfile,natoms,mask,mass,nrestr)
      if (allocated(coordx)) deallocate(coordx)
      if (allocated(coordy)) deallocate(coordy)
      if (allocated(coordz)) deallocate(coordz)
      if (allocated(rref)) deallocate(rref)
      allocate(coordx(nsteps),coordy(nsteps),coordz(nsteps),rref(3,natoms))

      call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin)
      call getavcoordanforces(iname,nsteps,natoms,spatial,coordx,coordy, coordz,&
                    coordall,nrestr,mask,kref,rav,fav,nrep,i,rref,wgrad,dontg,&
                    skip,wtemp,dt,mass,tempfilesize,temp)
    end do

    if (wtemp) then
      open(unit=2203280, file="temperature.dat")
      do i=start,end
        do j=wtempstart, wtempend
          if (mod(j,wtempfrec) .eq. 0) write(2203280,*) j, temp(j,i)
        end do
        write(2203280,*)
      end do
      close(2203280)
    end if



    if (rextrema) call getposforcesextrema(rav,fav,nrestr,nrep)

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
!----------- Puts reference values in a single array (rrefall). Currently not used.!TESTTTTTTTT
    test=.True.
    do i=start,end

      call getfilenames(i,chi,infile,reffile,outfile,iname,rname,oname) !rname = NAME_r_i.rst7 ; i=replica
      call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin)
      call getcoordextrema(rref,natoms,rrefall,nrestr,nrep,i,mask)

    end do
    test=.False.

!----------- Compute the free energy profile by umbrella integration
    allocate(profile(2,nrep-1))
    call getprofile(rav,fav,nrep,nrestr,profile)
    call gettang(rav,tang,nrestr,nrep)
    call getnebforce(rav,fav,tang,nrestr,nrep,kspring,maxforceband,ftol,relaxd,&
                    ftrue,ftang,fperp,fspring,.true.,dontg)
  ! fav ---> fneb





    call getbarrier(profile, nrep, barrier, minpoint, maxpoint)

    allocate(selfdist(2,nrestr,nrep-1))
    call getselfdist(rav, rrefall, nrep, nrestr, selfdist)

    open(unit=1644, file="selfdist_rref.dat")
    open(unit=1645, file="selfdist_rav.dat")
    do i=1,nrestr
      do n=1,nrep-1
        write(1644,'(2x, I6,2x, f20.10)') n, selfdist(1,i,n)
        write(1645,'(2x, I6,2x, f20.10)') n, selfdist(2,i,n)
      end do
    write(1644,*)
    write(1645,*)
    end do
    close(1644)
    close(1645)


!----------- moves the band
    if (.not. converged) then
       ! if (rpoint .eq. 1) rav=rrefall !TEST USAR RESTR POS
       ! if (rpoint .eq. 2) rav=rcorr
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

        ! maxpreached=.False.
        equispaced=.False.
        k=1
        do while ((k .le. nscycle) .and. (.not. equispaced))
          ! if (.not. maxpreached) then
          !Computes spring force and others
          call gettang(rav,tang,nrestr,nrep) !to test whether to recalculate tg or not

          call getnebforce(rav,fav,tang,nrestr,nrep,kspring,maxforceband,ftol,converged,&
                          ftrue,ftang,fperp,fspring,.false.,dontg)

          !como wrmforce es false, acá usa fspring para determinar maxforceband
          !Moves band using spring force only
          ! write(9999,'(A6,2x,I3,2x,A15,2x,f8.6,2x,A20,2x,f8.6)') "Step: ", k, "Step length", &

          !NO QUIERO TANTO OUTPUT
          !write(9999,'(A6,2x,I7,2x,A15,2x,f20.10,2x,A20,2x,f20.10)') "Step: ", k, "Step length", &
          ! steep_spring, "Band max fspringN: ", maxforceband
          dontg=0.d0
          ravprevsetp=rav
          maxforcebandprevsetp=maxforceband

            do i=2,nrep-1
              call steep(rav,fspring,nrep,i,steep_spring,maxforcebandprevsetp,nrestr,lastmforce,stepl,deltaA,dontg)
            end do

          call getdistrightminusleft(rav, nrep, nrestr, equispaced)
          if ((k .eq. nscycle) .or. equispaced) then
            write(9999,*) "-----------------------------------------------------"
            write(9999,*) "Band max fspringLast: ", maxforceband
            write(9999,*) "Total spring steps: ", k
            write(9999,*) "Equispaced: ", equispaced
            write(9999,*) "-----------------------------------------------------"
          end if
          k=k+1
        end do

        call getselfdist(rav, rrefall, nrep, nrestr, selfdist)

        open(unit=1646, file="selfdist_afterstep.dat")
        do i=1,nrestr
          do n=1,nrep-1
            write(1646,'(2x, I6,2x, f20.10)') n, selfdist(2,i,n)
          end do
          write(1646,*)
        end do

    !------------ Get coordinates for previously optimized extrema
    if (.not. rextrema) then
        !Reactants
        call getfilenames(1,chi,infile,reffile,outfile,iname,rname,oname)
        call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin)
        call getcoordextrema(rref,natoms,rav,nrestr,nrep,1,mask)
        call writenewcoord(oname,rref,boxinfo,natoms,nrestr,mask,per,velout,rav,nrep,1,test)
        !Products
        call getfilenames(nrep,chi,infile,reffile,outfile,iname,rname,oname)
        call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin)
        call getcoordextrema(rref,natoms,rav,nrestr,nrep,nrep,mask)
        call writenewcoord(oname,rref,boxinfo,natoms,nrestr,mask,per,velout,rav,nrep,nrep,test)
    end if
        do i=2,nrep-1
          ! call getfilenames(i,chi,infile,reffile,outfile,iname,rname,oname)
         call getfilenames(i,chi,infile,infile,outfile,iname,rname,oname) !toma ultima foto p/ siguiente paso
         call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,.True.)
         call writenewcoord(oname,rref,boxinfo,natoms,nrestr,mask,per,velout,rav,nrep,i,test)
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
