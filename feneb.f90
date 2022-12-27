!-----------------------------------------------------------------------------------------------!
!  This code performes nudged elastic bands on the free energy surface, using free energy 
!  gradients computed from biased molecular dynamics
!
!  References:
!   - Semelak, et. al. (2022).ChemRxiv. Cambridge: Cambridge Open Engage. 
!      This content is a preprint and has not been peer-reviewed.
!   - Bohner, et. al. (2014). Nudged-elastic band used to find reaction
!     coordinates based on the free energy. The Journal of Chemical Physics, 140(7), 074109.
!  
! Written by J. A. Semelak
!----------------------------------------------------------------------------------------------!

program feneb
use netcdf
use readandget
implicit none
character(len=50) :: infile, reffile, outfile,topfile, chi, iname, rname, oname, avname, tempname
integer :: nsteps, spatial, natoms, nrestr, nrep, nscycle,maxforceat, atj, maxstdat, tangoption, optoption
integer :: i, j, k, n, start, nend, skip, wtempstart, wtempend, wtempfrec, tempfilesize, minsegmentlenght, nevalfluc, nstepsexternal
integer, allocatable, dimension (:) :: mask
real(4) :: coordinate
real(4), allocatable, dimension (:) :: coordx,coordy,coordz, coordstat
integer, dimension (3) :: point
double precision :: kref, steep_size, ftol, maxforce, kspring, maxforceband, lastmforce, maxforcebandprevsetp, steep_spring
double precision :: stepl, rmsfneb, minpoint, maxpoint, barrier, dt, Z, goodrav, gooddevav, maxstd, maxdisp
double precision :: FIRE_dt_max, maxdist
integer :: FIRE_Ndescend, iteration
double precision, dimension(6) :: boxinfo
double precision, allocatable, dimension(:) :: rmsd, mass, FIRE_alpha, FIRE_dt
double precision, allocatable, dimension(:,:) :: rref, profile, temp
double precision, allocatable, dimension(:,:,:) :: rav, fav, tang, ftang, ftrue, fperp, rrefall, ravprevsetp, devav, ravout
double precision, allocatable, dimension(:,:,:) :: fspring, dontg, selfdist,coordall, FIRE_vel
logical ::  per, velin, velout, relaxd, converged, wgrad, wtemp, moved, maxpreached, equispaced, rextrema, test
logical ::  dostat, H0, H0T, rfromtraj, usensteps, smartstep, typicalneb, historyfound, tangrecalc

!------------ Read input
    call readinput(nrep,infile,reffile,outfile,topfile,mask,nrestr,lastmforce, &
                 rav,ravout,devav,fav,ftrue,ftang,fperp,fspring,tang,kref,kspring,steep_size,steep_spring, &
                 ftol,per,velin,velout,wgrad,wtemp,dt,wtempstart,wtempend,wtempfrec,mass,rrefall, &
                 nscycle,dontg,ravprevsetp,rextrema, skip,dostat, minsegmentlenght,nevalfluc,rfromtraj, &
                 usensteps,nstepsexternal,smartstep,typicalneb,tangoption,optoption,FIRE_dt_max, tangrecalc, maxdist)
!test
!------------ Read feneb.history if the optimizer is FIRE
if (optoption.eq.1) then
  allocate(FIRE_vel(3,nrestr,nrep),FIRE_dt(nrep),FIRE_alpha(nrep))
  INQUIRE(FILE="feneb.history", EXIST=historyfound)
  if(historyfound) then
    call readhistory(iteration,FIRE_Ndescend,FIRE_dt,FIRE_alpha,FIRE_vel,nrep,nrestr)
  else
    FIRE_Ndescend=0
    iteration=1
    FIRE_vel=0.d0
    FIRE_dt=0.1d0
  endif
end if
!------------

 test=.False.
 open(unit=9999, file="feneb.out") !Opten file for feneb output
!------------ Main loop
  if (nrep .eq. 1) then !FE opt only
    write(9999,*)
    write(9999,*) "Performing FE full optmization for a single replica"
    write(9999,*)

    call getfilenames(nrep,chi,infile,reffile,outfile,iname,rname,oname,avname)
    call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin) !saves all coordinates from rname (rst7 file) in rref
    call getcoordextrema(rref,natoms,rrefall,nrestr,nrep,nrep,mask) !puts mask coordinates from rref in rrefall selected replica
    if (rfromtraj) then
      if (usensteps) write(9999,*) " The option usensteps is not available if rfromtraj is True "
      if (usensteps) write(9999,*) " Reading nsteps from feneb.traj file "
      if (allocated(coordall)) deallocate(coordall,coordstat)
      call getcoordfromfenebtraj(nsteps,coordall,nrestr,nrep)
      allocate(coordstat(nsteps))
    else
      call getdims(iname,nsteps,spatial,natoms)
      if (usensteps) write(9999,*) "Using ", nstepsexternal, "out of ", nsteps
      if (usensteps) nsteps=nstepsexternal
      if (allocated(coordx)) deallocate(coordx)
      if (allocated(coordy)) deallocate(coordy)
      if (allocated(coordz)) deallocate(coordz)
      allocate(coordx(nsteps),coordy(nsteps),coordz(nsteps))
      allocate(coordall(3,nrestr,nsteps),coordstat(nsteps))
      call getcoordfromnetcdf(iname,nsteps,natoms,spatial,coordx,coordy,coordz,coordall,nrestr,mask)
    end if

    tempfilesize=(nsteps-1)
    allocate(temp(tempfilesize,nrep))
    call readtop(topfile,natoms,mask,mass,nrestr)
    write(9999,*) "Reading masses from file: ", trim(topfile)
    j=1
    do while (j .le. (nrestr/3)*3)
      write(9999,*) mass(j), mass(j+1), mass(j+2)
      j = j + 3
    enddo
    write(9999,*) mass(j:nrestr)
    write(9999,*)

    call getravfav(coordall,nsteps,natoms,nrestr,mask,kref,rav,devav,fav,nrep,nrep,rref,wgrad,skip,wtemp,dt,mass,tempfilesize,temp)

    if (dostat) then
      write(9999,*)
      write(9999,*) "Statistic stuff"
      write(9999,*)
      H0T=.True.
      do j=1,nrestr
        atj=mask(j)
        write(9999,*) "Atom:",j
        write(9999,*) "MK test H0"

        coordstat(1:nsteps)=coordall(1,j,1:nsteps)
        call getsstatistics(coordstat,nsteps,skip,nevalfluc,dt,Z,H0,minsegmentlenght,goodrav,gooddevav)
        write(9999,*) "coord x:",H0
        rav(1,j,nrep)=goodrav
        fav(1,j,nrep)=kref*(rav(1,j,nrep)-rref(1,atj))
        devav(1,j,nrep)=gooddevav
        H0T=(H0T.and.H0)

        coordstat(1:nsteps)=coordall(2,j,1:nsteps)
        call getsstatistics(coordstat,nsteps,skip,nevalfluc,dt,Z,H0,minsegmentlenght,goodrav,gooddevav)
        write(9999,*) "coord y:",H0
        rav(2,j,nrep)=goodrav
        fav(2,j,nrep)=kref*(rav(2,j,nrep)-rref(2,atj))
        devav(2,j,nrep)=gooddevav
        H0T=(H0T.and.H0)

        coordstat(1:nsteps)=coordall(3,j,1:nsteps)
        call getsstatistics(coordstat,nsteps,skip,nevalfluc,dt,Z,H0,minsegmentlenght,goodrav,gooddevav)
        write(9999,*) "coord z:",H0
        rav(3,j,nrep)=goodrav
        fav(3,j,nrep)=kref*(rav(3,j,nrep)-rref(3,atj))
        devav(3,j,nrep)=gooddevav
        H0T=(H0T.and.H0)
      end do
      write(*,*) "Trend free:", H0T
      write(9999,*)
    end if ! dostat

    if (wtemp) then
      open(unit=2203280, file="temperature.dat")
      do j=wtempstart, wtempend
        if (mod(j,wtempfrec) .eq. 0) write(2203280,*) j, temp(j,1)
      end do
      close(2203280)
    end if
    call writeposforces(rav,fav,nrestr,nrep,nrep)
    call writeposdev(rav,devav,nrestr,nrep,nrep)

    call getmaxforce(nrestr,nrep,nrep,fav,maxforce,ftol,relaxd,maxforceat,rmsfneb)
    call getmaxstd(nrestr,nrep,nrep,fav,devav,maxstd,maxstdat,.False.)
    call getmaxdisplacement(nrestr,nrep,rav,rrefall,maxdisp)

    write(9999,*) "Max force: ", maxforce
    write(9999,*) "Max STD: ", maxstd
    write(9999,*) "Max displacement due MD: ", maxdisp

    ravout=rav
    if (.not. relaxd) then
      if(optoption.eq.0) then
       call steep(rav,fav,nrep,nrep,steep_size,maxforce,nrestr,lastmforce,stepl,smartstep)
       if (smartstep) then
         write(9999,*) "Using smartstep option"
         write(9999,*) "Base step length: ", stepl
       endif
       if (stepl .lt. 1d-10) then
         write(9999,*)
         write(9999,*) "Warning: max precision reached on atomic displacement"
         write(9999,*) "step length has been set to zero"
         write(9999,*)
       end if
       elseif(optoption.eq.1)then
         call FIRE(nrestr, rav(:,:,nrep), fav(:,:,nrep), FIRE_vel(:,:,nrep), FIRE_dt(nrep),&
         FIRE_Ndescend, FIRE_dt_max, FIRE_alpha(nrep))
         call writehistory(iteration,FIRE_Ndescend,FIRE_dt,FIRE_alpha,FIRE_vel,nrep,nrestr)
       else
         STOP "optoption should be 0 or 1"
       endif
       call getfilenames(nrep,chi,infile,infile,outfile,iname,rname,oname,avname) !toma ultima foto p/ siguiente paso
       call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,.True.)
       call writenewcoord(oname,rref,boxinfo,natoms,nrestr,mask,per,velout,rav,nrep,nrep,test)
       call writenewcoord(avname,rref,boxinfo,natoms,nrestr,mask,per,velout,ravout,nrep,nrep,test)

       call getmaxdisplacement(nrestr,nrep,rav,rrefall,maxdisp)
       write(9999,*) "Max displacement due MD+steep: ", maxdisp

       write(9999,*) "System converged: F"
    else
       write(9999,*) "System converged: T"
       write(9999,*) "Convergence criteria of ", ftol, " (kcal/mol A) achieved"
    endif


  elseif (nrep .gt. 1) then !NEB on FE surface

    write(9999,*)
    write(9999,*) "Performing NEB on the FE surface"
    write(9999,*)

!------------ Set forces to zero

    fav=0.d0
    devav=0.d0

!------------ Band loop
    if (rextrema) then
      start=2
      nend=nrep-1
    else
      start=1
      nend=nrep
    endif
    do i=start,nend
      call getfilenames(i,chi,infile,reffile,outfile,iname,rname,oname,avname)
      call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin) !saves all coordinates from rname (rst7 file) in rref
      if (rfromtraj) then
        if (usensteps) write(9999,*) " The option usensteps is not available if rfromtraj is True "
        if (usensteps) write(9999,*) " Reading nsteps from feneb.traj file "
        if (allocated(coordall)) deallocate(coordall,coordstat)
        call getcoordfromfenebtraj(nsteps,coordall,nrestr,i)
        allocate(coordstat(nsteps))
      else
        call getdims(iname,nsteps,spatial,natoms)
        if (usensteps) write(9999,*) "Using ", nstepsexternal, "out of ", nsteps
        if (usensteps) nsteps=nstepsexternal
        if (allocated(coordx)) deallocate(coordx)
        if (allocated(coordy)) deallocate(coordy)
        if (allocated(coordz)) deallocate(coordz)
        if (allocated(coordall)) deallocate(coordall,coordstat)
        allocate(coordx(nsteps),coordy(nsteps),coordz(nsteps))
        allocate(coordall(3,nrestr,nsteps),coordstat(nsteps))
        call getcoordfromnetcdf(iname,nsteps,natoms,spatial,coordx,coordy,coordz,coordall,nrestr,mask)
      end if

      tempfilesize=(nsteps-1)
      if(i .eq. start) allocate(temp(tempfilesize,nrep))
      call readtop(topfile,natoms,mask,mass,nrestr)

      if(i .eq. start) then
        write(9999,*) "Reading masses from file: ", trim(topfile)
        j=1
        do while (j .le. (nrestr/3)*3)
          write(9999,*) mass(j), mass(j+1), mass(j+2)
          j = j + 3
        enddo
        write(9999,*) mass(j:nrestr)
        write(9999,*)
      end if

      call getravfav(coordall,nsteps,natoms,nrestr,mask,kref,rav,devav,fav,nrep,i,rref,wgrad,skip,wtemp,dt,mass,tempfilesize,temp)

      if (dostat) then
        write(9999,*)
        write(9999,*) "Statistic stuff"
        write(9999,*) "Replica:", i
        write(9999,*)
        H0T=.True.
        do j=1,nrestr
          atj=mask(j)
          write(9999,*) "Atom:",j
          write(9999,*) "MK test H0"

          coordstat(1:nsteps)=coordall(1,j,1:nsteps)
          call getsstatistics(coordstat,nsteps,skip,nevalfluc,dt,Z,H0,minsegmentlenght,goodrav,gooddevav)
          write(9999,*) "coord x:",H0
          write(111111,*) rav(1,j,i),goodrav
          rav(1,j,i)=goodrav
          fav(1,j,i)=kref*(rav(1,j,i)-rref(1,atj))
          devav(1,j,i)=gooddevav
          H0T=(H0T.and.H0)

          coordstat(1:nsteps)=coordall(2,j,1:nsteps)
          call getsstatistics(coordstat,nsteps,skip,nevalfluc,dt,Z,H0,minsegmentlenght,goodrav,gooddevav)
          write(9999,*) "coord y:",H0
          write(111111,*) rav(2,j,i),goodrav
          rav(2,j,i)=goodrav
          fav(2,j,i)=kref*(rav(2,j,i)-rref(2,atj))
          devav(2,j,i)=gooddevav
          H0T=(H0T.and.H0)

          coordstat(1:nsteps)=coordall(3,j,1:nsteps)
          call getsstatistics(coordstat,nsteps,skip,nevalfluc,dt,Z,H0,minsegmentlenght,goodrav,gooddevav)
          write(9999,*) "coord z:",H0
          write(111111,*) rav(3,j,i),goodrav
          rav(3,j,i)=goodrav
          fav(3,j,i)=kref*(rav(3,j,i)-rref(3,atj))
          devav(3,j,i)=gooddevav
          H0T=(H0T.and.H0)
        end do
        write(*,*) "Trend free:", H0T
        write(9999,*)
      end if ! dostat
    end do

    if (wtemp) then
      open(unit=2203280, file="temperature.dat")
      do i=start,nend
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

        do i=1,nrep
          call writeposdev(rav,devav,nrestr,i,nrep)
        end do

    !----------- Write RMSD
        allocate(rmsd(nrep))
        call getrmsd(fav, kref, nrep, nrestr,rmsd)
        open(unit=40000, file="rmsd.dat")
          write(40000,*) rmsd(1:nrep)
        close(40000)
!----------- Puts reference values in a single array (rrefall).
    do i=start,nend
      call getfilenames(i,chi,infile,reffile,outfile,iname,rname,oname,avname) !rname = NAME_r_i.rst7 ; i=replica
      call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin)
      call getcoordextrema(rref,natoms,rrefall,nrestr,nrep,i,mask)
    end do

!----------- Saves rav in ravout array
    ravout=rav

!----------- Compute the free energy profile by umbrella integration
    allocate(profile(2,nrep))
    if (dostat) call geterror(rav,devav,nrep,nrestr,profile)
    call geterror(rav,devav,nrep,nrestr,profile)
    call getprofile(rav,fav,nrep,nrestr,profile)
    call gettang(rav,tang,nrestr,nrep,tangoption,profile)
    write(9999,*) "Replica", "|", "Max force", "|", "Converged"
    call getnebforce(rav,devav,fav,tang,nrestr,nrep,kspring,maxforceband,ftol,relaxd,&
                    ftrue,ftang,fperp,fspring,.true.,dontg,typicalneb)
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

    call getmaxdisplacement(nrestr,nrep,rav,rrefall,maxdisp)
    write(9999,*) "Max displacement due MD: ", maxdisp

!----------- moves the band
    if (.not. converged) then
      if (.not. typicalneb) then
          if(optoption.eq.0) then
            if (smartstep) then
              write(9999,*) "Using smartstep option"
              write(9999,*) "Base step length: ", stepl
            endif
            do i=2,nrep-1
              if (.not. relaxd) call steep(rav,fperp,nrep,i,steep_size,maxforceband,nrestr,lastmforce,stepl,smartstep)
            end do
            call getmaxdisplacement(nrestr,nrep,rav,rrefall,maxdisp)
            write(9999,*) "Max displacement due MD+steepfperp: ", maxdisp
            write(9999,'(1x,a,f8.6)') "Step length: ", stepl
            if (stepl .lt. 1d-5) then
              write(9999,*)
              write(9999,*) "Warning: max precision reached on atomic displacement"
              write(9999,*) "step length has been set to zero"
              write(9999,*)
            end if
          elseif(optoption.eq.1) then
            do i=2,nrep-1
              if (.not. relaxd) call FIRE(nrestr, rav(:,:,i), fperp(:,:,i), FIRE_vel(:,:,i), FIRE_dt(i),&
              FIRE_Ndescend, FIRE_dt_max, FIRE_alpha(i))
            end do
            call writehistory(iteration,FIRE_Ndescend,FIRE_dt,FIRE_alpha,FIRE_vel,nrep,nrestr)
          else
            STOP "optoption should be 0 or 1"
          end if
          rmsfneb=0.d0
          do i=1,nrep
            call getmaxforce(nrestr,nrep,i,fav,maxforce,ftol,relaxd,maxforceat,rmsfneb)
          end do
          rmsfneb=dsqrt(rmsfneb/dble(nrep*nrestr))
          write(9999,'(1x,a,f8.6)') "rmsfneb(FNEB): ", rmsfneb/nrep


          if (nscycle .eq. 0) write(9999,*) "WARNING: Using only fperp to move the band!"
          if (nscycle .gt. 1) then

            write(9999,*)
            write(9999,*) "Performing extra optimization steps using fspring    "
            write(9999,*) "to get a better distribution of replicas.            "
            write(9999,'(1x,a,I8)') "Extra optmization movements: ", nscycle
            write(9999,*)
          end if

          equispaced=.False.
          k=1

          if (tangoption .eq. 2) then
            write(9999,*) "Tangoption=2 is not available for spring optimization, using tangoption=1 instead"
            write(9999,*) "(the optimization with Fperp will be performed with tangoption=2)"
          end if
          if(.not.tangrecalc) then
            if (tangoption .eq. 2) then
              call gettang(rav,tang,nrestr,nrep,1,profile)
            else
              call gettang(rav,tang,nrestr,nrep,tangoption,profile)
            endif
          endif
          do while ((k .le. nscycle) .and. (.not. equispaced))
            !Computes spring force and others

            if(tangrecalc) then
              if (tangoption .eq. 2) then
                call gettang(rav,tang,nrestr,nrep,1,profile)
              else
                call gettang(rav,tang,nrestr,nrep,tangoption,profile)
              endif
            endif

            call getnebforce(rav,devav,fav,tang,nrestr,nrep,kspring,maxforceband,ftol,converged,&
                            ftrue,ftang,fperp,fspring,.false.,dontg,typicalneb)
            !since wrmforce is false, maxforceband is determined with fsrping
            dontg=0.d0
            ravprevsetp=rav
            maxforcebandprevsetp=maxforceband

            do i=2,nrep-1
              call steep(rav,fspring,nrep,i,steep_spring,maxforcebandprevsetp,nrestr,lastmforce,stepl,.False.)
            end do

            call getdistrightminusleft(rav, nrep, nrestr, equispaced, maxdist)
            if ((k .eq. nscycle) .or. equispaced) then
              write(9999,*)
              write(9999,*) "Band max fspringLast: ", maxforceband
              write(9999,*) "Total spring steps: ", k
              write(9999,*) "Equispaced: ", equispaced
              write(9999,*)
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
          close(unit=1646)

          call getmaxdisplacement(nrestr,nrep,rav,rrefall,maxdisp)

          write(9999,*) "Max displacement due MD+steepfspring: ", maxdisp
      else
        write(9999,*) "Performing a conventional NEB optimization"

        if(optoption.eq.0) then
          if (smartstep) then
            write(9999,*) "Using smartstep option"
            write(9999,*) "Base step length: ", stepl
          endif
          do i=2,nrep-1
            if (.not. relaxd) call steep(rav,fav,nrep,i,steep_size,maxforceband,nrestr,lastmforce,stepl,smartstep)
          end do
          call getmaxdisplacement(nrestr,nrep,rav,rrefall,maxdisp)
          write(9999,*) "Max displacement due MD+steepfperp: ", maxdisp
          write(9999,'(1x,a,f8.6)') "Step length: ", stepl
          if (stepl .lt. 1d-5) then
            write(9999,*)
            write(9999,*) "Warning: max precision reached on atomic displacement"
            write(9999,*) "step length has been set to zero"
            write(9999,*)
          end if
        elseif(optoption.eq.1) then
          do i=2,nrep-1
            if (.not. relaxd) call FIRE(nrestr, rav(:,:,i), fav(:,:,i), FIRE_vel(:,:,i), FIRE_dt(i),&
            FIRE_Ndescend, FIRE_dt_max, FIRE_alpha(i))
          end do
          call writehistory(iteration,FIRE_Ndescend,FIRE_dt,FIRE_alpha,FIRE_vel,nrep,nrestr)
        else
          STOP "optoption should be 0 or 1"
        end if
        rmsfneb=0.d0
        do i=1,nrep
          call getmaxforce(nrestr,nrep,i,fav,maxforce,ftol,relaxd,maxforceat,rmsfneb)
        end do
        rmsfneb=dsqrt(rmsfneb/dble(nrep*nrestr))
        write(9999,'(1x,a,f8.6)') "rmsfneb(FNEB): ", rmsfneb/nrep
      end if

    !------------ Get coordinates for previously optimized extrema
    if (.not. rextrema) then
        !Reactants
        call getfilenames(1,chi,infile,reffile,outfile,iname,rname,oname,avname)
        call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin)
        call getcoordextrema(rref,natoms,rav,nrestr,nrep,1,mask)
        call writenewcoord(oname,rref,boxinfo,natoms,nrestr,mask,per,velout,rav,nrep,1,test)
        call writenewcoord(avname,rref,boxinfo,natoms,nrestr,mask,per,velout,ravout,nrep,1,test)
        !Products
        call getfilenames(nrep,chi,infile,reffile,outfile,iname,rname,oname,avname)
        call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin)
        call getcoordextrema(rref,natoms,rav,nrestr,nrep,nrep,mask)
        call writenewcoord(oname,rref,boxinfo,natoms,nrestr,mask,per,velout,rav,nrep,nrep,test)
        call writenewcoord(avname,rref,boxinfo,natoms,nrestr,mask,per,velout,ravout,nrep,nrep,test)
    end if
    do i=2,nrep-1
      call getfilenames(i,chi,infile,infile,outfile,iname,rname,oname,avname) !toma ultima foto p/ siguiente paso
      call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,.True.)
      call writenewcoord(oname,rref,boxinfo,natoms,nrestr,mask,per,velout,rav,nrep,i,test)
      call writenewcoord(avname,rref,boxinfo,natoms,nrestr,mask,per,velout,ravout,nrep,i,test)
    end do
      call getfilenames(1,chi,infile,infile,outfile,iname,rname,oname,avname)
      call writenewcoord(avname,rref,boxinfo,natoms,nrestr,mask,per,velout,ravout,nrep,1,test)
      call getfilenames(nrep,chi,infile,infile,outfile,iname,rname,oname,avname)
      call writenewcoord(avname,rref,boxinfo,natoms,nrestr,mask,per,velout,ravout,nrep,nrep,test)

    end if !converged
  write(9999,*)
  write(9999,*) "Extrema info"
  write(9999,*) "Barrier: ", barrier
  write(9999,*) "Minimum point: ", minpoint
  write(9999,*) "Maximum point: ", maxpoint
  write(9999,*)

  end if !nrep gt 1

  close(unit=9999)

end program feneb
