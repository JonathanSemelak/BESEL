 program feneb
use netcdf
use readandget
implicit none
character(len=50) :: infile, reffile, outfile, chi, iname, rname, oname, tempname
integer :: nsteps, spatial, natoms, nrestr, nrep, nscycle,maxforceat, rpoint, tgpoint, fpoint
integer :: i, j, k, n, start, end, skip
integer, allocatable, dimension (:) :: mask
real(4) :: coordinate
real(4), allocatable, dimension (:) :: coordx,coordy,coordz
integer, dimension (3) :: point
double precision :: kref, steep_size, ftol, maxforce, kspring, maxforceband, lastmforce, maxforcebandprevsetp, steep_spring
double precision :: stepl, deltaA, rmsfneb, minpoint, maxpoint, barrier
double precision, dimension(6) :: boxinfo
double precision, allocatable, dimension(:) :: rmsd
double precision, allocatable, dimension(:,:) :: rref, profile
double precision, allocatable, dimension(:,:,:) :: rav, fav, tang, ftang, ftrue, fperp, rrefall, ravprevsetp, rcorr, fonref
double precision, allocatable, dimension(:,:,:) :: fspring, dontg, selfdist
logical ::  per, velin, velout, relaxd, converged, wgrad, moved, maxpreached, equispaced, rextrema, test

!------------ Read input
    call readinput(nrep,infile,reffile,outfile,mask,nrestr,lastmforce, &
                 rav,fav,ftrue,ftang,fperp,fspring,tang,kref,kspring,steep_size,steep_spring, &
                 ftol,per,velin,velout,wgrad,rrefall,nscycle,dontg,ravprevsetp,rpoint, tgpoint, &
                 fpoint, rcorr, rextrema, skip)
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


    if (allocated(coordx)) deallocate(coordx)
    if (allocated(coordy)) deallocate(coordy)
    if (allocated(coordz)) deallocate(coordz)
    if (allocated(rref)) deallocate(rref)

    allocate(coordx(nsteps),coordy(nsteps),coordz(nsteps),rref(3,natoms))

    call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin)

    call getavcoordanforces(iname,nsteps,natoms,spatial,coordx,coordy,coordz,&
                        nrestr,mask,kref,rav,fav,nrep,nrep,rref,wgrad,dontg,skip)

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
!------------ Get coordinates for previously optimized extrema
!------------ And set forces to zero


    ! !Reactants
    ! call getfilenames(1,chi,infile,reffile,outfile,iname,rname,oname)
    ! call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin)
    ! call getcoordextrema(rref,natoms,rav,nrestr,nrep,1,mask)
    !
    !
    ! !Products
    ! call getfilenames(nrep,chi,infile,reffile,outfile,iname,rname,oname)
    ! call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin)
    ! call getcoordextrema(rref,natoms,rav,nrestr,nrep,nrep,mask)

    !Forces set to zero
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

      if (allocated(coordx)) deallocate(coordx)
      if (allocated(coordy)) deallocate(coordy)
      if (allocated(coordz)) deallocate(coordz)
      if (allocated(rref)) deallocate(rref)
      allocate(coordx(nsteps),coordy(nsteps),coordz(nsteps),rref(3,natoms))

      call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin)

      call getavcoordanforces(iname,nsteps,natoms,spatial,coordx,coordy, coordz,&
                    nrestr,mask,kref,rav,fav,nrep,i,rref,wgrad,dontg,skip)
    end do

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
      !-----------------test
      call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin)
      ! if (i .le. 9) write(chi,'(I1)') i
      ! if (i .gt. 9 .and. i .le. 99) write(chi,'(I2)') i
      ! if (i .gt. 99 .and. i .le. 999) write(chi,'(I3)') i
      ! tempname = "testrrefall"// trim(chi)
      ! tempname = trim(tempname) // ".rst7"
      ! call writenewcoord(tempname,rref,boxinfo,natoms,nrestr,mask,per,velout,rav,nrep,i,test)
      !----------------tests
      call getcoordextrema(rref,natoms,rrefall,nrestr,nrep,i,mask)
!      call getcoordextrema(rref,natoms,rav,nrestr,nrep,1,mask)

    end do
    test=.False.
    ! allocate(selfdist(2,nrestr,nrep-1))
    !
    ! call getselfdist(rav, rrefall, nrep, nrestr, selfdist)
    !
    ! open(unit=1648, file="selfdist_rrefCHECK.dat")
    ! do i=1,nrestr
    !   do n=1,nrep-1
    !     write(1648,'(2x, I6,2x, f20.10)') n, selfdist(1,i,n)
    !   end do
    !   write(1648,*)
    ! end do
    ! close(1648)
    ! do k=1,nrep
    ! do i=1,nrestr
    ! write(*,'(2x, I6,2x, 3(f20.10,2x))') i, rrefall(1:3,i,k)
    ! end do
    ! end do


!----------- Compute the free energy profile by umbrella integration
    allocate(profile(2,nrep-1))
!----------- Computes tangent and nebforce
    write(9999,*) "-----------------------------------------------------"
    write(9999,*) "r, tg, f points: ", rpoint, tgpoint, fpoint
    write(9999,*) "-----------------------------------------------------"
    if ((rpoint .eq. 0) .and. (tgpoint .eq. 0) .and. (fpoint .eq. 0)) then
      call getprofile(rav,fav,nrep,nrestr,profile)
      call gettang(rav,tang,nrestr,nrep)
      call getnebforce(rav,fav,tang,nrestr,nrep,kspring,maxforceband,ftol,relaxd,&
                          ftrue,ftang,fperp,fspring,.true.,dontg)
    elseif ((rpoint .eq. 1) .and. (tgpoint .eq. 0) .and. (fpoint .eq. 0)) then
      call getprofile(rav,fav,nrep,nrestr,profile)
      call gettang(rav,tang,nrestr,nrep)
      call getnebforce(rav,fav,tang,nrestr,nrep,kspring,maxforceband,ftol,relaxd,&
                          ftrue,ftang,fperp,fspring,.true.,dontg)
      rav=rrefall
    elseif ((rpoint .eq. 0) .and. (tgpoint .eq. 1) .and. (fpoint .eq. 0)) then
      call getprofile(rav,fav,nrep,nrestr,profile)
      call gettang(rrefall,tang,nrestr,nrep)
      call getnebforce(rav,fav,tang,nrestr,nrep,kspring,maxforceband,ftol,relaxd,&
                          ftrue,ftang,fperp,fspring,.true.,dontg)
    elseif ((rpoint .eq. 1) .and. (tgpoint .eq. 1) .and. (fpoint .eq. 0)) then
      call getprofile(rav,fav,nrep,nrestr,profile)
      call gettang(rrefall,tang,nrestr,nrep)
      call getnebforce(rav,fav,tang,nrestr,nrep,kspring,maxforceband,ftol,relaxd,&
                          ftrue,ftang,fperp,fspring,.true.,dontg)
      rav=rrefall
    elseif ((rpoint .eq. 1) .and. (tgpoint .eq. 1) .and. (fpoint .eq. 1)) then
      call gettang(rrefall,tang,nrestr,nrep)
      allocate(fonref(3,nrestr,nrep))
      fonref=0.d0
      do i=2,nrep-1
        do j=1,nrestr
          do k=1,3
            fonref(k,j,i)=fav(k,j,i)*(1-((fav(k,j,i+1)-fav(k,j,i-1))/(rav(k,j,i+1)-rav(k,j,i-1)))*(1.d0/kref))
          end do
        end do
      end do
      fav=fonref
      call getnebforce(rrefall,fav,tang,nrestr,nrep,kspring,maxforceband,ftol,relaxd,&
                          ftrue,ftang,fperp,fspring,.true.,dontg)
      call getprofile(rrefall,fonref,nrep,nrestr,profile)
      rav=rrefall
     elseif ((rpoint .eq. 2) .and. (tgpoint .eq. 2) .and. (fpoint .eq. 2)) then
      call gettang(rrefall,tang,nrestr,nrep)
      do i=2,nrep-1
        do j=1,nrestr
          do k=1,3
            rcorr(k,j,i)=((rav(1,j,i)-rrefall(1,j,i))*tang(1,j,i)+&
                         (rav(2,j,i)-rrefall(2,j,i))*tang(2,j,i)+&
                         (rav(3,j,i)-rrefall(3,j,i))*tang(3,j,i))*tang(k,j,i)
          end do
        end do
      end do
      rcorr=rav+rcorr
      allocate(fonref(3,nrestr,nrep))
      fonref=0.d0

      do i=2,nrep-1
        do j=1,nrestr
          do k=1,3
              fonref(k,j,i)=fav(k,j,i)+((fav(k,j,i+1)-fav(k,j,i-1))/(rav(k,j,i+1)-rav(k,j,i-1)))*(rcorr(k,j,i)-rav(k,j,i))
          end do
        end do
      end do
      ! fonref=fav
      ! write(*,*) fonref
      fav=fonref
      call getprofile(rcorr,fonref,nrep,nrestr,profile)
      call gettang(rcorr,tang,nrestr,nrep)
      call getnebforce(rcorr,fav,tang,nrestr,nrep,kspring,maxforceband,ftol,relaxd,&
                          ftrue,ftang,fperp,fspring,.true.,dontg)
      rav=rcorr
     else
       write(9999,*) "Combination of r,tg and f point not allowed"
       STOP
     end if

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


     ! if (tgpoint .eq. 2) call gettang(rcorr,tang,nrestr,nrep)


     ! if (fpoint .eq. 2) call getnebforce(rav,fav,tang,nrestr,nrep,kspring,maxforceband,ftol,relaxd,&
                      ! ftrue,ftang,fperp,fspring,.true.,dontg)
! fav ---> fneb

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
          write(9999,'(A6,2x,I7,2x,A15,2x,f20.10,2x,A20,2x,f20.10)') "Step: ", k, "Step length", &
          steep_spring, "Band max fspringN: ", maxforceband
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
        !
        ! !Products
        call getfilenames(nrep,chi,infile,reffile,outfile,iname,rname,oname)
        call getrefcoord(rname,nrestr,mask,natoms,rref,boxinfo,per,velin)
        call getcoordextrema(rref,natoms,rav,nrestr,nrep,nrep,mask)
        call writenewcoord(oname,rref,boxinfo,natoms,nrestr,mask,per,velout,rav,nrep,nrep,test)
        !
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
