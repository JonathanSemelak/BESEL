subroutine writeposforces(rav,fav,nrestr,rep,nrep)

implicit none
integer, intent(in) :: nrestr, rep, nrep
double precision, dimension(3,nrestr,nrep), intent(in) :: rav, fav
integer :: i

open(unit=1001, file="feneb.gradients", position='append')
do i=1,nrestr
write(1001,'(2x, I6,2x, 6(f20.10,2x))') i, rav(1:3,i,rep), fav(1:3,i,rep)
end do
write(1001,'(2x, I6,2x, 6(f20.10,2x))')
close(1001)

end subroutine writeposforces

subroutine writeposxyz(rav,nrestr,rep,nrep)

implicit none
integer, intent(in) :: nrestr, rep, nrep
double precision, dimension(3,nrestr,nrep), intent(in) :: rav
integer :: i

open(unit=10011, file="Pos_average.xyz", position='append')
do i=1,nrestr
write(10011,'(2x, I6,2x, 6(f20.10,2x))') i, rav(1:3,i,rep)
end do
write(10011,'(2x, I6,2x, 6(f20.10,2x))')
close(10011)

end subroutine writeposxyz

subroutine writeposdev(rav,devav,nrestr,rep,nrep)

implicit none
integer, intent(in) :: nrestr, rep, nrep
double precision, dimension(3,nrestr,nrep), intent(in) :: rav, devav
integer :: i

open(unit=1002, file="Pos_dev.dat", position='append')
do i=1,nrestr
write(1002,'(2x, I6,2x, 6(f20.10,2x))') i, rav(1:3,i,rep), devav(1:3,i,rep)
end do
write(1002,'(2x, I6,2x, 6(f20.10,2x))')
close(1002)

end subroutine writeposdev

subroutine writenewcoord(oname,rref,boxinfo,natoms,nrestr,mask,per,velout,rav,nrep,rep,test)

implicit none
character(len=500), intent(in) :: oname
integer, intent(in) :: natoms, nrestr, nrep, rep
double precision, dimension(3,nrestr,nrep) :: rav
double precision, dimension(3,natoms), intent(in) :: rref
double precision, dimension(3,natoms) :: rout
double precision, dimension(6), intent(in) :: boxinfo
integer, dimension(nrestr), intent(in) :: mask
logical, intent(in) :: per, velout, test
integer :: i, j, at, auxunit

rout=rref
auxunit=21000000+rep
if (.not. test) then
do i=1,nrestr
  do j=1,3
    at=mask(i)
    rout(j,at)=rav(j,i,rep)
  end do
end do
end if

open (unit=auxunit, file=oname)
write(auxunit,*) "FENEB restart, replica: ", rep
!write(auxunit,'(I8)') natoms
write(auxunit,'(I6)') natoms
i=1
do while (i .le. natoms/2)
  write(auxunit,'(6(f12.7))') rout(1,2*i-1), rout(2,2*i-1), rout(3,2*i-1), &
                          rout(1,2*i), rout(2,2*i), rout(3,2*i)
  i = i + 1
enddo
if (mod(natoms,2) .ne. 0) write(auxunit,'(3(f12.7))') rout(1:3,2*i-1)
if (velout) then
  i=1
  do while (i .le. natoms/2)
    write(auxunit,'(6(f12.7))') 0.d0,0.d0,0.d0,0.d0,0.d0,0.d0
    i = i + 1
  enddo
  if (mod(natoms,2) .ne. 0) write(auxunit,'(3(f12.7))') 0.d0, 0.d0, 0.d0
endif
if (per) write(auxunit,'(6(f12.7))') boxinfo(1:6)
close (unit=auxunit)

end subroutine writenewcoord

subroutine writehistory(iteration,FIRE_Ndescend,FIRE_dt,FIRE_alpha,FIRE_vel,nrep,nrestr)
implicit none
integer, intent(in) :: nrep, nrestr
integer, intent(in) :: FIRE_Ndescend, iteration
double precision, dimension(nrep), intent(in):: FIRE_dt, FIRE_alpha
double precision, dimension(3,nrestr,nrep), intent(in) :: FIRE_vel
integer :: i ,j

  write(9999,*) "Writing feneb.present"
  open (unit=22000112, file="feneb.present", position='append')
  write(22000112,*) iteration,FIRE_Ndescend
  do i=1,nrep
    write(22000112,*) FIRE_dt(i),FIRE_alpha(i)
  end do
  do i=1,nrestr
    do j=1,nrep
      write(22000112,*) FIRE_vel(1,i,j),FIRE_vel(2,i,j),FIRE_vel(3,i,j)
    end do
  end do
  close(unit=22000112)

end subroutine writehistory
