subroutine writeposforces(rav,fav,nrestr,rep,nrep)

implicit none
integer, intent(in) :: nrestr, rep, nrep
double precision, dimension(3,nrestr,nrep), intent(in) :: rav, fav
integer :: i

open(unit=1001, file="Pos_forces.dat", position='append')
do i=1,nrestr
write(1001,'(2x, I6,2x, 6(f20.10,2x))') i, rav(1:3,i,rep), fav(1:3,i,rep)
end do
write(1001,'(2x, I6,2x, 6(f20.10,2x))')
close(1001)

end subroutine writeposforces

subroutine writenewcoord(oname,rref,boxinfo,natoms,nrestr,mask,per,vel,rav,nrep,rep)

implicit none
character(len=50), intent(in) :: oname
integer, intent(in) :: natoms, nrestr, nrep, rep
double precision, dimension(3,nrestr,nrep) :: rav
double precision, dimension(3,natoms), intent(in) :: rref
double precision, dimension(3,natoms) :: rout
double precision, dimension(6), intent(in) :: boxinfo
integer, dimension(nrestr), intent(in) :: mask
logical, intent(in) :: per, vel
integer :: i, j, at, auxunit

rout=rref
auxunit=21000000+rep

do i=1,nrestr
  do j=1,3
    at=mask(i)
    rout(j,at)=rav(j,i,rep)
  end do
end do

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
if (mod(natoms,2) .ne. 0) write(auxunit,'(6(f12.7))') rout(1:3,2*i)
if (vel) then
  i=1
  do while (i .le. natoms/2)
    write(auxunit,'(6(f12.7))') 0.d0,0.d0,0.d0,0.d0,0.d0,0.d0
    i = i + 1
  enddo
 endif
  if (mod(natoms,2) .ne. 0) write(auxunit,'(6(f12.7))') 0.d0, 0.d0, 0.d0
  if (per) write(auxunit,'(6(f12.7))') boxinfo(1:6)
close (unit=auxunit)

end subroutine writenewcoord
