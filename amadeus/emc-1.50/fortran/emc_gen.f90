! Copyright (c) "2012, Alexandr Fonari
!                URL: https://github.com/alexandr-fonari/emc
!                License: MIT License"

program emc_gen ! version 1.0
use emc_functions
implicit none
real(kind=8), parameter :: b2a = 0.52917721092d0
real(kind=8), parameter :: pi = 4.d0*DATAN(1.d0)
integer(kind=4), parameter :: nkpoints = 61
integer(kind=4), parameter :: iunt = 10, ilog = 11, w = 1 ! k point weight for VASP
character(len=8), parameter :: version_number = '1.5f'

real(kind=8) :: kp(3), kpc(3), dk, E(-2:2,-2:2,-2:2), A(4), f(3,3), g(3,3)
integer(kind=4) :: i, i1, j, j1, band
character prg
real(kind=8) :: a1(3),a2(3),a3(3),b1(3),b2(3),b3(3),omega,xkpc(3)

open(unit=ilog,file='emc_gen.log',form='formatted')
write(ilog,*) "Effective Mass Calculator generator ", version_number
call print_time(ilog)
write(ilog,*)

open(unit=iunt,file='inp',form='formatted')
    read(iunt,fmt=*) (kp(i),i=1,size(kp))
    read(iunt,fmt=*) dk
    read(iunt,fmt=*) band
    read(iunt,fmt=*) prg
    read(iunt,*) ((f(i,j),j=1,3),i=1,3)
close(iunt)

write(ilog,"(A10,I5,F12.6)") "band, dk: ", band, dk

if(prg .eq. 'V') then
    dk = dk/(2.0D0*pi*b2a)
    write(ilog,"(A44,F12.6)") "dk will be converted to VASP units (2Pi/A): ", dk
end if
write(ilog,*)


! g = inverse(transpose(f)), will do it in two steps
g=f
g=transpose(g)
call inverse(g, 3)
g=g*2.d0*pi
!234567890
       i=1
       i=0
       if(i == 1)then
       a1(:)=f(1,:) ; a2(:)=f(2,:) ; a3(:)=f(3,:)
       call cross3(a2,a3,b1) ; call cross3(a3,a1,b2) ; call cross3(a1,a2,b3)
       omega=abs(dot_product(b1,a1))
       b1=b1*(2.d0*pi/omega) ; b2=b2*(2.d0*pi/omega) ; b3=b3*(2.d0*pi/omega)
       write(6,*) maxval(abs(g(1,:)-b1(:))),' should be zero'
       write(6,*) maxval(abs(g(2,:)-b2(:))),' should be zero'
       write(6,*) maxval(abs(g(3,:)-b3(:))),' should be zero'
                 endif
!234567890
write(ilog,*) "direct lattice vectors              reciprocal lattice vectors"
do i=1,3
    write(ilog,"(3F10.6,A5,3F10.6)") , (f(i,j), j=1,3), "     ", (g(i,j), j=1,3)
end do
write(ilog,*)

kpc = pureDEGMV(g, kp, 'T')
!234567890
       i=1
       i=0
       if(i == 1)then
       xkpc(:)=kp(1)*b1(:)+kp(2)*b2(:)+kp(3)*b3(:)
       write(6,*) kpc
       write(6,*) xkpc
       write(6,*) maxval(abs(xkpc(:)-kpc(:))),' should be zero'
                 endif
!234567890
write(ilog,*) "k-point in: reciprocal space       reciprocal Cartesian space"
write(ilog,"(3F10.6,A5,3F10.6)") (kp(j), j=1,3), "     ", (kpc(j), j=1,3)

! write KPOINTS file ###########################################

open(unit=iunt,file='KPOINTS',form='formatted')

write(unit=iunt,fmt='(A5,F9.4,F9.4,F9.4,A6,F9.5)') 'kpc: ', (kpc(j), j=1,3), ', dk: ', dk
write(unit=iunt,fmt='(I5)') nkpoints
if(prg .eq. 'C') then
    write(unit=iunt,fmt='(A)') 'Reciprocal'
else
    write(unit=iunt,fmt='(A)') 'Cartesian'
end if
!

call set_next_eigeval(iunt, g, kpc, 0, 0, 0, w, dk, prg) ! 0 0 0

! http://stackoverflow.com/questions/9791001/loop-in-fortran-from-a-list
A = (/ -2, -1, 1, 2 /)
! x
do i = 1, size(A)
    i1 = A(i)
    call set_next_eigeval(iunt, g, kpc, i1, 0, 0, w, dk, prg)
end do

! y
do i = 1, size(A)
    i1 = A(i)
    call set_next_eigeval(iunt, g, kpc, 0, i1, 0, w, dk, prg)
end do

! z
do i = 1, size(A)
    i1 = A(i)
    call set_next_eigeval(iunt, g, kpc, 0, 0, i1, w, dk, prg)
end do

! xy
do i = 1, size(A)
    i1 = A(i)
    do j=1, size(A)
        j1 = A(j)
        call set_next_eigeval(iunt, g, kpc, j1, i1, 0, w, dk, prg)
    end do
end do

! xz
do i = 1, size(A)
    i1 = A(i)
    do j=1, size(A)
        j1 = A(j)
        call set_next_eigeval(iunt, g, kpc, j1, 0, i1, w, dk, prg)
    end do
end do

! yz
do i = 1, size(A)
    i1 = A(i)
    do j=1, size(A)
        j1 = A(j)
        call set_next_eigeval(iunt, g, kpc, 0, j1, i1, w, dk, prg)
    end do
end do

close(iunt)
end program

subroutine set_next_eigeval(iunt, g, kp, i, j, k, w, dk, prg)
    use emc_functions
    implicit none
    integer(kind=4), intent(in) :: iunt, i, j, k, w
    real(kind=8), intent(in) :: g(3,3), kp(3), dk
    character, intent(in) :: prg
    real(kind=8) :: kp1(3)

    kp1(1) = kp(1)+dble(i)*dk
    kp1(2) = kp(2)+dble(j)*dk
    kp1(3) = kp(3)+dble(k)*dk

    if(prg .eq. 'C') then
        kp1 = pureDGESV(g, kp1, 'T') ! get reciprocal coords
    end if
    
    write(unit=iunt,fmt='(F22.12,F22.12,F22.12,2x,F12.8)') kp1(1), kp1(2), kp1(3), dble(w)
    ! write(unit=iunt,fmt='(I5,I5,I5)') i, j, k
    return
end subroutine
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine cross3(a,b,c)
       implicit none
       real(kind=8), intent(in) ::  a(3),b(3)
       real(kind=8), intent(out) :: c(3)

       c(1)=a(2)*b(3)-a(3)*b(2)
       c(2)=a(3)*b(1)-a(1)*b(3)
       c(3)=a(1)*b(2)-a(2)*b(1)
       return
       end
!234567890
