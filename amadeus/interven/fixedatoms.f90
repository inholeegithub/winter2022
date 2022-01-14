!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       module poscarrefmod
       implicit none
       public
       save
       integer mm
       logical lpbc
       integer nspecies,ndeg
       integer, allocatable :: nelements(:)
       character*2, allocatable :: symbl(:)
       real*8, allocatable :: qqq(:)
        
       end module poscarrefmod
!
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine poscarref(cname,lfault2)
       USE strings, ONLY : parse,value
       USE poscarrefmod, ONLY : mm,nspecies,nelements,symbl,qqq,lpbc,ndeg
       implicit none
       logical lfault2
       character*280 cname
       integer i,j,na,ish
       real*8 scale1,aa1(3),aa2(3),aa3(3),x,y,z,d1,d2,d3,cmatrix(3,3),r6(6),vtest
       character*9 ch9_1,ch9_2
       character*1 ch1
       logical lexist
       character*200 str1
       integer ios,nargs
       character*200 args(40)
       character*20 delims
  
       lfault2=.false.
       inquire(file=trim(cname),exist=lexist)
       if(.not. lexist)then
       write(6,*) 'somehow POSCAR is not present'
       lfault2=.true.
                       return
                       endif
       open(81,file=trim(cname),form='formatted')
!      read(81,*,err=911,end=999) string0
       read(81,'(a200)',err=911,end=999) str1
       delims=' '
       call parse(str1,delims,args,nargs)
       call value(args(1),mm,ios)
       if(ios /= 0) mm=0

       read(81,*,err=911,end=999) scale1
       read(81,*,err=911,end=999) aa1(1),aa1(2),aa1(3)
       read(81,*,err=911,end=999) aa2(1),aa2(2),aa2(3)
       read(81,*,err=911,end=999) aa3(1),aa3(2),aa3(3)
!
       if(scale1 > 0.d0)then
       aa1=aa1*scale1 ; aa2=aa2*scale1 ; aa3=aa3*scale1
                        endif
       if(scale1 < 0.d0)then
       cmatrix(1,:)=aa1(:) ; cmatrix(2,:)=aa2(:) ; cmatrix(3,:)=aa3(:)
       vtest=(cmatrix(1,2)*cmatrix(2,3)-cmatrix(1,3)*cmatrix(2,2))*cmatrix(3,1) &
            +(cmatrix(1,3)*cmatrix(2,1)-cmatrix(1,1)*cmatrix(2,3))*cmatrix(3,2) &
            +(cmatrix(1,1)*cmatrix(2,2)-cmatrix(1,2)*cmatrix(2,1))*cmatrix(3,3)
       vtest=(abs(scale1)/abs(vtest))**(1.d0/3.d0)
       cmatrix=cmatrix*vtest
       aa1(:)=cmatrix(1,:) ; aa2(:)=cmatrix(2,:) ; aa3(:)=cmatrix(3,:)
                        endif
       cmatrix(1,:)=aa1(:) ; cmatrix(2,:)=aa2(:) ; cmatrix(3,:)=aa3(:)
       call latmat00(r6,cmatrix,0)
!
       read(81,'(a200)',err=911,end=999) str1
       delims=' '
       call parse(str1,delims,args,nargs)
       nspecies=nargs
       allocate(symbl(nspecies))
       do i=1,nargs
       symbl(i)=trim(adjustl(args(i)))
       enddo
       allocate(nelements(nspecies))
!      read(81,*,err=911,end=999) 
       read(81,*,err=911,end=999) (nelements(j),j=1,nspecies) 
       na=sum(nelements)
       ndeg=na*3+6
       allocate(qqq(ndeg))
       ish=ndeg-6
       do i=1,6
       qqq(ish+i)=r6(i)
       enddo
       read(81,*,err=911,end=999) ch9_1
       ch1=ch9_1(1:1)
       if(ch1 == 'S') ch9_1='Selective'
       if(ch1 == 's') ch9_1='Selective'
       if(ch1 == 'C') ch9_1='Cartesian'
       if(ch1 == 'c') ch9_1='Cartesian'
       if(ch1 == 'D') ch9_1='Direct'
       if(ch1 == 'd') ch9_1='Direct'
       if(ch9_1 == 'Selective' .or. ch9_1 == 'selective')then
       read(81,*,err=911,end=999) ch9_2
       ch1=ch9_2(1:1)
       if(ch1 == 'C') ch9_2='Cartesian'
       if(ch1 == 'c') ch9_2='Cartesian'
       if(ch1 == 'D') ch9_2='Direct'
       if(ch1 == 'd') ch9_2='Direct'
       if(ch9_2 == 'Cartesian' .or. ch9_2 == 'cartesian')then
       do i=1,na
       read(81,*,err=911,end=999) x,y,z
       if(scale1 > 0.d0)then
       x=x*scale1 ; y=y*scale1 ; z=z*scale1
                        endif
       qqq(3*(i-1)+1)=x ; qqq(3*(i-1)+2)=y ; qqq(3*(i-1)+3)=z
       enddo
       if(lpbc) call tolatx00(ndeg,qqq)
       goto 900
                                                         else
       do i=1,na
       read(81,*,err=911,end=999) d1,d2,d3
       qqq(3*(i-1)+1)=d1 ; qqq(3*(i-1)+2)=d2 ; qqq(3*(i-1)+3)=d3
       enddo
       goto 900
                                                         endif
                                                         else
       if(ch9_1 == 'Cartesian' .or. ch9_1 == 'cartesian')then
       do i=1,na
       read(81,*,err=911,end=999) x,y,z
       if(scale1 > 0.d0)then
       x=x*scale1 ; y=y*scale1 ; z=z*scale1
                        endif
       qqq(3*(i-1)+1)=x ; qqq(3*(i-1)+2)=y ; qqq(3*(i-1)+3)=z
       enddo
       if(lpbc) call tolatx00(ndeg,qqq)
       goto 900
                                                         else
       do i=1,na
       read(81,*,err=911,end=999) d1,d2,d3
       qqq(3*(i-1)+1)=d1 ; qqq(3*(i-1)+2)=d2 ; qqq(3*(i-1)+3)=d3
       enddo
       goto 900
                                                         endif
                                                         endif
  911  continue
  999  continue
       lfault2=.true.
  900  continue
       close(81)
       if(.not. lfault2)then
       cmatrix(1,:)=aa1(:) ; cmatrix(2,:)=aa2(:) ; cmatrix(3,:)=aa3(:)
       call latmat00(r6,cmatrix,0)
       ish=ndeg-6
       do i=1,6
       qqq(ish+i)=r6(i)
       enddo
                        endif
       end
!
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine write_poscar(pname)
       USE poscarrefmod, ONLY : mm,nspecies,nelements,symbl,qqq,ndeg
       implicit none
       character*280 pname
       integer i,j,na,ish
       real*8 r6(6),cmat(3,3),a1(3),a2(3),a3(3)

       ish=ndeg-6
       do i=1,6
       r6(i)=qqq(ish+i)
       enddo
       call latmat00(r6,cmat,1)
       a1(:)=cmat(1,:) ; a2(:)=cmat(2,:) ; a3(:)=cmat(3,:)
       open(71,file=trim(pname),form='formatted')
       write(71,*) mm 
       write(71,'(a3)') '1.0'
       write(71,'(3f23.16)') a1(1),a1(2),a1(3)
       write(71,'(3f23.16)') a2(1),a2(2),a2(3)
       write(71,'(3f23.16)') a3(1),a3(2),a3(3)
       write(71,'(20(2x,a2,1x))') (symbl(i),i=1,nspecies)
       write(71,'(20(i4,1x))') (nelements(i),i=1,nspecies)
       write(71,'(a6)') "Direct"
       na=0
       do i=1,nspecies
       do j=1,nelements(i)
       na=na+1
       write(71,'(3f20.16)') qqq(3*(na-1)+1),qqq(3*(na-1)+2),qqq(3*(na-1)+3)
       enddo
       enddo
       close(71)
       end subroutine write_poscar
!
!234567890
!      Written by In-Ho Lee, KRISS, January 28, 2013.
       subroutine tolatx00(ndeg,qqq)
       implicit none
       integer ndeg
       real*8 qqq(ndeg)
       real*8 b(3,3),devid
       integer j,i,ish,natom
       real*8 t6(6),cmatrix(3,3),a1(3),a2(3),a3(3),d1,d2,d3

       ish=ndeg-6
       natom=ish/3
       do i=1,6
       t6(i)=qqq(ish+i)
       enddo
       call latmat00(t6,cmatrix,1)
       a1(:)=cmatrix(1,:) ; a2(:)=cmatrix(2,:) ; a3(:)=cmatrix(3,:)
!
       devid=a1(1)*a2(2)*a3(3)-a1(2)*a2(1)*a3(3)-a1(1)*a2(3)*a3(2)   &
            +a1(3)*a2(1)*a3(2)+a1(2)*a2(3)*a3(1)-a1(3)*a2(2)*a3(1)
       b(1,1)=-a2(3)*a3(2)+a2(2)*a3(3)
       b(2,1)= a1(3)*a3(2)-a1(2)*a3(3)
       b(3,1)=-a1(3)*a2(2)+a1(2)*a2(3)
       b(1,2)= a2(3)*a3(1)-a2(1)*a3(3)
       b(2,2)=-a1(3)*a3(1)+a1(1)*a3(3)
       b(3,2)= a1(3)*a2(1)-a1(1)*a2(3)
       b(1,3)=-a2(2)*a3(1)+a2(1)*a3(2)
       b(2,3)= a1(2)*a3(1)-a1(1)*a3(2)
       b(3,3)=-a1(2)*a2(1)+a1(1)*a2(2)
       b(:,:)=b(:,:)/devid
!
       do j=1,natom
       d1=b(1,1)*qqq(3*(j-1)+1)+b(1,2)*qqq(3*(j-1)+2)+b(1,3)*qqq(3*(j-1)+3)
       d2=b(2,1)*qqq(3*(j-1)+1)+b(2,2)*qqq(3*(j-1)+2)+b(2,3)*qqq(3*(j-1)+3)
       d3=b(3,1)*qqq(3*(j-1)+1)+b(3,2)*qqq(3*(j-1)+2)+b(3,3)*qqq(3*(j-1)+3)
       qqq(3*(j-1)+1)=d1
       qqq(3*(j-1)+2)=d2
       qqq(3*(j-1)+3)=d3
       enddo
       do j=1,natom
       qqq(3*(j-1)+1)=qqq(3*(j-1)+1)-anint(qqq(3*(j-1)+1))
       qqq(3*(j-1)+2)=qqq(3*(j-1)+2)-anint(qqq(3*(j-1)+2))
       qqq(3*(j-1)+3)=qqq(3*(j-1)+3)-anint(qqq(3*(j-1)+3))
       enddo
       do j=1,natom
       if(qqq(3*(j-1)+1) < 0.d0) qqq(3*(j-1)+1)=qqq(3*(j-1)+1)+1.d0
       if(qqq(3*(j-1)+2) < 0.d0) qqq(3*(j-1)+2)=qqq(3*(j-1)+2)+1.d0
       if(qqq(3*(j-1)+3) < 0.d0) qqq(3*(j-1)+3)=qqq(3*(j-1)+3)+1.d0
       enddo
       end
!
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine latmat00(rlat,wmat,ksign)
       implicit none
       integer ksign
       real*8 rlat(6),wmat(3,3)
       real*8 ra,rb,rc,cosinea,cosineb,cosinec
       real*8 epslat,tmr
       integer i,j

       epslat=1.0d-6
       if(ksign == 1)then
       wmat=0.0d0
       wmat(1,1)=rlat(1)
       wmat(2,1)=rlat(2)*cos(rlat(6))
       wmat(2,2)=rlat(2)*sin(rlat(6))
       wmat(3,1)=rlat(3)*cos(rlat(5))
       wmat(3,2)=rlat(3)*cos(rlat(4))*sin(rlat(6)) &
        -((rlat(3)*cos(rlat(5))-rlat(3)*cos(rlat(4))*cos(rlat(6)))/tan(rlat(6)))
       tmr=rlat(3)**2-wmat(3,1)**2-wmat(3,2)**2  ; if(tmr <= 1.0d-12) tmr=0.0d0
       wmat(3,3)=sqrt(tmr)
       do i=1,3
       do j=1,3
       if(abs(wmat(i,j)) < epslat) wmat(i,j)=0.0d0
       enddo
       enddo
                     else
       rlat=0.0d0
       ra=sqrt(wmat(1,1)**2+wmat(1,2)**2+wmat(1,3)**2)
       rb=sqrt(wmat(2,1)**2+wmat(2,2)**2+wmat(2,3)**2)
       rc=sqrt(wmat(3,1)**2+wmat(3,2)**2+wmat(3,3)**2)
       cosinea=(wmat(2,1)*wmat(3,1)+wmat(2,2)*wmat(3,2)+wmat(2,3)*wmat(3,3))/rb/rc
       cosineb=(wmat(1,1)*wmat(3,1)+wmat(1,2)*wmat(3,2)+wmat(1,3)*wmat(3,3))/rc/ra
       cosinec=(wmat(1,1)*wmat(2,1)+wmat(1,2)*wmat(2,2)+wmat(1,3)*wmat(2,3))/ra/rb
       rlat(1)=ra ; rlat(2)=rb ; rlat(3)=rc
       rlat(4)=acos(cosinea) ; rlat(5)=acos(cosineb) ; rlat(6)=acos(cosinec)
                     endif
       end
!
!234567890
       program fixedatoms
       USE poscarrefmod, ONLY : nelements,symbl,qqq,lpbc,ndeg
       implicit none
       integer i,ish,na0,ndeg0
       logical lfault2
       character*280 cname
       character*280 pname
       real*8, allocatable :: ref(:)

       lfault2=.false.
       lpbc=.true.
       cname='POSCAR_ref'
       cname='../POSCAR_ref'
       call poscarref(cname,lfault2)
       ish=ndeg-6
       do i=1,6
       write(6,*) qqq(ish+i)
       enddo
       na0=sum(nelements)
       ndeg0=na0*3+6
       if(ndeg0 /= ndeg)then
       write(6,*) 'system size mismatch between POSCAR_ref and POSCAR'
                        endif
       do i=1,na0
       write(6,*) qqq(3*(i-1)+1),qqq(3*(i-1)+2),qqq(3*(i-1)+3)
       enddo
       allocate(ref(ndeg0))
       ref=qqq

       pname='POSCAR'
! modification start {
       do i=1,na0
       qqq(3*(i-1)+1)=ref(3*(i-1)+1)
       qqq(3*(i-1)+2)=ref(3*(i-1)+2)
       qqq(3*(i-1)+3)=ref(3*(i-1)+3)
       enddo
! modification end }
! lattice vectors are fixed, here
       call write_poscar(pname)
       if(allocated(symbl)) deallocate(symbl)
       if(allocated(nelements)) deallocate(nelements)
       if(allocated(qqq)) deallocate(qqq)
       if(allocated(ref)) deallocate(ref)
       stop
       end program fixedatoms
!234567890
