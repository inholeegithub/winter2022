!234567890
!      awk '{if(NR ==1) {print "12.34", $0} else { print }}' POSCAR
!      ifort -o poscar2prdf.x strings.f90  poscar2prdf.f90 
!      ./poscar2prdf.x < POSCAR
!      Written by In-Ho Lee, KRISS, April 15, 2016.
       module xprdf
       implicit none
       private
       save
       integer nr,ndim
       real*8 r0,r1,dr,rmax
       real*8, allocatable :: prdf1(:,:),prdf0(:,:)
       integer, allocatable :: irow(:,:)
       integer nspecies
       integer, allocatable :: nelements(:)
       character*2, allocatable :: symbl(:)
       real*8, allocatable :: qqq(:)
     
       public :: init_prdf,final_prdf,do_cal_prdf,get_prdf_init00
       contains
!234567890
!      Written by In-Ho Lee, KRISS, April 15, 2016.
       subroutine init_prdf(rmax00)
       implicit none
       real*8 rmax00

       if(rmax00 <= 0.0d0) rmax00=10.5d0
       rmax=rmax00
       nr=1001 ; r0=0.0d0 ; r1=20.d0      ; dr=(r1-r0)/dble(nr-1)
       if(rmax > 20.d0)then
       nr=1001 ; r0=0.0d0 ; r1=rmax*2.0d0 ; dr=(r1-r0)/dble(nr-1)
                       endif
!      write(6,'(2f18.8,1x,a7)') rmax,'rmax,dr'
!      write(6,'(3f18.8)') r0,r1,dr
!      write(6,'(i5)') nr
       end subroutine init_prdf
!234567890
!      Written by In-Ho Lee, KRISS, April 15, 2016.
       subroutine final_prdf()
       implicit none

       if(allocated(prdf0)) deallocate(prdf0)
       if(allocated(prdf1)) deallocate(prdf1)
       end subroutine final_prdf
!234567890
!      Written by In-Ho Lee, KRISS, April 15, 2016.
       subroutine do_cal_prdf()
       implicit none
       real*8 value1,value2,rr
       integer j,ir,iiss

       call get_prdf_final() ; prdf1=prdf0
       do j=1,ndim
       call smooth(nr,prdf1(1,j))
       enddo
       end subroutine do_cal_prdf
!234567890
!      Written by In-Ho Lee, KRISS, April 13, 2016.
       real*8 function theta(x)
       implicit none
       real*8 x
      
       theta=0.0d0 
       if(x >= 0.0d0)theta=1.0d0 
       end function theta
!234567890
!      Written by In-Ho Lee, KRISS, April 15, 2016.
       subroutine smooth(n,a)
       implicit none
       integer n
       real*8 a(n)
       real*8, allocatable :: b(:)
       integer i

       allocate(b(n))
       b=0.d0
       do i=3,n-2
       b(i-2)=(69.d0*a(i-2)+4.d0*a(i-1)-6.d0*a(i)+4.d0*a(i+1)-a(i+2))/70.d0
       b(i+2)=(69.d0*a(i+2)+4.d0*a(i+1)-6.d0*a(i)+4.d0*a(i-1)-a(i-2))/70.d0
       b(i-1)=(2.d0*a(i-2)+27.d0*a(i-1)+12.d0*a(i)-8.d0*a(i+1)+2.d0*a(i+2))/35.d0
       b(i+1)=(2.d0*a(i+2)+27.d0*a(i+1)+12.d0*a(i)-8.d0*a(i-1)+2.d0*a(i-2))/35.d0
       b(i  )=(-3.d0*a(i-2)+12.d0*a(i-1)+17.d0*a(i)+12.d0*a(i+1)-3.d0*a(i+2))/35.d0
       enddo
       a=b
       deallocate(b)
       end subroutine smooth
!234567890
!      Written by In-Ho Lee, KRISS, April 13, 2016.
       subroutine get_prdf_final()
       implicit none
       real*8 a1(3),a2(3),a3(3),aa1(3),aa2(3),aa3(3),t6(6),cmatrix(3,3)
       real*8, allocatable :: dir(:,:),dirext(:,:)
       integer i,j,k,ie,je,ke,k0,kk0,natot,ncext(3)
       real*8 vec(3),wec(3),vtest
       integer, allocatable :: itype(:),nelementsext(:)
       integer i1,i2,i3,j1,j2,k1,k2,ir,ip,iq,ks,ksext,ndeg,ish
       real*8 rr,pi,tmr

       call get_prdf_init()
       ndeg=3*sum(nelements)+6
       ish=ndeg-6
       do i=1,6
       t6(i)=qqq(ish+i)
       enddo
       natot=ish/3
!      write(6,*) ndeg,natot
       call latmat(t6,cmatrix,1)
       a1(:)=cmatrix(1,:) ; a2(:)=cmatrix(2,:) ; a3(:)=cmatrix(3,:)
       call get_extension(a1,a2,a3,rmax,ncext)
       write(6,'(a1,1x,3i5)') '#',ncext(1),ncext(2),ncext(3)
       ncext=ncext*2
!      write(6,*) ncext
!
       allocate(dir(natot,3))
       do i=1,natot
       dir(i,1)=qqq(3*(i-1)+1)
       dir(i,2)=qqq(3*(i-1)+2)
       dir(i,3)=qqq(3*(i-1)+3)
       enddo
       deallocate(qqq) 
       aa1=a1*dble(ncext(1)) ; aa2=a2*dble(ncext(2)) ; aa3=a3*dble(ncext(3))
       k=sum(nelements)*(ncext(1)*ncext(2)*ncext(3))
       allocate(nelementsext(nspecies)) ; allocate(itype(k))
       nelementsext(:)=nelements(:)*(ncext(1)*ncext(2)*ncext(3))
       ksext=sum(nelementsext)
!      write(6,*)k,ksext
       allocate(dirext(ksext,3))
       cmatrix(1,:)=aa1(:) ; cmatrix(2,:)=aa2(:) ; cmatrix(3,:)=aa3(:)
       vtest=(cmatrix(1,2)*cmatrix(2,3)-cmatrix(1,3)*cmatrix(2,2))*cmatrix(3,1) &
            +(cmatrix(1,3)*cmatrix(2,1)-cmatrix(1,1)*cmatrix(2,3))*cmatrix(3,2) &
            +(cmatrix(1,1)*cmatrix(2,2)-cmatrix(1,2)*cmatrix(2,1))*cmatrix(3,3)
       vtest=abs(vtest)
       pi=4.0d0*atan(1.0d0)
       k=0 ; kk0=0
       do i=1,nspecies
       k0=0
       do j=1,nelements(i)
       k0=k0+1 ; kk0=kk0+1
       do ie=0,ncext(1)-1
       do je=0,ncext(2)-1
       do ke=0,ncext(3)-1
       k=k+1
       dirext(k,1)=(dble(ie)+dir(kk0,1))/dble(ncext(1))
       dirext(k,2)=(dble(je)+dir(kk0,2))/dble(ncext(2))
       dirext(k,3)=(dble(ke)+dir(kk0,3))/dble(ncext(3))
       itype(k)=i
       enddo
       enddo
       enddo
       enddo
       if(k0 /= nelements(i))then
       write(6,*) 'error, k0',i
                             stop
                             endif
       enddo
       if(k /= ksext)then
       write(6,*) 'error ksext',k,ksext
                     stop
                     endif
       deallocate(symbl) ; deallocate(nelements)  
       deallocate(dir) 
       i3=0
       do i1=1,nspecies
       do i2=1,nspecies
       i3=i3+1
       if(ndim < i3)then
       write(6,*) 'check ndim',ndim
                    stop
                    endif
       irow(i1,i2)=i3
       prdf0(:,i3)=0.0d0
       enddo
       enddo
       k1=0
       do i1=1,nspecies
       do j1=1,nelementsext(i1)
       k1=k1+1
       if(k1 > ksext) stop
       k2=0
       do i2=1,nspecies
       i3=irow(i1,i2)
       do j2=1,nelementsext(i2)
       k2=k2+1
       if(k2 > ksext) stop
       if(k1 == k2) cycle
       if(itype(k1) == i1 .and. itype(k2) == i2)then
       vec(:)=dirext(k1,:)-dirext(k2,:) 
       vec(1)=vec(1)-anint(vec(1))
       vec(2)=vec(2)-anint(vec(2))
       vec(3)=vec(3)-anint(vec(3))
       wec(:)=vec(1)*aa1(:)+vec(2)*aa2(:)+vec(3)*aa3(:)
       tmr=sqrt(dot_product(wec,wec))
       ir=(1+int(tmr/dr))-2 ; ip=max(ir,2) ; ir=(1+int(tmr/dr))+2 ; iq=min(ir,nr-1)
       do ir=ip,iq
       rr=r0+dr*float(ir-1)
       prdf0(ir,i3)=prdf0(ir,i3)+theta(tmr-rr)*theta(rr+dr-tmr)*dble(ksext)/(dble(ksext)/vtest) &
                              /(dble(nelementsext(i1)))/(dble(nelementsext(i2)))
       enddo
                                                endif
       enddo
       enddo 
       enddo
       enddo
       deallocate(dirext) ; deallocate(nelementsext) ; deallocate(itype)
       do i1=1,nspecies
       do i2=1,nspecies
       i3=irow(i1,i2)
       write(6,'(a1,1x,2i4)') '#', i1,i2
       do ir=2,nr
       rr=r0+dr*float(ir-1)
       prdf0(ir,i3)=prdf0(ir,i3)/(dr*4.0d0*pi*rr**2)
       if(rr <= rmax) write(6,'(f12.6,f22.12)') rr,prdf0(ir,i3)
       enddo
       write(6,*) '&'
       enddo
       enddo
       deallocate(irow)
       end subroutine get_prdf_final
!234567890
!      Written by In-Ho Lee, KRISS, April 15, 2016.
       subroutine get_prdf_init00(nspecies0,nelements0,symbl0,ndeg0,qqq0,rmax0)
       implicit none
       integer nspecies0,ndeg0,nelements0(nspecies0)
       real*8 qqq0(ndeg0),rmax0
       character*2 symbl0(nspecies0)
       integer i,j,ks

       if(rmax0 <= 0.0d0) rmax0=10.0d0
       call init_prdf(rmax0)
       nspecies=nspecies0
       allocate(nelements(nspecies)) ; allocate(symbl(nspecies))
       allocate(irow(nspecies,nspecies))
       ndim=0
       do i=1,nspecies
       do j=1,nspecies
       ndim=ndim+1
       irow(i,j)=ndim
       enddo
       enddo
       if(.not. allocated(prdf0)) allocate(prdf0(nr,ndim))
       if(.not. allocated(prdf1)) allocate(prdf1(nr,ndim))
       do i=1,nspecies
       symbl(i)=trim(adjustl(symbl0(i)))
       nelements(i)=nelements0(i)
       enddo
       ks=sum(nelements)
       if(ks*3+6 /= ndeg0)then
       write(6,*) 'ks*3+6, ndeg0'
       write(6,*)  ks*3+6, ndeg0
                          stop
                          endif
       allocate(qqq(ndeg0))
       do j=1,ndeg0
       qqq(j)=qqq0(j)
       enddo
       end subroutine get_prdf_init00
!234567890
!      Written by In-Ho Lee, KRISS, April 15, 2016.
       subroutine get_prdf_init()
       USE strings, ONLY : parse,value
       implicit none
       real*8 scale0,amatrix(3,3),t6(6),a1(3),a2(3),a3(3),rmax00,vtest
       real*8, allocatable :: dir(:,:),car(:,:)
       integer i,j,ks
       integer ios,nargs
       character*200 str1
       character*200 args(40)
       character*20 delims
       logical lfault

       lfault=.false.
       read(5,*) rmax00
       rmax=rmax00
       call init_prdf(rmax00)
       read(5,*) scale0
       read(5,*) amatrix(1,1),amatrix(1,2),amatrix(1,3)
       read(5,*) amatrix(2,1),amatrix(2,2),amatrix(2,3)
       read(5,*) amatrix(3,1),amatrix(3,2),amatrix(3,3)
       if(scale0 < 0.d0)then
       vtest=(amatrix(1,2)*amatrix(2,3)-amatrix(1,3)*amatrix(2,2))*amatrix(3,1) &
            +(amatrix(1,3)*amatrix(2,1)-amatrix(1,1)*amatrix(2,3))*amatrix(3,2) &
            +(amatrix(1,1)*amatrix(2,2)-amatrix(1,2)*amatrix(2,1))*amatrix(3,3)
       vtest=abs(vtest)
       vtest=abs(scale0)/vtest ; vtest=vtest**(1.d0/3.d0)
       amatrix=amatrix*vtest
                        endif
       if(scale0 > 0.d0)then
       amatrix=amatrix*scale0 
                        endif
       read(5,'(a200)',err=911,end=999) str1
       delims=' '
       call parse(str1,delims,args,nargs)
       nspecies=nargs
       allocate(nelements(nspecies)) ; allocate(symbl(nspecies))
       allocate(irow(nspecies,nspecies))
!      write(6,*) nspecies
       ndim=0
       do i=1,nspecies
       do j=1,nspecies
       ndim=ndim+1
       irow(i,j)=ndim
       enddo
       enddo
       if(.not. allocated(prdf0)) allocate(prdf0(nr,ndim))
       if(.not. allocated(prdf1)) allocate(prdf1(nr,ndim))
       do i=1,nspecies
       symbl(i)=trim(adjustl(args(i)))
       enddo
!      write(6,'(10(1x,a2))') (symbl(i),i=1,nspecies)
       read(5,'(a200)',err=911,end=999) str1
       delims=' '
       call parse(str1,delims,args,nargs)
       do i=1,nspecies
       call value(args(i),nelements(i),ios)
       enddo
       ks=sum(nelements)
!      write(6,'(10(1x,i4))') (nelements(i),i=1,nspecies)
!      write(6,*) ks
       allocate(qqq(ks*3+6))
       allocate(dir(ks,3)) ; allocate(car(ks,3))
       read(5,'(a200)',err=911,end=999) str1
       delims=' '
       call parse(str1,delims,args,nargs)
       if(nargs > 0)then
       if(args(1) == 'DIR' .or. args(1) == 'dir' .or. args(1) == 'D' .or. args(1) == 'd' .or. &
          args(1) == 'direct' .or. args(1) == 'Direct' .or. args(1) == 'Dir')then
       do j=1,ks
       read(5,*) dir(j,1),dir(j,2),dir(j,3)
       dir(j,1)=dir(j,1)-anint(dir(j,1))
       dir(j,2)=dir(j,2)-anint(dir(j,2))
       dir(j,3)=dir(j,3)-anint(dir(j,3))
       if(dir(j,1) < 0.d0) dir(j,1)=dir(j,1)+1.d0
       if(dir(j,2) < 0.d0) dir(j,2)=dir(j,2)+1.d0
       if(dir(j,3) < 0.d0) dir(j,3)=dir(j,3)+1.d0
       car(j,1)=dir(j,1)*amatrix(1,1)+dir(j,2)*amatrix(2,1)+dir(j,3)*amatrix(3,1)
       car(j,2)=dir(j,1)*amatrix(1,2)+dir(j,2)*amatrix(2,2)+dir(j,3)*amatrix(3,2)
       car(j,3)=dir(j,1)*amatrix(1,3)+dir(j,2)*amatrix(2,3)+dir(j,3)*amatrix(3,3)
       enddo
                                                                             else
       do j=1,ks
       read(5,*) car(j,1),car(j,2),car(j,3)
       enddo
       dir=car
       a1(:)=amatrix(1,:) ; a2(:)=amatrix(2,:) ; a3(:)=amatrix(3,:)
       call tolaty(dir,a1,a2,a3,ks)
                                                                             endif
                    endif
  911  continue
       lfault=.true.
  999  continue
       do j=1,ks
       qqq(3*(j-1)+1)=dir(j,1)
       qqq(3*(j-1)+2)=dir(j,2)
       qqq(3*(j-1)+3)=dir(j,3)
       enddo
       call latmat(t6,amatrix,0)
       do j=1,6
       qqq(ks*3+j)=t6(j)
       enddo
       deallocate(dir,car)
       end subroutine get_prdf_init

       end module xprdf
!234567890
!      Written by In-Ho Lee, KRISS, April 13, 2016.
       subroutine get_extension(a1,a2,a3,rmax0,ncext)
       implicit none
       real*8 a1(3),a2(3),a3(3),rmax0
       integer ncext(3)
       real*8 v(3),h(3)

       call cross3(a1,a2,v) ; v=v/sqrt(sum(v*v)) ; h(3)=abs(sum(v*a3))
       call cross3(a3,a1,v) ; v=v/sqrt(sum(v*v)) ; h(2)=abs(sum(v*a2))
       call cross3(a2,a3,v) ; v=v/sqrt(sum(v*v)) ; h(1)=abs(sum(v*a1))
       v=rmax0/h+0.5d0
       ncext=nint(v)
       end
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine cross3(a,b,c)
       implicit none
       real*8 a(3),b(3),c(3)

       c(1)=a(2)*b(3)-a(3)*b(2)
       c(2)=a(3)*b(1)-a(1)*b(3)
       c(3)=a(1)*b(2)-a(2)*b(1)
       end
!234567890
!      Written by In-Ho Lee, KRISS, March 27, 2016
       subroutine tolaty(xyz,a1,a2,a3,ktot)
       implicit none
       integer ktot
       real*8 xyz(ktot,3),a1(3),a2(3),a3(3)
       real*8 b(3,3),devid
       integer j,i
       real*8 d1,d2,d3

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
       do j=1,ktot
       d1=b(1,1)*xyz(j,1)+b(1,2)*xyz(j,2)+b(1,3)*xyz(j,3)
       d2=b(2,1)*xyz(j,1)+b(2,2)*xyz(j,2)+b(2,3)*xyz(j,3)
       d3=b(3,1)*xyz(j,1)+b(3,2)*xyz(j,2)+b(3,3)*xyz(j,3)
       xyz(j,1)=d1
       xyz(j,2)=d2
       xyz(j,3)=d3
       enddo
       do j=1,ktot
       xyz(j,1)=xyz(j,1)-anint(xyz(j,1))
       xyz(j,2)=xyz(j,2)-anint(xyz(j,2))
       xyz(j,3)=xyz(j,3)-anint(xyz(j,3))
       enddo
       do j=1,ktot
       if(xyz(j,1) < 0.d0) xyz(j,1)=xyz(j,1)+1.d0
       if(xyz(j,2) < 0.d0) xyz(j,2)=xyz(j,2)+1.d0
       if(xyz(j,3) < 0.d0) xyz(j,3)=xyz(j,3)+1.d0
       enddo
       end
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine latmat(rlat,wmat,ksign)
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
!234567890
!      Written by In-Ho Lee, KRISS, April 15, 2016.
       program prdf_cal
       USE xprdf
       implicit none

       call do_cal_prdf()
       call final_prdf()
       end program prdf_cal
!234567890
