!234567890
!      Written by In-Ho Lee, KRISS, January 20, 2019.
       program makevertical
       implicit none
       integer n1,n2,n3
       real*8 a1(3),a2(3),a3(3),cmat0(3,3),cmat2(3,3),cellvol0,vtest
       real*8 b1(3),b2(3),b3(3),height,pi,tmp,tmq
       real*8 cmat3(3,3),scale0,uec(3),vec(3),wec(3),rlat(6),zmat(3,3)
       real*8, allocatable :: dir(:,:)
       real*8, allocatable :: cartext(:,:),pppext(:,:)
       real*8 tmpv(100)
       integer nspecies,itmpv(100)
       integer iatom,natom0,i,j,k,ipos,natomext
       integer, allocatable :: nelements(:)
       character*2, allocatable :: symbl(:)
       character*20 ich20
       character*1 ich1
       character*280 string,char1(100)
       character*280 str0

       pi=4.d0*atan(1.d0)
!      original cell  [
       read(5,'(a280)') str0
       read(5,*) scale0
       read(5,*) a1(1),a1(2),a1(3)
       read(5,*) a2(1),a2(2),a2(3)
       read(5,*) a3(1),a3(2),a3(3)
       if(scale0 > 0.d0)then
       cmat0(1,:)=a1(:) ; cmat0(2,:)=a2(:) ; cmat0(3,:)=a3(:)
       cmat0=cmat0*scale0
       a1(:)=cmat0(1,:) ; a2(:)=cmat0(2,:) ; a3(:)=cmat0(3,:)
                        endif
       if(scale0 < 0.d0)then
       cmat0(1,:)=a1(:) ; cmat0(2,:)=a2(:) ; cmat0(3,:)=a3(:)
       vtest=(cmat0(1,2)*cmat0(2,3)-cmat0(1,3)*cmat0(2,2))*cmat0(3,1) &
            +(cmat0(1,3)*cmat0(2,1)-cmat0(1,1)*cmat0(2,3))*cmat0(3,2) &
            +(cmat0(1,1)*cmat0(2,2)-cmat0(1,2)*cmat0(2,1))*cmat0(3,3)
       vtest=(abs(scale0)/abs(vtest))**(1.d0/3.d0)
       cmat0=cmat0*vtest
       a1(:)=cmat0(1,:) ; a2(:)=cmat0(2,:) ; a3(:)=cmat0(3,:)
                        endif
       read(5,'(a280)') string
       j=0
       ipos=1
       do i=1,100
       call readnext_c(string,ipos,char1(i))
!      write(6,*) trim(char1(i))
       if(char1(i)(1:1) == ' ')then
       j=i-1
       exit
                               endif
       enddo
!      write(6,*) j
       nspecies=j
       allocate(symbl(nspecies)) ; allocate(nelements(nspecies))
       do i=1,nspecies
       symbl(i)=adjustl(char1(i))
!      write(6,*) symbl(i)
       enddo
       read(5,'(a280)') string
       ipos=1
       do i=1,nspecies
       call readnext_i4(string,ipos,itmpv(i))
       nelements(i)=itmpv(i)
!      write(6,*) nelements(i)
       enddo
       natom0=0
       do i=1,nspecies
       natom0=natom0+nelements(i)
       enddo
       allocate(dir(natom0,3))
       read(5,*) ich20
       ich1='d'
       if(ich20(1:1) == 'C') ich1='c'
       if(ich20(1:1) == 'c') ich1='c'
       if(ich20(1:1) == 'K') ich1='c'
       if(ich20(1:1) == 'k') ich1='c'
!      natom0=sum(nelements)
       do i=1,natom0
       read(5,*) dir(i,1),dir(i,2),dir(i,3)
       if(ich1 == 'c')then
       dir(i,1)=dir(i,1)*scale0
       dir(i,2)=dir(i,2)*scale0
       dir(i,3)=dir(i,3)*scale0
       uec(:)=dir(i,:)
       call cart2direct(uec,vec,a1,a2,a3)
       dir(i,:)=vec(:)
                      endif
       enddo
       do i=1,natom0
       vec(:)=dir(i,:)
       if(vec(1) < 0.0d0) vec(1)=vec(1)+1.d0
       if(vec(2) < 0.0d0) vec(2)=vec(2)+1.d0
       if(vec(3) < 0.0d0) vec(3)=vec(3)+1.d0
       dir(i,:)=vec(:)
       enddo
       cmat0(1,:)=a1(:) ; cmat0(2,:)=a2(:) ; cmat0(3,:)=a3(:)
       cellvol0=(cmat0(1,2)*cmat0(2,3)-cmat0(1,3)*cmat0(2,2))*cmat0(3,1) &
               +(cmat0(1,3)*cmat0(2,1)-cmat0(1,1)*cmat0(2,3))*cmat0(3,2) &
               +(cmat0(1,1)*cmat0(2,2)-cmat0(1,2)*cmat0(2,1))*cmat0(3,3)
       cellvol0=abs(cellvol0)
!      write(6,*) cellvol0
!      original cell  ]
       zmat=cmat0
       call latmat(rlat,zmat,-1)
!      write(6,*) (rlat(j),j=1,3),(rlat(j)*180.d0/pi,j=4,6)
!      n1=2 ; n2=2 ; n3=1
       n1=1 ; n2=1 ; n3=1
!      write(6,*) n1,n2,n3
       allocate(cartext(n1*n2*n3*natom0,3),pppext(n1*n2*n3*natom0,3))
       natomext=0
       do i=1,n1
       do j=1,n2
       do k=1,n3
       uec(:)=dble(i-1)*a1(:)+dble(j-1)*a2(:)+dble(k-1)*a3(:)
       do iatom=1,natom0
       vec(1)=dir(iatom,1)*a1(1)+dir(iatom,2)*a2(1)+dir(iatom,3)*a3(1)
       vec(2)=dir(iatom,1)*a1(2)+dir(iatom,2)*a2(2)+dir(iatom,3)*a3(2)
       vec(3)=dir(iatom,1)*a1(3)+dir(iatom,2)*a2(3)+dir(iatom,3)*a3(3)
       wec(:)=vec(:)+uec(:)
       natomext=natomext+1
       cartext(natomext,:)=wec(:)      
       enddo
       enddo
       enddo
       enddo
!
       a1(:)=cmat0(1,:) ; a2(:)=cmat0(2,:) ; a3(:)=cmat0(3,:)
       cmat2(1,:)=n1*a1(:) ; cmat2(2,:)=n2*a2(:) ; cmat2(3,:)=n3*a3(:)
       vtest=(cmat2(1,2)*cmat2(2,3)-cmat2(1,3)*cmat2(2,2))*cmat2(3,1) &
            +(cmat2(1,3)*cmat2(2,1)-cmat2(1,1)*cmat2(2,3))*cmat2(3,2) &
            +(cmat2(1,1)*cmat2(2,2)-cmat2(1,2)*cmat2(2,1))*cmat2(3,3)
       vtest=abs(vtest)
!      write(6,*) vtest,vtest/cellvol0,n1,n2,n3
       a1(:)=cmat2(1,:) ; a2(:)=cmat2(2,:) ; a3(:)=cmat2(3,:)
       do iatom=1,natomext
       uec(:)=cartext(iatom,:)
       call cart2direct(uec,vec,a1,a2,a3)
       if(vec(1) < 0.0d0) vec(1)=vec(1)+1.d0
       if(vec(2) < 0.0d0) vec(2)=vec(2)+1.d0
       if(vec(3) < 0.0d0) vec(3)=vec(3)+1.d0
       pppext(iatom,:)=vec(:)
       enddo
!      write(6,*)
!      write(6,*) trim(str0)
!      write(6,*) '1.'
!      write(6,'(3f25.12)') cmat2(1,1),cmat2(1,2),cmat2(1,3)
!      write(6,'(3f25.12)') cmat2(2,1),cmat2(2,2),cmat2(2,3)
!      write(6,'(3f25.12)') cmat2(3,1),cmat2(3,2),cmat2(3,3)
!      write(6,'(10(1x,a2))') (symbl(i),i=1,nspecies)
!      write(6,'(10(i4))') (n1*n2*n3*nelements(i),i=1,nspecies)
!      write(6,*) 'direct'
!      do iatom=1,natomext
!      write(6,'(3f22.12)') pppext(iatom,1),pppext(iatom,2),pppext(iatom,3)
!      enddo
!      write(6,*)

       a1(:)=cmat0(1,:) ; a2(:)=cmat0(2,:) ; a3(:)=cmat0(3,:)
       open(91,file='../csa.in',form='formatted')
       read(91,*)
       read(91,*)
       read(91,*)
       read(91,*)
       read(91,*)
       read(91,*)
       read(91,*) tmp,tmq,height
       close(91)
       a1(3)=0.d0 ; a2(3)=0.d0 ; a3(1)=0.d0 ; a3(2)=0.d0 ; a3(3)=height

       cmat3(1,:)=a1(:) ; cmat3(2,:)=a2(:) ; cmat3(3,:)=a3(:)
       zmat=cmat3
       call latmat(rlat,zmat,-1)
!      write(6,*) (rlat(j),j=1,3),(rlat(j)*180.d0/pi,j=4,6)
       call cross3(a2,a3,b1)
       call cross3(a3,a1,b2)
       call cross3(a1,a2,b3)
       b1=2.d0*pi*b1
       b2=2.d0*pi*b2
       b3=2.d0*pi*b3
!      write(6,*) b1(1),b1(2),b1(3)
!      write(6,*) b2(1),b2(2),b2(3)
!      write(6,*) b3(1),b3(2),b3(3)
       do iatom=1,natomext
       uec(:)=cartext(iatom,:)
       call cart2direct(uec,vec,a1,a2,a3)
       if(vec(1) < 0.0d0) vec(1)=vec(1)+1.d0
       if(vec(2) < 0.0d0) vec(2)=vec(2)+1.d0
       if(vec(3) < 0.0d0) vec(3)=vec(3)+1.d0
       pppext(iatom,:)=vec(:)
       enddo
!      write(6,*)
!      write(6,*)
       write(6,*) trim(str0)
       write(6,*) '1.0'
       write(6,'(3f25.12)') cmat3(1,1),cmat3(1,2),cmat3(1,3)
       write(6,'(3f25.12)') cmat3(2,1),cmat3(2,2),cmat3(2,3)
       write(6,'(3f25.12)') cmat3(3,1),cmat3(3,2),cmat3(3,3)
       write(6,'(10(1x,a2))') (symbl(i),i=1,nspecies)
       write(6,'(10(i4))') (n1*n2*n3*nelements(i),i=1,nspecies)
       write(6,*) 'Direct'
       do iatom=1,natomext
       write(6,'(3f22.12)') pppext(iatom,1),pppext(iatom,2),pppext(iatom,3)
       enddo
       deallocate(symbl) ; deallocate(nelements)
       deallocate(dir) ; deallocate(cartext,pppext)
       stop
       end program makevertical
!234567890
!      Written by In-Ho Lee, KRISS, January 28, 2013.
       subroutine tocar(qqq,na,a1,a2,a3)
       implicit none
       integer na
       real*8 qqq(na,3),a1(3),a2(3),a3(3)
       integer j
       real*8 x,y,z
!
       do j=1,na
       x=a1(1)*qqq(j,1)+a2(1)*qqq(j,2)+a3(1)*qqq(j,3)
       y=a1(2)*qqq(j,1)+a2(2)*qqq(j,2)+a3(2)*qqq(j,3)
       z=a1(3)*qqq(j,1)+a2(3)*qqq(j,2)+a3(3)*qqq(j,3)
       qqq(j,1)=x
       qqq(j,2)=y
       qqq(j,3)=z
       enddo
       end
!234567890
!      Written by In-Ho Lee, KRISS, January 28, 2013.
       subroutine tolat(qqq,na,a1,a2,a3)
       implicit none
       integer na
       real*8 qqq(na,3),a1(3),a2(3),a3(3)
       integer j
       real*8 b(3,3),devid
       real*8 d1,d2,d3
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
       do j=1,na
       d1=b(1,1)*qqq(j,1)+b(1,2)*qqq(j,2)+b(1,3)*qqq(j,3)
       d2=b(2,1)*qqq(j,1)+b(2,2)*qqq(j,2)+b(2,3)*qqq(j,3)
       d3=b(3,1)*qqq(j,1)+b(3,2)*qqq(j,2)+b(3,3)*qqq(j,3)
       qqq(j,1)=d1
       qqq(j,2)=d2
       qqq(j,3)=d3
       enddo
       do j=1,na
       qqq(j,1)=qqq(j,1)-anint(qqq(j,1))
       qqq(j,2)=qqq(j,2)-anint(qqq(j,2))
       qqq(j,3)=qqq(j,3)-anint(qqq(j,3))
       enddo
       do j=1,na
       if(qqq(j,1) < 0.d0) qqq(j,1)=qqq(j,1)+1.d0
       if(qqq(j,2) < 0.d0) qqq(j,2)=qqq(j,2)+1.d0
       if(qqq(j,3) < 0.d0) qqq(j,3)=qqq(j,3)+1.d0
       enddo
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
!      Written by In-Ho Lee, KRISS, May 25 (2008).
       subroutine cart2direct(uec,vec,a1,a2,a3)
       implicit none
       real*8 uec(3),vec(3),a1(3),a2(3),a3(3)
       real*8 devid,b(3,3)

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
       vec(1)=b(1,1)*uec(1)+b(1,2)*uec(2)+b(1,3)*uec(3)
       vec(2)=b(2,1)*uec(1)+b(2,2)*uec(2)+b(2,3)*uec(3)
       vec(3)=b(3,1)*uec(1)+b(3,2)*uec(2)+b(3,3)*uec(3)
       vec=vec-anint(vec)
       end
!234567890
!      Written by In-Ho Lee, KRISS, January 20, 2019.
       subroutine readnext_r8(string,ipos,r8value)
       implicit none
       character(len=*), intent(in)    :: string
       integer,          intent(inout) :: ipos
       real*8,           intent(out)   :: r8value
       integer                         :: i1,i2

       i2=len_trim(string)
       if(ipos > i2)then
       ipos=0 ; r8value=0.0d0
       return
                    endif
       i1=ipos
       do
       if(string(i1:i1) /= ' ') exit
       i1=i1+1
       enddo
       read(string(i1:i2), *) r8value
       ipos=scan(string(i1:i2), ' ')
       if(ipos == 0)then
       ipos=i2+1
                    else
       ipos=ipos+i1-1
                    endif
       end subroutine readnext_r8
!234567890
!      Written by In-Ho Lee, KRISS, January 20, 2019.
       subroutine readnext_i4(string,ipos,i4value)
       implicit none
       character(len=*), intent(in)    :: string
       integer,          intent(inout) :: ipos
       integer,          intent(out)   :: i4value
       integer                         :: i1,i2

       i2=len_trim(string)
       if(ipos > i2)then
       ipos=0 ; i4value=0
       return
                    endif
       i1=ipos
       do
       if(string(i1:i1) /= ' ')exit
       i1=i1+1
       enddo
       read(string(i1:i2), *) i4value
       ipos=scan(string(i1:i2), ' ')
       if(ipos == 0)then
       ipos=i2+1
                    else
       ipos=ipos+i1-1
                    endif
       end subroutine readnext_i4
!234567890
!      Written by In-Ho Lee, KRISS, January 20, 2019.
       subroutine readnext_c(string,ipos,char1)
       implicit none
       character(len=*), intent(in)    :: string
       integer,          intent(inout) :: ipos
       character(len=*), intent(out)   :: char1
       integer                         :: i1,i2

       i2=len_trim(string)
       if(ipos > i2)then
       ipos=0 ; char1=' '
       return
                    endif
       i1=ipos
       do
       if(string(i1:i1) /= ' ') exit
       i1=i1+1
       enddo
       read(string(i1:i2), *) char1
       char1=adjustl(char1)
       char1=trim(char1)
       ipos=scan(string(i1:i2), ' ')
       if(ipos == 0)then
       ipos=i2+1
                    else
       ipos=ipos+i1-1
                    endif
       end subroutine readnext_c
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine latmatvol(rlat,wmat,volume)
       implicit none
       real*8 rlat(6),wmat(3,3),volume
       real*8 ylat(6),slat(6),zmat(3,3),tmq,tmp

       slat=rlat
       call latmat(slat,zmat,1)
       tmp=(zmat(1,2)*zmat(2,3)-zmat(1,3)*zmat(2,2))*zmat(3,1)  &
          +(zmat(1,3)*zmat(2,1)-zmat(1,1)*zmat(2,3))*zmat(3,2)  &
          +(zmat(1,1)*zmat(2,2)-zmat(1,2)*zmat(2,1))*zmat(3,3)
       tmq=volume/tmp ; tmq=tmq**(1.0d0/3.0d0)
       ylat(1)=rlat(1)*tmq ; ylat(2)=rlat(2)*tmq ; ylat(3)=rlat(3)*tmq
       ylat(4)=rlat(4) ; ylat(5)=rlat(5) ; ylat(6)=rlat(6)
       call latmat(ylat,wmat,1)
       return
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
       return
       end
!234567890
