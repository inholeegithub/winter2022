!234567890
!      Written by In-Ho Lee, KRISS, January 20, 2019.
       program new_bravais_lattice
       implicit none
       integer n1,n2,n3
       real*8 a1(3),a2(3),a3(3),cmat0(3,3),cmat2(3,3),cellvol0,vtest,vtest3
       real*8 cmat3(3,3),tmp,scale0,uec(3),vec(3),wec(3),pi,tmpv(100)
       real*8 dir(1000,3),cart(1000,3),ppp(1000,3)
       real*8 cart3(1000,3),ppp3(1000,3)
       integer iatom,natom0,i,j,k,ipos,natomext
       integer nspecies,itmpv(100)
       integer, allocatable :: nelements(:)
       character*2, allocatable :: symbl(:)
       character*20 ich20
       character*1 ich1
       character*280 string,char1(100)




!      original cell  [
       read(5,*)
       read(5,*) scale0
       read(5,*)  a1(1),a1(2),a1(3)
       read(5,*)  a2(1),a2(2),a2(3)
       read(5,*)  a3(1),a3(2),a3(3)
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
       call readnext_c(string, ipos, char1(i))
!      write(6,*) trim(char1(i))
       if(char1(i)(1:1) ==  ' ')then
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
       call readnext_i4(string, ipos, itmpv(i))
       nelements(i)=itmpv(i)
!      write(6,*) nelements(i)
       enddo
       read(5,*) ich20
       ich1='d'
       if(ich20(1:1) == 'C') ich1='c'
       if(ich20(1:1) == 'c') ich1='c'
       if(ich20(1:1) == 'K') ich1='c'
       if(ich20(1:1) == 'k') ich1='c'
       natom0=sum(nelements)

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
       write(6,*) cellvol0
!      original cell  ]





       n1=2 ; n2=2 ; n3=1
       n1=1 ; n2=1 ; n3=1
       write(6,*) n1,n2,n3




       pi=4.d0*atan(1.d0)
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
       cart(natomext,:)=wec(:)
       cart3(natomext,:)=wec(:)
       enddo
       enddo
       enddo
       enddo




       a1(:)=cmat0(1,:) ; a2(:)=cmat0(2,:) ; a3(:)=cmat0(3,:)
       cmat2(1,:)=n1*a1(:) ; cmat2(2,:)=n2*a2(:) ; cmat2(3,:)=n3*a3(:)
       vtest=(cmat2(1,2)*cmat2(2,3)-cmat2(1,3)*cmat2(2,2))*cmat2(3,1) &
            +(cmat2(1,3)*cmat2(2,1)-cmat2(1,1)*cmat2(2,3))*cmat2(3,2) &
            +(cmat2(1,1)*cmat2(2,2)-cmat2(1,2)*cmat2(2,1))*cmat2(3,3)
       vtest=abs(vtest)
       write(6,*) vtest,vtest/cellvol0,n1,n2,n3




       a1(:)=cmat2(1,:) ; a2(:)=cmat2(2,:) ; a3(:)=cmat2(3,:)
       do iatom=1,natomext
       uec(:)=cart(iatom,:)
       call cart2direct(uec,vec,a1,a2,a3)
       if(vec(1) < 0.0d0) vec(1)=vec(1)+1.d0
       if(vec(2) < 0.0d0) vec(2)=vec(2)+1.d0
       if(vec(3) < 0.0d0) vec(3)=vec(3)+1.d0
       ppp(iatom,:)=vec(:)
       enddo
       write(6,*)
       write(6,*) 'ext-simple'
       write(6,*) '1.'
       write(6,'(3f25.12)') cmat2(1,1),cmat2(1,2),cmat2(1,3)
       write(6,'(3f25.12)') cmat2(2,1),cmat2(2,2),cmat2(2,3)
       write(6,'(3f25.12)') cmat2(3,1),cmat2(3,2),cmat2(3,3)
       write(6,'(10(1x,a2))')  (symbl(i),i=1,nspecies)
       write(6,'(10(i4))')  (n1*n2*n3*nelements(i),i=1,nspecies)
       write(6,*) 'direct'
       do iatom=1,natomext
       write(6,'(3f22.12)') ppp(iatom,1),ppp(iatom,2),ppp(iatom,3)
       enddo
       write(6,*)










       a1(:)=cmat0(1,:) ; a2(:)=cmat0(2,:) ; a3(:)=cmat0(3,:)
       a3(3)=20.d0





       cmat3(1,:)=a1(:) ; cmat3(2,:)=a2(:) ; cmat3(3,:)=a3(:)
       do iatom=1,natomext
       uec(:)=cart3(iatom,:)
       call cart2direct(uec,vec,a1,a2,a3)
       if(vec(1) < 0.0d0) vec(1)=vec(1)+1.d0
       if(vec(2) < 0.0d0) vec(2)=vec(2)+1.d0
       if(vec(3) < 0.0d0) vec(3)=vec(3)+1.d0
       ppp3(iatom,:)=vec(:)
       enddo




       write(6,*)
       write(6,*)
       write(6,*) 'ext-new-b'
       write(6,*) '1.'
       write(6,'(3f25.12)') cmat3(1,1),cmat3(1,2),cmat3(1,3)
       write(6,'(3f25.12)') cmat3(2,1),cmat3(2,2),cmat3(2,3)
       write(6,'(3f25.12)') cmat3(3,1),cmat3(3,2),cmat3(3,3)
       write(6,'(10(1x,a2))')  (symbl(i),i=1,nspecies)
       write(6,'(10(i4))')  (n1*n2*n3*nelements(i),i=1,nspecies)
       write(6,*) 'direct'
       do iatom=1,natomext
       write(6,'(3f22.12)') ppp3(iatom,1),ppp3(iatom,2),ppp3(iatom,3)
       enddo




       deallocate(symbl) ; deallocate(nelements)
       stop
       end program new_bravais_lattice
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

