!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       program psocar2xyz
       USE strings, ONLY : parse,value
       implicit none
       integer nspecies,natom
       real*8 a1(3),a2(3),a3(3),scale0,aa1(3),aa2(3),aa3(3),amatrix(3,3),vtest
       integer, allocatable :: nelements(:),itype(:)
       real*8, allocatable :: qqq(:),extqqq(:)
       integer, allocatable :: itype_ext(:)
       character*2, allocatable :: symbl(:)
       character*1 ch1
       integer n1,n2,n3,i,j,k,m,nb,natom_ext
       character*200 str1
       integer ios,nargs
       character*200 args(40)
       character*20 delims

       open(1,file='CONTCAR',form='formatted')
       read(1,*) str1
       read(1,*) scale0
       read(1,*) a1(1),a1(2),a1(3)
       read(1,*) a2(1),a2(2),a2(3)
       read(1,*) a3(1),a3(2),a3(3)
       if(scale0 < 0.d0)then
       amatrix(1,:)=a1(:) ; amatrix(2,:)=a2(:) ; amatrix(3,:)=a3(:)
       vtest=(amatrix(1,2)*amatrix(2,3)-amatrix(1,3)*amatrix(2,2))*amatrix(3,1) &
            +(amatrix(1,3)*amatrix(2,1)-amatrix(1,1)*amatrix(2,3))*amatrix(3,2) &
            +(amatrix(1,1)*amatrix(2,2)-amatrix(1,2)*amatrix(2,1))*amatrix(3,3)
       vtest=abs(vtest)
       vtest=abs(scale0)/vtest ; vtest=vtest**(1.d0/3.d0)
       amatrix=amatrix*vtest
       a1(:)=amatrix(1,:) ; a2(:)=amatrix(2,:) ; a3(:)=amatrix(3,:)
                        endif
       if(scale0 > 0.d0)then
       a1=a1*scale0 ; a2=a2*scale0 ; a3=a3*scale0
                        endif
!      read(1,*) (symbl(i),i=1,nspecies)
!      read(1,*) (nelements(i),i=1,nspecies)
       read(1,'(a200)') str1
       delims=' '
       call parse(str1,delims,args,nargs)
       nspecies=nargs
       allocate(symbl(nspecies))
       allocate(nelements(nspecies))
       do j=1,nargs
       symbl(j)=trim(adjustl(args(j)))
       enddo
!print*, nargs, (symbl(j),j=1,nargs)
       read(1,'(a200)') str1
       delims=' '
       call parse(str1,delims,args,nargs)
       do j=1,nargs
       call value(args(j),i,ios)
       nelements(j)=i
       enddo
!print*, nargs, (nelements(j),j=1,nargs)
!
       natom=0
       do i=1,nspecies
       do j=1,nelements(i)
       natom=natom+1
       enddo
       enddo
       allocate(itype(natom)) ; allocate(qqq(3*natom))
       natom=0
       do i=1,nspecies
       do j=1,nelements(i)
       natom=natom+1
       itype(natom)=i
       enddo
       enddo
       read(1,*) ch1
       if(ch1 =='D' .or. ch1=='d')then
       do i=1,natom
       read(1,*) qqq(3*(i-1)+1),qqq(3*(i-1)+2),qqq(3*(i-1)+3)
       enddo
!
       n1=2 ; n2=2 ; n3=2
       j=natom*(n1*n2*n3)
       allocate(extqqq(3*j)) ; allocate(itype_ext(3*j))
       aa1=a1*n1 ; aa2=a2*n2 ; aa3=a3*n3
       nb=0
       do m=1,natom
       do i=0,n1-1
       do j=0,n2-1
       do k=0,n3-1
       nb=nb+1
       extqqq(3*(nb-1)+1)=qqq(3*(m-1)+1)+float(i)
       extqqq(3*(nb-1)+2)=qqq(3*(m-1)+2)+float(j)
       extqqq(3*(nb-1)+3)=qqq(3*(m-1)+3)+float(k)
       itype_ext(nb)=itype(m)
       enddo
       enddo
       enddo
       enddo
       do i=1,nb
       extqqq(3*(i-1)+1)=extqqq(3*(i-1)+1)/float(n1)
       extqqq(3*(i-1)+2)=extqqq(3*(i-1)+2)/float(n2)
       extqqq(3*(i-1)+3)=extqqq(3*(i-1)+3)/float(n3)
       enddo
       natom_ext=nb
!
       call tocarx(natom,a1,a2,a3,qqq)
       write(6,*) natom
       write(6,*) trim(str1)
       do i=1,natom
       write(6,'(a2,2x,3f24.16)') symbl(itype(i)),qqq(3*(i-1)+1),qqq(3*(i-1)+2),qqq(3*(i-1)+3)
       enddo
!
       call tocarx(natom_ext,aa1,aa2,aa3,extqqq)
       write(6,*) natom_ext
       write(6,*) trim(str1), n1,n2,n3
       do i=1,natom_ext
       write(6,'(a2,2x,3f24.16)') symbl(itype_ext(i)),extqqq(3*(i-1)+1),extqqq(3*(i-1)+2),extqqq(3*(i-1)+3)
       enddo
                                  endif
       close(1)

       deallocate(itype_ext) ; deallocate(extqqq)
       deallocate(qqq) ; deallocate(itype)
       deallocate(symbl) ; deallocate(nelements)
       stop
       end
!234567890
       subroutine tocarx(natom,a1,a2,a3,qqq)
!      Written by In-Ho Lee, KRISS, January 28, 2013.
       implicit none
       integer natom
       real*8 qqq(3*natom),a1(3),a2(3),a3(3)
       integer j,i
       real*8 x,y,z

       do j=1,natom
       x=(a1(1)*qqq(3*(j-1)+1)+a2(1)*qqq(3*(j-1)+2)+a3(1)*qqq(3*(j-1)+3))
       y=(a1(2)*qqq(3*(j-1)+1)+a2(2)*qqq(3*(j-1)+2)+a3(2)*qqq(3*(j-1)+3))
       z=(a1(3)*qqq(3*(j-1)+1)+a2(3)*qqq(3*(j-1)+2)+a3(3)*qqq(3*(j-1)+3))
       qqq(3*(j-1)+1)=x
       qqq(3*(j-1)+2)=y
       qqq(3*(j-1)+3)=z
       enddo
       end subroutine tocarx
!234567890
       subroutine tolatx(natom,a1,a2,a3,qqq)
!      Written by In-Ho Lee, KRISS, January 28, 2013.
       implicit none
       integer natom
       real*8 qqq(3*natom),a1(3),a2(3),a3(3)
       real*8 b(3,3),devid
       real*8 d1,d2,d3
       integer j,i

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
       do j=1,natom
       d1=(b(1,1)*qqq(3*(j-1)+1)+b(1,2)*qqq(3*(j-1)+2)+b(1,3)*qqq(3*(j-1)+3))
       d2=(b(2,1)*qqq(3*(j-1)+1)+b(2,2)*qqq(3*(j-1)+2)+b(2,3)*qqq(3*(j-1)+3))
       d3=(b(3,1)*qqq(3*(j-1)+1)+b(3,2)*qqq(3*(j-1)+2)+b(3,3)*qqq(3*(j-1)+3))
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
       if(qqq(3*(j-1)+1) <0.d0) qqq(3*(j-1)+1)=qqq(3*(j-1)+1)+1.d0
       if(qqq(3*(j-1)+2) <0.d0) qqq(3*(j-1)+2)=qqq(3*(j-1)+2)+1.d0
       if(qqq(3*(j-1)+3) <0.d0) qqq(3*(j-1)+3)=qqq(3*(j-1)+3)+1.d0
       enddo
       end subroutine tolatx
!234567890
