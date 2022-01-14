!234567890
!      awk '{if(NR ==1) {print "12.34", $0} else { print }}' POSCAR
!      ifort -o poscar2xyz.x strings.f90  poscar2xyz.f90 
!      ./poscar2xyz.x < POSCAR
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       program psocar2xyz
       USE strings, ONLY : parse,value
       implicit none
       integer nspecies,natot
       real*8 a1(3),a2(3),a3(3),scale0,aa1(3),aa2(3),aa3(3),vtest,amatrix(3,3)
       real*8 rmax
       integer, allocatable :: nelements(:),itype(:)
       real*8, allocatable :: qqq(:),extqqq(:),car(:,:)
       integer, allocatable :: itype_ext(:)
       character*2, allocatable :: symbl(:)
       character*1 ch1
       integer n1,n2,n3,i,j,k,m,nb,natot_ext,ncext(3)
       character*200 str1,str2
       integer ios,nargs
       character*200 args(40)
       character*20 delims

       rmax=10.5d0
       read(5,*) str2
       delims=' '
       call parse(str2,delims,args,nargs)
       if(nargs > 0)then
       call value(args(1),rmax,ios)
       if(ios /= 0)then
       rmax=10.5d0
                   endif
                    endif
       read(5,*) scale0
       read(5,*) a1(1),a1(2),a1(3)
       read(5,*) a2(1),a2(2),a2(3)
       read(5,*) a3(1),a3(2),a3(3)
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
       read(5,'(a200)') str1
       delims=' '
       call parse(str1,delims,args,nargs)
       nspecies=nargs
       allocate(symbl(nspecies))
       allocate(nelements(nspecies))
       do j=1,nargs
       symbl(j)=trim(adjustl(args(j)))
       enddo
       read(5,'(a200)') str1
       delims=' '
       call parse(str1,delims,args,nargs)
       do j=1,nargs
       call value(args(j),i,ios)
       nelements(j)=i
       enddo
       natot=0
       do i=1,nspecies
       do j=1,nelements(i)
       natot=natot+1
       enddo
       enddo
       allocate(itype(natot)) ; allocate(qqq(3*natot))
       natot=0
       do i=1,nspecies
       do j=1,nelements(i)
       natot=natot+1
       itype(natot)=i
       enddo
       enddo
       read(5,*) ch1
       if(ch1 =='D' .or. ch1=='d')then
       do i=1,natot
       read(5,*) qqq(3*(i-1)+1),qqq(3*(i-1)+2),qqq(3*(i-1)+3)
       enddo
                                  else
       allocate(car(natot,3))
       do i=1,natot
       read(5,*) car(i,1),car(i,2),car(i,3)
       enddo
       call tolaty(car,a1,a2,a3,natot)
       do i=1,natot
       qqq(3*(i-1)+1)=car(i,1)
       qqq(3*(i-1)+2)=car(i,2)
       qqq(3*(i-1)+3)=car(i,3)
       enddo
       deallocate(car)
                                  endif
       do i=1,natot
       qqq(3*(i-1)+1)=qqq(3*(i-1)+1)-anint(qqq(3*(i-1)+1))
       qqq(3*(i-1)+2)=qqq(3*(i-1)+2)-anint(qqq(3*(i-1)+2))
       qqq(3*(i-1)+3)=qqq(3*(i-1)+3)-anint(qqq(3*(i-1)+3))
       if(qqq(3*(i-1)+1) < 0.d0) qqq(3*(i-1)+1)=qqq(3*(i-1)+1)+1.d0
       if(qqq(3*(i-1)+2) < 0.d0) qqq(3*(i-1)+2)=qqq(3*(i-1)+2)+1.d0
       if(qqq(3*(i-1)+3) < 0.d0) qqq(3*(i-1)+3)=qqq(3*(i-1)+3)+1.d0
       enddo
!
       call get_extension(a1,a2,a3,rmax,ncext)
       n1=ncext(1) ; n2=ncext(2) ; n3=ncext(3)
       j=natot*(n1*n2*n3)
       allocate(extqqq(3*j)) ; allocate(itype_ext(3*j))
       aa1=a1*n1 ; aa2=a2*n2 ; aa3=a3*n3
       nb=0
       do m=1,natot
       do i=0,n1-1
       do j=0,n2-1
       do k=0,n3-1
       nb=nb+1
       extqqq(3*(nb-1)+1)=qqq(3*(m-1)+1)+dble(i)
       extqqq(3*(nb-1)+2)=qqq(3*(m-1)+2)+dble(j)
       extqqq(3*(nb-1)+3)=qqq(3*(m-1)+3)+dble(k)
       itype_ext(nb)=itype(m)
       enddo
       enddo
       enddo
       enddo
       do i=1,nb
       extqqq(3*(i-1)+1)=extqqq(3*(i-1)+1)/dble(n1)
       extqqq(3*(i-1)+2)=extqqq(3*(i-1)+2)/dble(n2)
       extqqq(3*(i-1)+3)=extqqq(3*(i-1)+3)/dble(n3)
       enddo
       natot_ext=nb
!
!      call tocarx(natot,a1,a2,a3,qqq)
!      write(6,*) natot
!      write(6,*) trim(str2)
!      do i=1,natot
!      write(6,'(a2,2x,3f24.16)') symbl(itype(i)),qqq(3*(i-1)+1),qqq(3*(i-1)+2),qqq(3*(i-1)+3)
!      enddo
!
       call tocarx(natot_ext,aa1,aa2,aa3,extqqq)
       write(6,'(i7)') natot_ext
!      write(6,*) trim(str2), n1,n2,n3
       write(6,'(3i5,1x,f18.8)') n1,n2,n3,rmax
       do i=1,natot_ext
       write(6,'(a2,2x,3f24.16)') symbl(itype_ext(i)),extqqq(3*(i-1)+1),extqqq(3*(i-1)+2),extqqq(3*(i-1)+3)
       enddo

       deallocate(itype_ext) ; deallocate(extqqq)
       deallocate(qqq) ; deallocate(itype)
       deallocate(symbl) ; deallocate(nelements)
       stop
       end program psocar2xyz
!234567890
!      Written by In-Ho Lee, KRISS, April 13, 2016.
       subroutine get_extension(a1,a2,a3,rmax,ncext)
       implicit none
       real*8 a1(3),a2(3),a3(3),rmax
       integer ncext(3)
       real*8 v(3),h(3)

       call cross3(a1,a2,v) ; v=v/sqrt(sum(v*v)) ; h(3)=abs(sum(v*a3))
       call cross3(a3,a1,v) ; v=v/sqrt(sum(v*v)) ; h(2)=abs(sum(v*a2))
       call cross3(a2,a3,v) ; v=v/sqrt(sum(v*v)) ; h(1)=abs(sum(v*a1))
       v=rmax/h+0.5d0
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
       return
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
!      Written by In-Ho Lee, KRISS, January 28, 2013.
       subroutine tocarx(natom,a1,a2,a3,qqq)
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
       end
!234567890
!      Written by In-Ho Lee, KRISS, January 28, 2013.
       subroutine tolatx(natom,a1,a2,a3,qqq)
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
       end
!234567890
