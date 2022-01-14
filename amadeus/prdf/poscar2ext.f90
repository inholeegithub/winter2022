!234567890
!      Written by In-Ho Lee, KRISS, April 13, 2016.
!      awk '{if(NR ==1) {print "12.34", $0} else { print }}' POSCAR
!      ifort -o poscar2ext.x strings.f90  poscar2ext.f90 
!      ./poscar2ext.x < POSCAR
       program poscar2ext
       USE strings, ONLY : parse,value
       implicit none
       integer nspecies,ncext(3)
       real*8 scale0,a1(3),a2(3),a3(3),tmp,aa1(3),aa2(3),aa3(3),rmax,vtest,amatrix(3,3)
       integer, allocatable :: nelements(:)
       real*8, allocatable :: dir(:,:),car(:,:),dirext(:,:),carext(:,:)
       character*2, allocatable :: symbl(:)
       integer i,j,k,ie,je,ke,k0,kk0,ks
       integer ios,nargs
       character*200 str1,str2
       character*200 args(40)
       character*20 delims
       logical lfault

       rmax=10.5d0
!      write(6,*) rmax,' rmax'
       lfault=.false.
       read(5,'(a200)',err=911,end=999) str2
       delims=' '
       call parse(str2,delims,args,nargs)
       if(nargs > 0)then
       call value(args(1),rmax,ios)
       if(ios /= 0)then
       rmax=10.5d0
       write(6,*) 'default rmax',rmax
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
       read(5,'(a200)',err=911,end=999) str1
       delims=' '
       call parse(str1,delims,args,nargs)
       nspecies=nargs
       allocate(nelements(nspecies)) ; allocate(symbl(nspecies))
       do i=1,nspecies
       symbl(i)=trim(adjustl(args(i)))
       enddo
       read(5,'(a200)',err=911,end=999) str1
       delims=' '
       call parse(str1,delims,args,nargs)
       do i=1,nspecies
       call value(args(i),nelements(i),ios)
!      write(6,*) nelements(i)
       enddo
       k=sum(nelements)
       ks=k
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
       car(j,1)=dir(j,1)*a1(1)+dir(j,2)*a2(1)+dir(j,3)*a3(1)
       car(j,2)=dir(j,1)*a1(2)+dir(j,2)*a2(2)+dir(j,3)*a3(2)
       car(j,3)=dir(j,1)*a1(3)+dir(j,2)*a2(3)+dir(j,3)*a3(3)
       enddo
                                                                             else
       do j=1,ks
       read(5,*) car(j,1),car(j,2),car(j,3)
       enddo
       dir=car
       call tolaty(dir,a1,a2,a3,ks)
                                                                             endif
                    endif
  911  continue
       lfault=.true.
  999  continue
       call get_extension(a1,a2,a3,rmax,ncext)
       aa1=a1*dble(ncext(1)) ; aa2=a2*dble(ncext(2)) ; aa3=a3*dble(ncext(3))
       k=sum(nelements)*(ncext(1)*ncext(2)*ncext(3))
       allocate(dirext(k,3)) ; allocate(carext(k,3))
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
       enddo
       enddo
       enddo
       enddo
       if(k0 /= nelements(i))then
       write(6,*) 'error, k0'
                             stop
                             endif
       enddo
       if(k /= ks*ncext(1)*ncext(2)*ncext(3))then
       write(6,*) 'error ks'
                                             stop
                                             endif
       do i=1,sum(nelements)*ncext(1)*ncext(2)*ncext(3)
       carext(i,1)=dirext(i,1)*aa1(1)+dirext(i,2)*aa2(1)+dirext(i,3)*aa3(1)
       carext(i,2)=dirext(i,1)*aa1(2)+dirext(i,2)*aa2(2)+dirext(i,3)*aa3(2)
       carext(i,3)=dirext(i,1)*aa1(3)+dirext(i,2)*aa2(3)+dirext(i,3)*aa3(3)
       enddo
!      dirext=carext
!      call tolaty(dirext,aa1,aa2,aa3,k)
       write(6,'(3i4,1x,f18.8)') ncext,rmax
       write(6,*) '1.0'
       write(6,'(3f22.12)') aa1(1),aa1(2),aa1(3)
       write(6,'(3f22.12)') aa2(1),aa2(2),aa2(3)
       write(6,'(3f22.12)') aa3(1),aa3(2),aa3(3)
       write(6,'(10(1x,a2))') (symbl(i),i=1,nspecies)
       write(6,'(10(i6))') (ncext(1)*ncext(2)*ncext(3)*nelements(i),i=1,nspecies)
       write(6,*) 'direct'
       k=0
       do i=1,nspecies
       do j=1,nelements(i)*ncext(1)*ncext(2)*ncext(3)
       k=k+1
       write(6,'(3f22.12)') dirext(k,1),dirext(k,2),dirext(k,3)
       enddo
       enddo
       deallocate(dirext,carext) ; deallocate(symbl,nelements)
       deallocate(dir,car)
       stop
       end program poscar2ext
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
