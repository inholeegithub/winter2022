!234567890
!      ifort -o poscar2xsf.x  string2.f90 poscar2xsf.f90
!      Written by In-Ho Lee, KRISS, March 27, 2016
       program poscar2xsf
       USE strings, ONLY : parse,value
       implicit none
       integer nspecies
       real*8 scale0,a1(3),a2(3),a3(3),tmp,amatrix(3,3),vtest
       integer, allocatable :: natoms(:)
       real*8, allocatable :: dir(:,:),car(:,:)
       character*2, allocatable :: symbl(:)
       integer i,j,k
       integer ios,nargs
       character*200 str1
       character*200 args(40)
       character*20 delims
       logical lfault

       lfault=.false.
       read(5,*) tmp
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
       allocate(natoms(nspecies)) ; allocate(symbl(nspecies))
       do i=1,nspecies
       symbl(i)=trim(adjustl(args(i)))
       enddo
       read(5,'(a200)',err=911,end=999) str1
       delims=' '
       call parse(str1,delims,args,nargs)
       do i=1,nspecies
       call value(args(i),natoms(i),ios)
       enddo
       k=sum(natoms)
       allocate(dir(k,3)) ; allocate(car(k,3))
       read(5,'(a200)',err=911,end=999) str1
       delims=' '
       call parse(str1,delims,args,nargs)
       if(nargs > 0)then
       if(args(1) == 'DIR' .or. args(1) == 'dir' .or. args(1) == 'D' .or. args(1) == 'd' .or. &
          args(1) == 'direct' .or. args(1) == 'Direct' .or. args(1) == 'Dir')then
       do j=1,k
       read(5,*) dir(j,1),dir(j,2),dir(j,3)
       car(j,1)=dir(j,1)*a1(1)+dir(j,2)*a2(1)+dir(j,3)*a3(1)
       car(j,2)=dir(j,1)*a1(2)+dir(j,2)*a2(2)+dir(j,3)*a3(2)
       car(j,3)=dir(j,1)*a1(3)+dir(j,2)*a2(3)+dir(j,3)*a3(3)
       enddo
                                                                             else
       do j=1,k
       read(5,*) car(j,1),car(j,2),car(j,3)
       enddo
                                                                             endif
                    endif
  911  continue
       lfault=.true.
  999  continue
!      aenet-xsf format
       write(6,'(a17,1x,f22.8,1x,a2)') '# total energy = ',tmp,'eV'
       write(6,*) 'CRYSTAL'
       write(6,*) 'PRIMVEC'
       write(6,'(3f22.12)') a1(1),a1(2),a1(3)
       write(6,'(3f22.12)') a2(1),a2(2),a2(3)
       write(6,'(3f22.12)') a3(1),a3(2),a3(3)
       write(6,*) 'PRIMCOORD'
       write(6,'(i3,1x,i1)') k, 1
       k=0
       do i=1,nspecies
       do j=1,natoms(i)
       k=k+1
       write(6,'(1x,a2,1x,3f22.12)') symbl(i),car(k,1),car(k,2),car(k,3)
       enddo
       enddo
       deallocate(symbl,natoms)
       deallocate(dir,car)
       stop
       end program poscar2xsf
