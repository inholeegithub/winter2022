!234567890
!      ifort strings.f90 transform.f90
!      Written by In-Ho Lee, KRISS, December 30, 2015.
       implicit none
       real*8 cmat2(3,3)
       logical lfault

       write(6,*) 'Enter Bravais lattice'
       read(5,*) cmat2(1,1),cmat2(1,2),cmat2(1,3)
       read(5,*) cmat2(2,1),cmat2(2,2),cmat2(2,3)
       read(5,*) cmat2(3,1),cmat2(3,2),cmat2(3,3)
       call read_poscar_ref(lfault,cmat2)
       stop
       end
!234567890
!      Written by In-Ho Lee, KRISS, December 30, 2015.
       subroutine read_poscar_ref(lfault,cmat2)
       USE strings, ONLY : parse,value
       implicit none
       logical lfault
       real*8 cmat2(3,3)
       integer mm0,nspecies1
       real*8 cmat1(3,3),vec(3),zheight,scale1
       character*280 fname
       integer i,j,natot
       logical ldirect
       integer, allocatable :: nelements1(:)
       real*8, allocatable :: dir(:,:)
       character*2, allocatable :: symbl1(:)
       integer ios,nargs
       character*200 str1
       character*200 args(40)
       character*20 delims

       zheight=cmat2(3,3)


       lfault=.false.
       fname='POSCAR'
       mm0=0
       open(17,file=trim(fname),form='formatted')
       read(17,'(a200)',err=711,end=799) str1
!      write(6,*) len_trim(str1)
       delims=' ' ; call parse(str1,delims,args,nargs)
       call value(args(1),mm0,ios) ; if(ios /= 0) mm0=0
       write(6,*) mm0
       read(17,'(a200)',err=711,end=799) str1
       delims=' ' ; call parse(str1,delims,args,nargs)
       call value(args(1),scale1,ios)
       write(6,*) scale1
       read(17,'(a200)',err=711,end=799) str1
       delims=' ' ; call parse(str1,delims,args,nargs)
       call value(args(1),cmat1(1,1),ios)
       call value(args(2),cmat1(1,2),ios)
       call value(args(3),cmat1(1,3),ios)
       read(17,'(a200)',err=711,end=799) str1
       delims=' ' ; call parse(str1,delims,args,nargs)
       call value(args(1),cmat1(2,1),ios)
       call value(args(2),cmat1(2,2),ios)
       call value(args(3),cmat1(2,3),ios)
       read(17,'(a200)',err=711,end=799) str1
       delims=' ' ; call parse(str1,delims,args,nargs)
       call value(args(1),cmat1(3,1),ios)
       call value(args(2),cmat1(3,2),ios)
       call value(args(3),cmat1(3,3),ios)
       cmat1=cmat1*scale1
       write(6,*) cmat1(1,1),cmat1(1,2),cmat1(1,3)
       write(6,*) cmat1(2,1),cmat1(2,2),cmat1(2,3)
       write(6,*) cmat1(3,1),cmat1(3,2),cmat1(3,3)
       read(17,'(a200)',err=711,end=799) str1
       delims=' ' ; call parse(str1,delims,args,nargs)
       nspecies1=nargs
       if(.not. allocated(symbl1)) allocate(symbl1(nspecies1))
       if(.not. allocated(nelements1)) allocate(nelements1(nspecies1))
!
       do j=1,nargs
       symbl1(j)=adjustl(trim(args(j)))
!      write(6,*) symbl1(j)
       enddo
       write(6,'(20(a2,1x))') (adjustl(symbl1(j)),j=1,nspecies1)
       read(17,'(a200)',err=711,end=799) str1
       delims=' ' ; call parse(str1,delims,args,nargs)
       natot=0
       do j=1,nargs
       call value(args(j),nelements1(j),ios)
       natot=natot+nelements1(j)
!      write(6,*) nelements1(j)
       enddo
       write(6,'(20(i4))') (nelements1(j),j=1,nspecies1)
       if(.not. allocated(dir)) allocate(dir(natot,3))
       read(17,'(a200)',err=711,end=799) str1
       delims=' ' ; call parse(str1,delims,args,nargs)
       ldirect=.false.
       write(6,*) trim(args(1))
       if(args(1) == 'Direct') ldirect=.true.
       if(args(1) == 'direct') ldirect=.true.
       if(args(1) == 'D') ldirect=.true.
       if(args(1) == 'd') ldirect=.true.
       if(args(1) == 'DIR') ldirect=.true.
       if(args(1) == 'Dir') ldirect=.true.
       if(args(1) == 'dir') ldirect=.true.
       if(args(1) == 'K') ldirect=.false.
       if(args(1) == 'k') ldirect=.false.
       if(args(1) == 'C') ldirect=.false.
       if(args(1) == 'c') ldirect=.false.
       if(args(1) == 'Kar') ldirect=.false.
       if(args(1) == 'kar') ldirect=.false.
       if(args(1) == 'Car') ldirect=.false.
       if(args(1) == 'CAR') ldirect=.false.
       if(args(1) == 'KAR') ldirect=.false.
       if(args(1) == 'car') ldirect=.false.
       if(args(1) == 'Cartesian') ldirect=.false.
       if(args(1) == 'cartesian') ldirect=.false.
       if(args(1) == 'Kartesian') ldirect=.false.
       if(args(1) == 'kartesian') ldirect=.false.
!      write(6,*) ldirect
       do j=1,natot
!      read(17,'(a200)',err=711,end=799) str1
       read(17,'(a200)',err=711) str1
       delims=' ' ; call parse(str1,delims,args,nargs)
       call value(args(1),dir(j,1),ios)
       call value(args(2),dir(j,2),ios)
       call value(args(3),dir(j,3),ios)
       enddo
       if(.not. ldirect) call tolatx1(cmat1,dir,natot)
       do j=1,natot
       write(6,*) dir(j,1),dir(j,2),dir(j,3)
       enddo
!      write(6,*) 'lfault',lfault
       goto 799
  711  continue
       lfault=.true.
  799  continue
       close(17)
       call print2dlattice(cmat1)
!      write(6,*) 'lfault',lfault
!
       cmat2(1,1)=cmat1(1,1) ; cmat2(1,2)=cmat1(1,2) ; cmat2(2,1)=cmat1(2,1) ; cmat2(2,2)=cmat1(2,2)
       cmat2(3,3)=zheight
       cmat2(1,3)=0.0d0 ; cmat2(2,3)=0.0d0 ; cmat2(3,1)=0.0d0 ; cmat2(3,2)=0.0d0
       do j=1,natot
       vec(:)=dir(j,1)*cmat1(1,:)+dir(j,2)*cmat1(2,:)+dir(j,3)*cmat1(3,:)
       dir(j,:)=vec(:)
       enddo
       call tolatx1(cmat2,dir,natot)
!
       fname='POSCAR_2d'
       open(71,file=trim(fname),form='formatted')
       write(71,*) mm0
       write(71,'(a3)') '1.0'
       write(71,'(3f23.16)') cmat2(1,1),cmat2(1,2),cmat2(1,3)
       write(71,'(3f23.16)') cmat2(2,1),cmat2(2,2),cmat2(2,3)
       write(71,'(3f23.16)') cmat2(3,1),cmat2(3,2),cmat2(3,3)
       write(71,'(20(2x,a2,1x))') (symbl1(i),i=1,nspecies1)
       write(71,'(20(i4,1x))') (nelements1(i),i=1,nspecies1)
       write(71,'(a6)') "Direct"
       natot=0
       do i=1,nspecies1
       do j=1,nelements1(i)
       natot=natot+1
       write(71,'(3f20.16)') dir(natot,1),dir(natot,2),dir(natot,3)
       enddo
       enddo
       close(71)
       call print2dlattice(cmat2)
!23456789
       deallocate(symbl1,nelements1,dir)
       end
!234567890
!      Written by In-Ho Lee, KRISS, December 30, 2015.
       subroutine tolatx1(amat,dir,natot)
       implicit none
       integer natot
       real*8 amat(3,3),dir(natot,3)
       real*8 b(3,3),devid
       integer j
       real*8 a1(3),a2(3),a3(3),d1,d2,d3

       a1(:)=amat(1,:) ; a2(:)=amat(2,:) ; a3(:)=amat(3,:)
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
       do j=1,natot
       d1=b(1,1)*dir(j,1)+b(1,2)*dir(j,2)+b(1,3)*dir(j,3)
       d2=b(2,1)*dir(j,1)+b(2,2)*dir(j,2)+b(2,3)*dir(j,3)
       d3=b(3,1)*dir(j,1)+b(3,2)*dir(j,2)+b(3,3)*dir(j,3)
       dir(j,1)=d1
       dir(j,2)=d2
       dir(j,3)=d3
       enddo
       do j=1,natot
       dir(j,1)=dir(j,1)-anint(dir(j,1))
       dir(j,2)=dir(j,2)-anint(dir(j,2))
       dir(j,3)=dir(j,3)-anint(dir(j,3))
       enddo
       do j=1,natot
       if(dir(j,1) < 0.d0) dir(j,1)=dir(j,1)+1.d0
       if(dir(j,2) < 0.d0) dir(j,2)=dir(j,2)+1.d0
       if(dir(j,3) < 0.d0) dir(j,3)=dir(j,3)+1.d0
       enddo
       end 
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine print2dlattice(cmat)
       implicit none
       real*8 cmat(3,3)
       real*8 a1(3),a2(3),a3(3),areap
       real*8 ra,rb,rc,alpha,beta,gama,cosinea,cosineb,cosinec,tmp,pi

       ra=sqrt(cmat(1,1)**2+cmat(1,2)**2+cmat(1,3)**2)
       rb=sqrt(cmat(2,1)**2+cmat(2,2)**2+cmat(2,3)**2)
       rc=sqrt(cmat(3,1)**2+cmat(3,2)**2+cmat(3,3)**2)
       cosinea=(cmat(2,1)*cmat(3,1)+cmat(2,2)*cmat(3,2)+cmat(2,3)*cmat(3,3))/rb/rc
       cosineb=(cmat(1,1)*cmat(3,1)+cmat(1,2)*cmat(3,2)+cmat(1,3)*cmat(3,3))/rc/ra
       cosinec=(cmat(1,1)*cmat(2,1)+cmat(1,2)*cmat(2,2)+cmat(1,3)*cmat(2,3))/ra/rb  
       pi=4.0d0*atan(1.0d0)
       tmp=180.0d0/pi
       alpha=tmp*acos(cosinea) ; beta=tmp*acos(cosineb) ; gama=tmp*acos(cosinec)
       write(6,'(6f16.5)') ra,rb,rc,alpha,beta,gama
       a1(:)=cmat(1,:) ; a2(:)=cmat(2,:) ; a3(:)=cmat(3,:)
       tmp=(a1(2)*a2(3)-a1(3)*a2(2))*a3(1) &
          +(a1(3)*a2(1)-a1(1)*a2(3))*a3(2) &
          +(a1(1)*a2(2)-a1(2)*a2(1))*a3(3)
       tmp=abs(tmp)
       areap=tmp/abs(cmat(3,3)) 
       write(6,'(f20.6,1x,a5,2x,f20.6,1x,a5)') tmp, '(A^3)', areap, '(A^2)'
       end
!234567890
