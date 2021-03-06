!
!      ######   ######     ###            ##     ##    ###     ######  ########
!     ##    ## ##    ##   ## ##           ##     ##   ## ##   ##    ## ##     ##
!     ##       ##        ##   ##          ##     ##  ##   ##  ##       ##     ##
!     ##        ######  ##     ## ####### ##     ## ##     ##  ######  ########
!     ##             ## #########          ##   ##  #########       ## ##     
!     ##    ## ##    ## ##     ##           ## ##   ##     ## ##    ## ##    
!      ######   ######  ##     ##            ###    ##     ##  ######  ##   
!           
!              #     #     #     #     ######   #######  #     #   ##### 
!             # #    ##   ##    # #    #     #  #        #     #  #     #
!            #   #   # # # #   #   #   #     #  #        #     #  #  
!           #     #  #  #  #  #     #  #     #  #####    #     #   ##### 
!           #######  #     #  #######  #     #  #        #     #        #
!           #     #  #     #  #     #  #     #  #        #     #  #     #
!           #     #  #     #  #     #  ######   #######   #####    #####
!                       Ab initio MAterials DEsign Using cSa           
!
!       ifort strings.f90 rmarin.f90 timestamp.f90 rlx_design.f90
!234567890
!      Written by In-Ho Lee, KRISS, December 30, 2015.
       module celldesign
       implicit none
       public
       save
       integer mm0,nspecies1
       real*8 cmat1(3,3),pp1(3,3),z0frac,z1frac,zheight
       real*8 cellvol0,voltol,refa1(3),refa2(3),refa3(3)
       integer, allocatable :: nelements1(:)
       real*8, allocatable :: dir(:,:),fxyz(:,:),sigmamatrix1(:,:),dirref(:,:)
       character*2, allocatable :: symbl1(:)
       logical ldesign,lrlxion,lrlxcell,lfirst,llast
       logical lvcs
       logical, allocatable :: lmove(:,:)
      
       contains

       subroutine setfree_celldesign()
       implicit none

       if(allocated(sigmamatrix1)) deallocate(sigmamatrix1)
       if(allocated(nelements1)) deallocate(nelements1)
       if(allocated(symbl1)) deallocate(symbl1)
       if(allocated(fxyz)) deallocate(fxyz)
       if(allocated(dir)) deallocate(dir)
       if(allocated(lmove)) deallocate(lmove)
       if(allocated(dirref)) deallocate(dirref)
       end subroutine setfree_celldesign

       end module celldesign
!234567890
!      Written by In-Ho Lee, KRISS, December 30, 2015.
       program rlx_design
       USE celldesign, ONLY : lfirst,llast,setfree_celldesign
       implicit none
       logical lfault1,lfault2,lfault_stdout,lexist
       integer ij,kl
       real*8 tmp,tmq

       call csa_vasp_banner1()
       call timestamp()
       call init_seed()
       call random_number(tmp) ; ij=tmp*31328.d0+1.d0
       call random_number(tmq) ; kl=tmq*30081.d0+1.d0
       call rmarin(ij,kl)
       inquire(file='KPOINTS',exist=lexist)
       if(lexist)then
       open(11,file='KPOINTS') ; close(11,status='delete')
                 endif
       inquire(file='LEXIT',exist=lexist)
       if(lexist)then
       open(11,file='LEXIT') ; close(11,status='delete')
                 endif
       llast=.false. 
       lfirst=.true. 
       inquire(file='POSCAR_design',exist=lexist)
       if(lexist)then
       lfirst=.false.
                 endif
       if(lfirst)then
       inquire(file='OUTCAR',exist=lexist)
       if(lexist)then
       open(11,file='OUTCAR') ; close(11,status='delete')
                 endif
       inquire(file='stdout.log',exist=lexist)
       if(lexist)then
       open(11,file='stdout.log') ; close(11,status='delete')
                 endif
       write(6,*) 'cell design, rectification/vacuum, no relaxation'
!      discard previous calculation by removing OUTCAR and stdout.log
                 else
       write(6,*) 'iterative relaxation for the designed cell'
!      utilize previous calculation by readding OUTCAR and stdout.log
                 endif
!      POSCAR file must be present. in addition ../csa.in and/or ../POSCAR_ref
       call read_poscar_design(lfault1)
       call read_stdout_design(lfault_stdout)
       call read_outcar_design(lfault2)
!
       call setfree_celldesign()
       call timestamp()
       stop
       end program rlx_design
!234567890
!      Written by In-Ho Lee, KRISS, December 30, 2015.
       subroutine write_lexit()
       USE celldesign, ONLY : lfirst
       implicit none
       real*8 dtinsec

       if(.not. lfirst)then
       open(11,file='LEXIT',form='formatted')
       write(11,*) .true.
       close(11)
       dtinsec=0.5d0
       call f90sleep(dtinsec)
                       endif
       end
!234567890
!      Written by In-Ho Lee, KRISS, December 30, 2015.
       subroutine read_poscar_design(lfault)
       USE celldesign, ONLY : mm0,cmat1,symbl1,nelements1,nspecies1,dir,fxyz,sigmamatrix1
       USE celldesign, ONLY : cellvol0,voltol,lvcs,lmove,refa1,refa2,refa3,dirref,z0frac,z1frac,zheight
       USE strings, ONLY : parse,value
       implicit none
       logical lfault
       character*280 fname
       character*1 ch1
       integer i,j,natot
       real*8 scale1,tmp,vtest,cmat0(3,3)
       logical ldirect,lexist,lselective
       integer ios,nargs
       character*280 str1
       character*280 args(40)
       character*20 delims

       lfault=.false.
       fname='POSCAR'
       inquire(file=trim(fname),exist=lexist)
       if(.not. lexist)then
       write(6,*) 'POSCAR is not present'
       lfault=.true.
       if(lfault) call write_lexit()
       return
                       endif
       open(17,file=trim(fname),form='formatted')
       read(17,'(a280)',err=711,end=799) str1
!      write(6,*) len_trim(str1)
       delims=' ' ; call parse(str1,delims,args,nargs)
       call value(args(1),mm0,ios) ; if(ios /= 0) mm0=0
       write(6,*) mm0
       read(17,'(a280)',err=711,end=799) str1
       delims=' ' ; call parse(str1,delims,args,nargs)
       call value(args(1),scale1,ios)
       write(6,*) scale1
       read(17,'(a280)',err=711,end=799) str1
       delims=' ' ; call parse(str1,delims,args,nargs)
       call value(args(1),cmat1(1,1),ios)
       call value(args(2),cmat1(1,2),ios)
       call value(args(3),cmat1(1,3),ios)
       read(17,'(a280)',err=711,end=799) str1
       delims=' ' ; call parse(str1,delims,args,nargs)
       call value(args(1),cmat1(2,1),ios)
       call value(args(2),cmat1(2,2),ios)
       call value(args(3),cmat1(2,3),ios)
       read(17,'(a280)',err=711,end=799) str1
       delims=' ' ; call parse(str1,delims,args,nargs)
       call value(args(1),cmat1(3,1),ios)
       call value(args(2),cmat1(3,2),ios)
       call value(args(3),cmat1(3,3),ios)
       if(scale1 > 0.d0)then
       cmat1=cmat1*scale1
                        endif
       if(scale1 < 0.d0)then
       cmat0=cmat1
       vtest=(cmat0(1,2)*cmat0(2,3)-cmat0(1,3)*cmat0(2,2))*cmat0(3,1) &
            +(cmat0(1,3)*cmat0(2,1)-cmat0(1,1)*cmat0(2,3))*cmat0(3,2) &
            +(cmat0(1,1)*cmat0(2,2)-cmat0(1,2)*cmat0(2,1))*cmat0(3,3)
       vtest=(abs(scale1)/abs(vtest))**(1.d0/3.d0)
       cmat1=cmat0*vtest
                        endif
       write(6,*) cmat1(1,1),cmat1(1,2),cmat1(1,3)
       write(6,*) cmat1(2,1),cmat1(2,2),cmat1(2,3)
       write(6,*) cmat1(3,1),cmat1(3,2),cmat1(3,3)
       read(17,'(a280)',err=711,end=799) str1
       delims=' ' ; call parse(str1,delims,args,nargs)
       nspecies1=nargs
       if(.not. allocated(sigmamatrix1)) allocate(sigmamatrix1(nspecies1,nspecies1))
       if(.not. allocated(symbl1)) allocate(symbl1(nspecies1))
       if(.not. allocated(nelements1)) allocate(nelements1(nspecies1))
!
       do j=1,nargs
       symbl1(j)=adjustl(trim(args(j)))
!      write(6,*) symbl1(j)
       enddo
       write(6,'(20(a2,1x))') (adjustl(symbl1(j)),j=1,nspecies1)
       read(17,'(a280)',err=711,end=799) str1
       delims=' ' ; call parse(str1,delims,args,nargs)
       natot=0
       do j=1,nargs
       call value(args(j),nelements1(j),ios)
       natot=natot+nelements1(j)
!      write(6,*) nelements1(j)
       enddo
       write(6,'(20(i4))') (nelements1(j),j=1,nspecies1)
       if(.not. allocated(dir)) allocate(dir(natot,3))
       if(.not. allocated(fxyz)) allocate(fxyz(natot,3))
       if(.not. allocated(lmove)) allocate(lmove(natot,3))
       if(.not. allocated(dirref)) allocate(dirref(natot,3))
       read(17,'(a280)',err=711,end=799) str1
       delims=' ' ; call parse(str1,delims,args,nargs)
       ldirect=.false.
       write(6,*) trim(args(1))
       ch1=args(1)(1:1)
       if(ch1 == 'K') args(1)='Cartesian'
       if(ch1 == 'k') args(1)='Cartesian'
       if(ch1 == 'C') args(1)='Cartesian'
       if(ch1 == 'c') args(1)='Cartesian'
       if(ch1 == 'D') args(1)='Direct'
       if(ch1 == 'd') args(1)='Direct'
       if(args(1) == 'Direct') ldirect=.true.
       if(args(1) == 'Cartesian') ldirect=.false.
!      write(6,*) ldirect
       do j=1,natot
!      read(17,'(a280)',err=711,end=799) str1
       read(17,'(a280)',err=711) str1
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
       lmove=.true. ; dirref=dir
!      write(6,*) 'lfault',lfault
       if(lfault) call write_lexit()
!
       open(7,file='../csa.in',form='formatted')
       read(7,*) 
       read(7,*)
       read(7,*)
       read(7,*) cellvol0,tmp,voltol
       read(7,*) refa1(1),refa1(2),refa1(3)
       read(7,*) refa2(1),refa2(2),refa2(3)
       read(7,*) refa3(1),refa3(2),refa3(3)
       zheight=refa3(3)
       write(6,*) 'zheight == refa3(3) from ../csa.in ',zheight
       write(6,*) 'cellvol0 (A^3), voltol from ../csa.in ',cellvol0,voltol
       do i=1,nspecies1
       read(7,*) (sigmamatrix1(i,j),j=1,nspecies1)
       enddo
       read(7,*) lvcs
       close(7)
       write(6,*) 'sigmamatrix1 (A) from ../csa.in'
       do i=1,nspecies1
       write(6,*) (sigmamatrix1(i,j),j=1,nspecies1)
       enddo
       if(.not. lvcs)then
       cmat1(1,:)=refa1(:) ; cmat1(2,:)=refa2(:) ; cmat1(3,:)=refa3(:)
                     endif
!
       z0frac=0.3d0 ; z1frac=0.7d0
       inquire(file='../POSCAR_ref',exist=lexist)
       if(.not. lexist)then
       write(6,*) z0frac,z1frac,' z0frac,z1frac, default'
       write(6,*) 'can be modified by ../POSCAR_ref'
                       endif
       if(lexist)then
       write(6,*) '../POSCAR_ref is present'
       open(8,file='../POSCAR_ref',form='formatted')
       read(8,'(a280)') str1
       delims=' ' ; call parse(str1,delims,args,nargs)
       if(nargs >= 2)then
       call value(args(1),z0frac,ios) ; if(ios /= 0) z0frac=0.3d0
       call value(args(2),z1frac,ios) ; if(ios /= 0) z1frac=0.7d0
       write(6,'(2f14.6,2x,a32)') z0frac,z1frac,'z0frac,z1frac from ../POSCAR_ref'
                     endif
       read(8,*) scale1
       read(8,*) refa1(1),refa1(2),refa1(3)
       read(8,*) refa2(1),refa2(2),refa2(3)
       read(8,*) refa3(1),refa3(2),refa3(3)
       if(scale1 > 0.d0)then
       refa1=refa1*scale1 ; refa2=refa2*scale1 ; refa3=refa3*scale1
                        endif
       if(scale1 < 0.d0)then
       cmat0(1,:)=refa1(:) ; cmat0(2,:)=refa2(:) ; cmat0(3,:)=refa3(:)
       vtest=(cmat0(1,2)*cmat0(2,3)-cmat0(1,3)*cmat0(2,2))*cmat0(3,1) &
            +(cmat0(1,3)*cmat0(2,1)-cmat0(1,1)*cmat0(2,3))*cmat0(3,2) &
            +(cmat0(1,1)*cmat0(2,2)-cmat0(1,2)*cmat0(2,1))*cmat0(3,3)
       vtest=(abs(scale1)/abs(vtest))**(1.d0/3.d0)
       cmat0=cmat0*vtest
       refa1(:)=cmat0(1,:) ; refa2(:)=cmat0(2,:) ; refa3(:)=cmat0(3,:)
                        endif
       if(.not. lvcs)then
       cmat1(1,:)=refa1(:) ; cmat1(2,:)=refa2(:) ; cmat1(3,:)=refa3(:)
                     endif
       read(8,'(a280)') str1
       delims=' ' ; call parse(str1,delims,args,nargs)
       if(nspecies1 /= nargs)then
       write(6,*) 'problem and/or inconsistency between ../csa.in and ../POSCAR_ref'
                             endif
       read(8,'(a280)') str1
       delims=' ' ; call parse(str1,delims,args,nargs)
       do i=1,nargs
       call value(args(i),j,ios)
       if(j /= nelements1(i))then
       write(6,*) 'problem and/or inconsistency between ../csa.in and ../POSCAR_ref'
                             endif
       enddo
       lselective=.false.
       read(8,'(a280)') str1
       delims=' ' ; call parse(str1,delims,args,nargs)
       ch1=args(1)(1:1)
       if(ch1 == 'S') args(1)='Selective'
       if(ch1 == 's') args(1)='Selective'
       if(args(1) == 'Selective' .or. args(1) == 'selective' .or. args(1) == 'S' .or. args(1) == 's')then
       lselective=.true.
       read(8,'(a280)') str1
       delims=' ' ; call parse(str1,delims,args,nargs)
                                                                                                     endif
       ldirect=.false.
       ch1=args(1)(1:1)
       if(ch1 == 'K') args(1)='Cartesian'
       if(ch1 == 'k') args(1)='Cartesian'
       if(ch1 == 'C') args(1)='Cartesian'
       if(ch1 == 'c') args(1)='Cartesian'
       if(ch1 == 'D') args(1)='Direct'
       if(ch1 == 'd') args(1)='Direct'
       if(args(1) == 'Direct') ldirect=.true.
       if(args(1) == 'Cartesian') ldirect=.false.
       lmove=.true.
       if(ldirect)then
       do i=1,natot
       read(8,'(a280)') str1
       delims=' ' ; call parse(str1,delims,args,nargs)
       if(.not. lselective)then
       call value(args(1),dirref(i,1),ios)
       call value(args(2),dirref(i,2),ios)
       call value(args(3),dirref(i,3),ios)
                           else
       call value(args(1),dirref(i,1),ios)
       call value(args(2),dirref(i,2),ios)
       call value(args(3),dirref(i,3),ios)
       if(args(4) == 'F' .or. args(4) == 'f') lmove(i,1)=.false.
       if(args(5) == 'F' .or. args(5) == 'f') lmove(i,2)=.false.
       if(args(6) == 'F' .or. args(6) == 'f') lmove(i,3)=.false.
                           endif
       enddo
                  else
       do i=1,natot
       read(8,'(a280)') str1
       delims=' ' ; call parse(str1,delims,args,nargs)
       if(.not. lselective)then
       call value(args(1),dirref(i,1),ios)
       call value(args(2),dirref(i,2),ios)
       call value(args(3),dirref(i,3),ios)
                           else
       call value(args(1),dirref(i,1),ios)
       call value(args(2),dirref(i,2),ios)
       call value(args(3),dirref(i,3),ios)
       if(args(4) == 'F' .or. args(4) == 'f') lmove(i,1)=.false.
       if(args(5) == 'F' .or. args(5) == 'f') lmove(i,2)=.false.
       if(args(6) == 'F' .or. args(6) == 'f') lmove(i,3)=.false.
                           endif
       enddo
       if(.not. lvcs)then
       cmat1(1,:)=refa1(:) ; cmat1(2,:)=refa2(:) ; cmat1(3,:)=refa3(:)
                     endif
       call tolatx1(cmat1,dirref,natot)
                  endif
       close(8)
       write(6,*) minval(abs(dirref(:,3))),'the lowest z value in ../POSCAR_ref'
       write(6,*) maxval(abs(dirref(:,3))),'the highest z value in ../POSCAR_ref'
                 endif
       end
!234567890
!      Written by In-Ho Lee, KRISS, December 30, 2015.
       subroutine mkcell_1()
       USE celldesign, ONLY : cmat1,dir,z0frac,z1frac,zheight,nspecies1,nelements1
       USE celldesign, ONLY : lvcs,cellvol0,voltol,refa1,refa2,refa3,dirref,lmove,lfirst
       implicit none
       real*8 cmat2(3,3),vec(3),tmp,tmq,tmr,tmz,pi,delx(3,3),rsclv(3,3)
       real*8 a1(3),a2(3),a3(3),areap,areat,vtest,xxr,xxi,yyr,yyi
       integer j,natot
       logical llattice,lomega
       real ranmar
       
       natot=0
       do j=1,nspecies1
       natot=natot+nelements1(j)
       enddo
!      you can use dirref, lmove, refa1, refa2, refa3, basically come from ../POSCAR_ref
!      z0frac=0.3d0 ; z1frac=0.7d0 ; zheight=15.0d0
       write(6,'(1x,a29,1x,f18.8,1x,a3)') '2D setting; zheight==refa3(3)',zheight,'(A)'
       write(6,'(1x,a10,f18.8,1x,a20)') '2D setting',(z0frac+(1.0d0-z1frac))*zheight,'vacuum thickness (A)'
!      reserve vacuum region by necessitation, do the relaxation by Monte Carlo
       do j=1,natot
       if(dir(j,3) > z1frac) dir(j,3)=ranmar()*(z1frac-z0frac)+z0frac
       if(dir(j,3) < z0frac) dir(j,3)=ranmar()*(z1frac-z0frac)+z0frac
       enddo
       do j=1,natot
       if((.not. lmove(j,1)) .and. (.not. lmove(j,2)) .and. (.not. lmove(j,3))) dir(j,:)=dirref(j,:)
       enddo
       if(.not. lvcs)then
       cmat2(1,:)=refa1(:) ; cmat2(2,:)=refa2(:) ; cmat2(3,:)=refa3(:)
                     else
       vec(1)=cmat1(1,1) ; vec(2)=cmat1(1,2) ; vec(3)=cmat1(1,3) ; xxr=sqrt(dot_product(vec,vec)) 
       vec(1)=cmat1(1,1) ; vec(2)=cmat1(1,2) ; vec(3)=0.0d0      ; yyr=sqrt(dot_product(vec,vec)) 
       cmat2(1,1)=cmat1(1,1)*xxr/yyr
       cmat2(1,2)=cmat1(1,2)*xxr/yyr
       vec(1)=cmat1(2,1) ; vec(2)=cmat1(2,2) ; vec(3)=cmat1(2,3) ; xxi=sqrt(dot_product(vec,vec)) 
       vec(1)=cmat1(2,1) ; vec(2)=cmat1(2,2) ; vec(3)=0.0d0      ; yyi=sqrt(dot_product(vec,vec)) 
       cmat2(2,1)=cmat1(2,1)*xxi/yyi
       cmat2(2,2)=cmat1(2,2)*xxi/yyi
       cmat2(3,3)=zheight
       cmat2(1,3)=0.0d0 ; cmat2(2,3)=0.0d0 ; cmat2(3,1)=0.0d0 ; cmat2(3,2)=0.0d0
!
  300  continue
       a1(:)=cmat2(1,:) ; a2(:)=cmat2(2,:) ; a3(:)=cmat2(3,:)
       tmp=(a1(2)*a2(3)-a1(3)*a2(2))*a3(1) &
          +(a1(3)*a2(1)-a1(1)*a2(3))*a3(2) &
          +(a1(1)*a2(2)-a1(2)*a2(1))*a3(3)
       tmp=abs(tmp) ; tmz=tmp
       lomega=.true.
       if(tmp > cellvol0*(1.d0+voltol*(1.0-0.5)*2.d0)) lomega=.false.
       if(tmp < cellvol0*(1.d0+voltol*(0.0-0.5)*2.d0)) lomega=.false.
       if(.not. lomega)then
       vtest=cellvol0*(1.d0+voltol*(ranmar()-0.5)*2.d0)
       areap=tmp/abs(cmat2(3,3)) ; areat=vtest/abs(cmat2(3,3)) ; tmq=areat/areap ; tmp=sqrt(tmq)
       cmat2(1,:)=cmat2(1,:)*tmp ; cmat2(2,:)=cmat2(2,:)*tmp
       write(6,*) 'volume is out of range, adjustment',tmz,vtest
                       else
       write(6,*) 'volume is within the range, no adjustment',tmz
                       endif
       rsclv=cmat2 
!
       call test_2dlattice(cmat2,llattice)
       if(.not. llattice)then
       delx=cmat2
       cmat2(3,1)=0.0d0 ; cmat2(3,2)=0.0d0 ; cmat2(3,3)=zheight ; pi=4.0d0*atan(1.0d0) 
!      vec(:)=cmat1(1,:) ; tmp=sqrt(dot_product(vec,vec)) 
       vec(:)=rsclv(1,:) ; tmp=sqrt(dot_product(vec,vec)) 
       tmz=0.15d0*ranmar() ; tmp=tmp*(1.0d0+tmz*(ranmar()-0.5)*2.0d0)
       if(tmp < 1.2d0) tmp=ranmar()+1.4d0
       tmq=pi*(ranmar()-0.5)*2.d0
       cmat2(1,1)=tmp*cos(tmq) ; cmat2(1,2)=tmp*sin(tmq) ; cmat2(1,3)=0.0d0
!      vec(:)=cmat1(2,:) ; tmp=sqrt(dot_product(vec,vec)) 
       vec(:)=rsclv(2,:) ; tmp=sqrt(dot_product(vec,vec)) 
       tmz=0.15d0*ranmar() ; tmp=tmp*(1.0d0+tmz*(ranmar()-0.5)*2.0d0)
       if(tmp < 1.2d0) tmp=ranmar()+1.4d0
       tmq=pi*(ranmar()-0.5)*2.d0
       cmat2(2,1)=tmp*cos(tmq) ; cmat2(2,2)=tmp*sin(tmq) ; cmat2(2,3)=0.0d0
       write(6,*) 'cell-edges are out of range, adjustment'
       write(6,'(3f22.12)') delx(1,1),delx(1,2),delx(1,3)
       write(6,'(3f22.12)') delx(2,1),delx(2,2),delx(2,3)
       write(6,'(3f22.12)') delx(3,1),delx(3,2),delx(3,3)
       goto 300
                         else
       write(6,*) 'cell-edges are within the range, no adjustment'
       write(6,'(3f22.12)') cmat2(1,1),cmat2(1,2),cmat2(1,3)
       write(6,'(3f22.12)') cmat2(2,1),cmat2(2,2),cmat2(2,3)
       write(6,'(3f22.12)') cmat2(3,1),cmat2(3,2),cmat2(3,3)
                         endif
       do j=1,natot
       vec(:)=dir(j,1)*cmat1(1,:)+dir(j,2)*cmat1(2,:)+dir(j,3)*cmat1(3,:)
       dir(j,:)=vec(:)
       enddo
       call tolatx1(cmat2,dir,natot)
                     endif
       write(6,'(2f20.8,2x,a32)') cellvol0,voltol*100.d0,'cellvol0 (A^3), voltol (percent)'
       a1(:)=cmat2(1,:) ; a2(:)=cmat2(2,:) ; a3(:)=cmat2(3,:)
       tmp=(a1(2)*a2(3)-a1(3)*a2(2))*a3(1) &
          +(a1(3)*a2(1)-a1(1)*a2(3))*a3(2) &
          +(a1(1)*a2(2)-a1(2)*a2(1))*a3(3)
       vtest=abs(tmp)
       tmp=cellvol0*(1.d0+voltol*(0.0-0.5)*2.d0)
       tmq=cellvol0*(1.d0+voltol*(1.0-0.5)*2.d0)
       tmr=(tmq-tmp)/(2.0d0*cellvol0)
       write(6,*) tmr*100.d0,' percent'
       write(6,'(f20.8,1x,a2,1x,f20.8,1x,a2,1x,f20.8)') tmp,'<=',vtest,'<=',tmq
       write(6,'(f20.8,1x,a2,1x,f20.8,1x,a2,1x,f20.8)') tmp/abs(cmat2(3,3)),'<=',vtest/abs(cmat2(3,3)),'<=',tmq/abs(cmat2(3,3))
       tmr=(tmq/abs(cmat2(3,3))-tmp/abs(cmat2(3,3)))/(2.0d0*cellvol0/abs(cmat2(3,3)))
       write(6,*) tmr*100.d0,' percent'
       tmr=(tmq/abs(cmat2(3,3))+tmp/abs(cmat2(3,3)))/2.0d0
       write(6,*) tmr,' (A^2)'
       cmat1=cmat2
       call corecore2d()
       if(lfirst)then
       write(6,*) 'not relaxation, but rectification/vacuum, lfirst',lfirst
                 else
       write(6,*) 'relaxation, rectification/vacuum, lfirst',lfirst
                 endif
       end
!234567890
!      Written by In-Ho Lee, KRISS, December 30, 2015.
       subroutine corecore2d()
       USE celldesign, ONLY : cmat1,dir,nspecies1,nelements1,sigmamatrix1,z0frac,z1frac,lmove,dirref
       implicit none
       real*8 vec(3),wec(3),dij,dv1,dv2,dv3,urep,urep0,urep_best,amp1
       integer i,j,iter,natot
       integer, allocatable :: is(:)
       real*8, allocatable :: dir_best(:,:),dir1(:,:),dir2(:,:)
       real ranmar

       natot=0
       do i=1,nspecies1
       natot=natot+nelements1(i)
       enddo
       allocate(is(natot)) ; allocate(dir_best(natot,3),dir1(natot,3),dir2(natot,3))
!
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
!
       do i=1,natot
       if((.not. lmove(i,1)) .and. (.not. lmove(i,2)) .and. (.not. lmove(i,3))) dir(i,:)=dirref(i,:)
       enddo
       dir_best=dir
       natot=0
       do i=1,nspecies1
       do j=1,nelements1(i)
       natot=natot+1 ; is(natot)=i
       enddo
       enddo
       urep0=0.0d0
       do i=1,natot-1
       do j=i+1,natot
       vec(1)=dir_best(i,1)-dir_best(j,1)
       vec(2)=dir_best(i,2)-dir_best(j,2)
       vec(3)=dir_best(i,3)-dir_best(j,3)
       wec(1)=vec(1)*cmat1(1,1)+vec(2)*cmat1(2,1)+vec(3)*cmat1(3,1)
       wec(2)=vec(1)*cmat1(1,2)+vec(2)*cmat1(2,2)+vec(3)*cmat1(3,2)
       wec(3)=vec(1)*cmat1(1,3)+vec(2)*cmat1(2,3)+vec(3)*cmat1(3,3)
       dij=sqrt(dot_product(wec,wec))
       if(dij < sigmamatrix1(is(i),is(j))) urep0=urep0+1.0d0
       enddo
       enddo
       do i=1,natot
       if(dir_best(i,3) > z1frac) urep0=urep0+1.0d0
       if(dir_best(i,3) < z0frac) urep0=urep0+1.0d0
       enddo
!      write(6,*) urep0,' urep0'
       urep_best=urep0
       wec(:)=cmat1(1,:) ; dv1=sqrt(dot_product(wec,wec))
       wec(:)=cmat1(2,:) ; dv2=sqrt(dot_product(wec,wec))
       wec(:)=cmat1(3,:) ; dv3=sqrt(dot_product(wec,wec))
       do iter=1,5000000*natot
       dir2=dir_best 
       if(ranmar() < 0.01)then
       amp1=ranmar()*0.1
       if(iter >  200) amp1=ranmar()*0.2
       if(iter >  300) amp1=ranmar()*0.3
       if(iter >  400) amp1=ranmar()*0.4
       if(iter >  500) amp1=ranmar()*0.5
       if(iter >  600) amp1=ranmar()*0.6
       if(iter >  700) amp1=ranmar()*0.7
       if(iter >  800) amp1=ranmar()*0.8
       if(iter >  900) amp1=ranmar()*0.9
       if(iter > 1000) amp1=ranmar()*2.0
       do i=1,natot
       dir2(i,1)=dir2(i,1)+amp1*(ranmar()-0.5)/dv1
       dir2(i,2)=dir2(i,2)+amp1*(ranmar()-0.5)/dv2
       dir2(i,3)=dir2(i,3)+amp1*(ranmar()-0.5)/dv3
       if(dir2(i,3) > z1frac) dir2(i,3)=ranmar()*(z1frac-z0frac)+z0frac
       if(dir2(i,3) < z0frac) dir2(i,3)=ranmar()*(z1frac-z0frac)+z0frac
       enddo
       do i=1,natot
       if((.not. lmove(i,1)) .and. (.not. lmove(i,2)) .and. (.not. lmove(i,3))) dir2(i,:)=dirref(i,:)
       enddo
                          endif
       do j=1,natot
       dir2(j,1)=dir2(j,1)-anint(dir2(j,1))
       dir2(j,2)=dir2(j,2)-anint(dir2(j,2))
       dir2(j,3)=dir2(j,3)-anint(dir2(j,3))
       enddo
       do j=1,natot
       if(dir2(j,1) < 0.d0) dir2(j,1)=dir2(j,1)+1.d0
       if(dir2(j,2) < 0.d0) dir2(j,2)=dir2(j,2)+1.d0
       if(dir2(j,3) < 0.d0) dir2(j,3)=dir2(j,3)+1.d0
       enddo
       dir1=dir2
       do i=1,natot-1
       do j=i+1,natot
       vec(1)=dir2(i,1)-dir2(j,1)
       vec(2)=dir2(i,2)-dir2(j,2)
       vec(3)=dir2(i,3)-dir2(j,3)
       wec(1)=vec(1)*cmat1(1,1)+vec(2)*cmat1(2,1)+vec(3)*cmat1(3,1)
       wec(2)=vec(1)*cmat1(1,2)+vec(2)*cmat1(2,2)+vec(3)*cmat1(3,2)
       wec(3)=vec(1)*cmat1(1,3)+vec(2)*cmat1(2,3)+vec(3)*cmat1(3,3)
       dij=sqrt(dot_product(wec,wec))
       if(dij < sigmamatrix1(is(i),is(j)))then
       amp1=ranmar()*2.5
       if(iter > 1000) amp1=ranmar()*3.5
       dir1(i,1)=dir2(i,1)+amp1*(ranmar()-0.5)/dv1
       dir1(i,2)=dir2(i,2)+amp1*(ranmar()-0.5)/dv2
       dir1(i,3)=dir2(i,3)+amp1*(ranmar()-0.5)/dv3
       dir1(j,1)=dir2(j,1)+amp1*(ranmar()-0.5)/dv1
       dir1(j,2)=dir2(j,2)+amp1*(ranmar()-0.5)/dv2
       dir1(j,3)=dir2(j,3)+amp1*(ranmar()-0.5)/dv3
       if(dir1(i,3) > z1frac) dir1(i,3)=ranmar()*(z1frac-z0frac)+z0frac
       if(dir1(i,3) < z0frac) dir1(i,3)=ranmar()*(z1frac-z0frac)+z0frac
       if(dir1(j,3) > z1frac) dir1(j,3)=ranmar()*(z1frac-z0frac)+z0frac
       if(dir1(j,3) < z0frac) dir1(j,3)=ranmar()*(z1frac-z0frac)+z0frac
                                          endif
       enddo
       enddo
       do i=1,natot
       if(dir1(i,3) > z1frac) dir1(i,3)=ranmar()*(z1frac-z0frac)+z0frac
       if(dir1(i,3) < z0frac) dir1(i,3)=ranmar()*(z1frac-z0frac)+z0frac
       enddo
       do i=1,natot
       if((.not. lmove(i,1)) .and. (.not. lmove(i,2)) .and. (.not. lmove(i,3))) dir1(i,:)=dirref(i,:)
       enddo
       do j=1,natot
       dir1(j,1)=dir1(j,1)-anint(dir1(j,1))
       dir1(j,2)=dir1(j,2)-anint(dir1(j,2))
       dir1(j,3)=dir1(j,3)-anint(dir1(j,3))
       enddo
       do j=1,natot
       if(dir1(j,1) < 0.d0) dir1(j,1)=dir1(j,1)+1.d0
       if(dir1(j,2) < 0.d0) dir1(j,2)=dir1(j,2)+1.d0
       if(dir1(j,3) < 0.d0) dir1(j,3)=dir1(j,3)+1.d0
       enddo
       urep=0.0d0
       do i=1,natot-1
       do j=i+1,natot
       vec(1)=dir1(i,1)-dir1(j,1)
       vec(2)=dir1(i,2)-dir1(j,2)
       vec(3)=dir1(i,3)-dir1(j,3)
       wec(1)=vec(1)*cmat1(1,1)+vec(2)*cmat1(2,1)+vec(3)*cmat1(3,1)
       wec(2)=vec(1)*cmat1(1,2)+vec(2)*cmat1(2,2)+vec(3)*cmat1(3,2)
       wec(3)=vec(1)*cmat1(1,3)+vec(2)*cmat1(2,3)+vec(3)*cmat1(3,3)
       dij=sqrt(dot_product(wec,wec))
       if(dij < sigmamatrix1(is(i),is(j))) urep=urep+1.0d0
       enddo
       enddo
       do i=1,natot
       if(dir1(i,3) > z1frac) urep=urep+1.0d0
       if(dir1(i,3) < z0frac) urep=urep+1.0d0
       enddo
       if(urep <= urep_best)then
       urep_best=urep
       dir_best=dir1
                            endif
       if(urep_best < 1.0d0) exit
       enddo
       dir=dir_best
       deallocate(dir_best,dir1,dir2) ; deallocate(is)
!      write(6,*) urep,' urep'
!      write(6,*) iter,' iter'
       write(6,'(f12.1,1x,a3,1x,f12.1,1x,i12)') urep0,'-->',urep,iter
       end
!234567890
       subroutine test_2dlattice(cmat,llattice)
       implicit none
       real*8 cmat(3,3)
       logical llattice
       real*8 ra,rb,rc,alpha,beta,gama,cosinea,cosineb,cosinec,tmp,pi

       ra=sqrt(cmat(1,1)**2+cmat(1,2)**2+cmat(1,3)**2)
       rb=sqrt(cmat(2,1)**2+cmat(2,2)**2+cmat(2,3)**2)
       rc=sqrt(cmat(3,1)**2+cmat(3,2)**2+cmat(3,3)**2)
       cosinea=(cmat(2,1)*cmat(3,1)+cmat(2,2)*cmat(3,2)+cmat(2,3)*cmat(3,3))/rb/rc
       cosineb=(cmat(1,1)*cmat(3,1)+cmat(1,2)*cmat(3,2)+cmat(1,3)*cmat(3,3))/rc/ra
       cosinec=(cmat(1,1)*cmat(2,1)+cmat(1,2)*cmat(2,2)+cmat(1,3)*cmat(2,3))/ra/rb
       pi=4.0d0*atan(1.0d0) ; tmp=180.0d0/pi
       alpha=tmp*acos(cosinea) ; beta=tmp*acos(cosineb) ; gama=tmp*acos(cosinec)
       llattice=.true.
       if(ra    < 1.2d0 .or. rb     < 1.2d0 ) llattice=.false.
       if(alpha < 20.0d0 .or. alpha > 160.0d0 ) llattice=.false.
       if(beta  < 20.0d0 .or. beta  > 160.0d0 ) llattice=.false.
       if(gama  < 20.0d0 .or. gama  > 160.0d0 ) llattice=.false.
       if(ra/rb > 6.0d0 .or. ra/rb  < 0.3d0 ) llattice=.false.
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
!      Written by In-Ho Lee, KRISS, December 30, 2015.
       subroutine read_outcar_design(lfault)
       USE celldesign, ONLY : mm0,cmat1,pp1,symbl1,dir,fxyz,nspecies1,nelements1
       USE celldesign, ONLY : lrlxion,lrlxcell,ldesign,lmove,lvcs,lfirst,llast
       USE strings, ONLY : parse,value
       implicit none
       logical lfault
       real*8 energy00_read,efermi,enthalpy,etot,prlxion,prlxcell,tmp,tmq
       real*8 a1(3),a2(3),a3(3),tma(3,3),tmb(3,3)
       real*8 extpress,dk
       logical lenth,lexist
       integer i,j,natot,iprint
       character*280 fname,kname
       integer ios,nargs
       character*280 str1
       character*280 args(40)
       character*20 delims

       ldesign=.true.
       lrlxion=.true.
       lrlxcell=.true. ; if(.not. lvcs) lrlxcell=.false.
       extpress=0.0d0
       prlxion=1.d-1
       prlxcell=1.d-3
       pp1=0.0d0 ; fxyz=0.0d0
       efermi=1.d8
       etot=1.d8
       enthalpy=1.d8
       lenth=.false.
       lfault=.false.
       fname='OUTCAR'
       inquire(file=trim(fname),exist=lexist)
       if(.not. lexist)then
       pp1=0.0d0 ; fxyz=0.0d0
       write(6,*) 'there is no calculation, there are no relaxations'
       goto 100
                       endif
       write(6,*) 'there is calculation, there are relaxations'
       open(18,file=trim(fname),form='formatted')
       do
       read(18,'(a280)',err=911,end=999) str1
!      write(6,*) len_trim(str1)
       delims=' ' ; call parse(str1,delims,args,nargs)
       if(nargs == 3)then
       if(args(1) == 'POSITION')    then
       if(args(2) == 'TOTAL-FORCE') then
       if(args(3) == '(eV/Angst)')  then
!      write(6,*) trim(args(3))
       read(18,'(a280)',err=911,end=999) str1
       natot=0
       do j=1,nspecies1
       natot=natot+nelements1(j)
       enddo
!      write(6,*) natot
       natot=0
       do 
       read(18,'(a280)',err=911,end=999) str1
       delims=' ' ; call parse(str1,delims,args,nargs)
       if(nargs == 6)then
       natot=natot+1
       call value(args(1),dir(natot,1),ios)  ; if(ios /= 0) dir(natot,1)=1.d-8
       call value(args(2),dir(natot,2),ios)  ; if(ios /= 0) dir(natot,2)=1.d-8
       call value(args(3),dir(natot,3),ios)  ; if(ios /= 0) dir(natot,3)=1.d-8
       call value(args(4),fxyz(natot,1),ios) ; if(ios /= 0) fxyz(natot,1)=1.d-8
       call value(args(5),fxyz(natot,2),ios) ; if(ios /= 0) fxyz(natot,2)=1.d-8
       call value(args(6),fxyz(natot,3),ios) ; if(ios /= 0) fxyz(natot,3)=1.d-8
                     else
                     exit
                     endif
     
       enddo
                                    endif
                                    endif
                                    endif
                     endif
!      write(6,*) natot
       if(nargs == 8)then
       if(args(1) == 'in')    then
       if(args(2) == 'kB')    then
       call value(args(3),pp1(1,1),ios) ; if(ios /= 0) pp1(1,1)=1.d-8
       call value(args(4),pp1(2,2),ios) ; if(ios /= 0) pp1(2,2)=1.d-8
       call value(args(5),pp1(3,3),ios) ; if(ios /= 0) pp1(3,3)=1.d-8
       call value(args(6),pp1(1,2),ios) ; if(ios /= 0) pp1(1,2)=1.d-8
       call value(args(7),pp1(2,3),ios) ; if(ios /= 0) pp1(2,3)=1.d-8
       call value(args(8),pp1(3,1),ios) ; if(ios /= 0) pp1(3,1)=1.d-8
       pp1(2,1)=pp1(1,2) ; pp1(3,2)=pp1(2,3) ; pp1(1,3)=pp1(3,1)
!      write(6,'(6f18.8)') pp1(1,1),pp1(2,2),pp1(3,3),pp1(1,2),pp1(2,3),pp1(3,1)
                              endif
                              endif
                     endif
       if(nargs == 6)then
       if(args(1) == 'direct')    then
       if(args(2) == 'lattice')   then
       if(args(3) == 'vectors')   then
       if(args(4) == 'reciprocal')then
       if(args(5) == 'lattice')   then
       if(args(6) == 'vectors')   then
       read(18,*,err=911,end=999) cmat1(1,1),cmat1(1,2),cmat1(1,3)
       read(18,*,err=911,end=999) cmat1(2,1),cmat1(2,2),cmat1(2,3)
       read(18,*,err=911,end=999) cmat1(3,1),cmat1(3,2),cmat1(3,3)
                                  endif
                                  endif
                                  endif
                                  endif
                                  endif
                                  endif
                     endif
       if(nargs == 7)then
       if(args(1) == 'energy'  )then
       if(args(2) == 'without' )then
       if(args(3) == 'entropy=')then
       call value(args(4),etot,ios) ; if(ios /= 0) etot=1.d8
!      print*, str1
!      print*, etot
                                endif
                                endif
                                endif
                     endif
       if(nargs == 9)then
       if(args(1) == 'enthalpy')then
       if(args(2) == 'is'      )then
       if(args(3) == 'TOTEN'   )then
       lenth=.true.
       call value(args(5),enthalpy,ios) ; if(ios /= 0) enthalpy=1.d8
!      print*, str1
!      print*, enthalpy
                                endif
                                endif
                                endif
                     endif
       if(nargs >= 7)then
       if(args(1) == 'E-fermi')then
       if(args(2) == ':'      )then
       call value(args(3),efermi,ios) ; if(ios /= 0) efermi=1.d8
!      print*, str1
!      print*, efermi
                               endif
                               endif
                     endif
       enddo
  911  continue
       lfault=.true.
  999  continue
       close(18)
       if(lfault)then
       write(6,*) 'there is a fault sign from stdout.log'
       call sleep(1)
                 endif
       if(lfault) call write_lexit()
       energy00_read=etot
       if(lenth) energy00_read=enthalpy
       if(efermi >= 1.d8 .or. energy00_read >= 1.d8) lfault=.true.
!
       iprint=1
       iprint=0
       if(iprint == 1)then
       write(6,*) 
       write(6,'(3f18.8,2x,3f18.8)') cmat1(1,1),cmat1(1,2),cmat1(1,3),pp1(1,1),pp1(1,2),pp1(1,3)
       write(6,'(3f18.8,2x,3f18.8)') cmat1(2,1),cmat1(2,2),cmat1(2,3),pp1(2,1),pp1(2,2),pp1(2,3)
       write(6,'(3f18.8,2x,3f18.8)') cmat1(3,1),cmat1(3,2),cmat1(3,3),pp1(3,1),pp1(3,2),pp1(3,3)
       write(6,'(20(a2,1x))') (adjustl(symbl1(j)),j=1,nspecies1)
       write(6,'(20(i4))') (nelements1(j),j=1,nspecies1)
                      endif
!      relaxation for ionic sites, note dir is in cartesian
       if(lrlxion) call ion_rlx(fxyz,dir,lmove,prlxion,natot,nspecies1,nelements1)
!      relaxation for cell-shape
       tma=cmat1
       if(lrlxcell) call cell_rlx(pp1,cmat1,extpress,prlxcell)
       tmb=cmat1-tma
       tmp=maxval(abs(fxyz)) ; tmq=maxval(abs(tmb))
!
       iprint=1
       iprint=0
       if(iprint == 1)then
       natot=0
       do i=1,nspecies1
       do j=1,nelements1(i)
       natot=natot+1
       write(6,'(3f18.8,2x,3f18.8)') dir(natot,1),dir(natot,2),dir(natot,3),fxyz(natot,1),fxyz(natot,2),fxyz(natot,3)
       enddo
       enddo
                      endif
       call tolatx1(cmat1,dir,natot)
!      note dir is in direct, the reference lattice vectors are given by matrix cmat1
!      cell design
  100  continue
!      a new trial structure that is relaxed and/or designed
!      here, a quasi-two dimensional system is considered
       if(.not. lfirst)then
       tmp=maxval(abs(fxyz)) ; tmq=maxval(abs(pp1))
       write(6,*) 'convergence test, because it is not the first step',lfirst
       write(6,'(2f14.4,2x,a7)') tmp,tmq,'tmp,tmq'
       if(tmp < 0.01d0 .and. tmq < 1.50d0)then
       write(6,'(2f18.8,1x,a29)') tmp,tmq,' ---converged calculation----'
       llast=.true.
       call write_lexit()
                                          endif
                       else
       write(6,*) 'no convergence test, because it is the first step',lfirst
                       endif
!      a new trial structure that is relaxed and/or designed
!      if(lfirst) call mkcell_1()
       call mkcell_1()
!
       fname='POSCAR_design'
       open(71,file=trim(fname),form='formatted')
       write(71,*) mm0
       write(71,'(a3)') '1.0'
       write(71,'(3f23.16)') cmat1(1,1),cmat1(1,2),cmat1(1,3)
       write(71,'(3f23.16)') cmat1(2,1),cmat1(2,2),cmat1(2,3)
       write(71,'(3f23.16)') cmat1(3,1),cmat1(3,2),cmat1(3,3)
       write(71,'(20(2x,a2,1x))') (symbl1(i),i=1,nspecies1)
       write(71,'(20(i4,1x))') (nelements1(i),i=1,nspecies1)
       write(71,'(a6)') "Direct"
       if(lfirst) write(6,*) 'designed cell'
       if(.not. lfirst) write(6,*) 'relaxation in the designed cell'
       write(6,*) mm0
       write(6,'(a3)') '1.0'
       write(6,'(3f23.16)') cmat1(1,1),cmat1(1,2),cmat1(1,3)
       write(6,'(3f23.16)') cmat1(2,1),cmat1(2,2),cmat1(2,3)
       write(6,'(3f23.16)') cmat1(3,1),cmat1(3,2),cmat1(3,3)
       write(6,'(20(2x,a2,1x))') (symbl1(i),i=1,nspecies1)
       write(6,'(20(i4,1x))') (nelements1(i),i=1,nspecies1)
       write(6,'(a6)') "Direct"
       natot=0
       do i=1,nspecies1
       do j=1,nelements1(i)
       natot=natot+1
       write(71,'(3f20.16)') dir(natot,1),dir(natot,2),dir(natot,3)
       write(6,'(3f20.16)') dir(natot,1),dir(natot,2),dir(natot,3)
       enddo
       enddo
       close(71)
       write(6,*) 
       call print2dlattice(cmat1)
       a1(:)=cmat1(1,:) ; a2(:)=cmat1(2,:) ; a3(:)=cmat1(3,:)
       dk=0.70d0 ; kname='KPOINTS_070' ; call write_kpoints2d(dk,a1,a2,a3,kname)
       dk=0.50d0 ; kname='KPOINTS_050' ; call write_kpoints2d(dk,a1,a2,a3,kname)
       dk=0.40d0 ; kname='KPOINTS_040' ; call write_kpoints2d(dk,a1,a2,a3,kname)
       dk=0.30d0 ; kname='KPOINTS_030' ; call write_kpoints2d(dk,a1,a2,a3,kname)
       dk=0.25d0 ; kname='KPOINTS_025' ; call write_kpoints2d(dk,a1,a2,a3,kname)
       dk=0.20d0 ; kname='KPOINTS_020' ; call write_kpoints2d(dk,a1,a2,a3,kname)
       dk=0.10d0 ; kname='KPOINTS_010' ; call write_kpoints2d(dk,a1,a2,a3,kname)
       end
!23456789
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine write_kpoints2d(dk,a1,a2,a3,kname)
       implicit none
       character*280 kname
       real*8 dk,a1(3),a2(3),a3(3)
       real*8 b1(3),b2(3),b3(3),omega,dum,pi,ga,gb,gc
       integer ka,kb,kc,iswitch

       iswitch=0
       if(dk < 0.d0) dk=0.06d0
       if(dk == 0.d0)then
       dk=0.015d0
       iswitch=1
                     endif
!
       call cross3(a2,a3,b1) ; call cross3(a3,a1,b2) ; call cross3(a1,a2,b3)
       omega=abs(dot_product(b1,a1)) ; pi=4.0d0*atan(1.0d0)
       b1=b1*(2.d0*pi/omega) ; b2=b2*(2.d0*pi/omega) ; b3=b3*(2.d0*pi/omega)
       ga=sqrt(dot_product(b1,b1)) ; gb=sqrt(dot_product(b2,b2)) ; gc=sqrt(dot_product(b3,b3))
       call kmesh2d(dk,ga,dum,ka) ; call kmesh2d(dk,gb,dum,kb) ; call kmesh2d(dk,gc,dum,kc)
       if(ka > 15) ka=15
       if(kb > 15) kb=15
       if(kc > 15) kc=15
       if(iswitch == 1)then
       ka=ka+1+ka/2
       kb=kb+1+kb/2
       kc=kc+1+kc/2
                       endif
       kc=1
       open(71,file=trim(kname),form='formatted')
       write(71,"(a1)") 'A'
       write(71,"(i1)")  0
       write(71,"(a1)") 'G'
       write(71,"(1x,i3,2x,i3,2x,i3)") ka,kb,kc
       write(71,"(1x,i2,2x,i2,2x,i2)") 0,0,0
!      write(71,"(1x,3f4.2)") 0.,0.,0.
       close(71)
       end
!23456789
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine kmesh2d(dk,ga,dum,ka)
       implicit none
       real*8 dk,ga,dum
       integer ka
       real*8 pi,gatmp
       integer i

       pi=4.0d0*atan(1.0d0)
!      gatmp=ga/(2.d0*pi)
       gatmp=ga
       ka=int(gatmp/dk) ; if(ka == 0) ka=1
       dum=gatmp/dble(ka)
       if(dum >= dk)then
       do i=1,15
       ka=ka+1 ; dum=gatmp/dble(ka)
       if(dum <= dk) exit
       enddo
                    endif
       end
!234567890
!      Written by In-Ho Lee, KRISS, December 30, 2015.
       subroutine cell_rlx(pp,amat,extpress,prlxcell)
       implicit none
       real*8 pp(3,3),amat(3,3),extpress,prlxcell
       real*8 a1(3),a2(3),a3(3),tms(3,3),tmt(3,3)
       real*8 h0(3,3),sg(3,3),vec(3),ei(3,3),ss(3,3),rr(3,3),rlxmax,tmp
       integer i,j

       write(6,*) 'cell-shape relaxtion, quasi 2D version'
!      write(6,*) pp,' pp'
       write(6,*)
!      write(6,*) amat,' before'
       tms=amat
       rlxmax=0.10d0
       a1(:)=amat(1,:) ; a2(:)=amat(2,:) ; a3(:)=amat(3,:)
       call h0aa(1,h0,a1,a2,a3)
       call cross3(a2,a3,vec) ; sg(:,1)=vec(:)
       call cross3(a3,a1,vec) ; sg(:,2)=vec(:)
       call cross3(a1,a2,vec) ; sg(:,3)=vec(:)
       ei=0.0d0 ; ei(1,1)=1.0d0 ; ei(2,2)=1.0d0 ; ei(3,3)=1.0d0 ; rr=ei*extpress
       ei=pp-rr ; ss=matmul(ei,sg)
!      write(6,*)
!      write(6,*) sg,' sg'
!      write(6,*)
!      write(6,*) ss,' ss'
       do j=1,3
       do i=1,3
       if(i == 3 .and. j == 1) cycle
       if(i == 3 .and. j == 2) cycle
       if(i == 3 .and. j == 3) cycle
       if(i == 1 .and. j == 3) cycle
       if(i == 2 .and. j == 3) cycle
       tmp=prlxcell*ss(j,i) ; if(tmp < -rlxmax) tmp=-rlxmax ; if(tmp > rlxmax) tmp=rlxmax
       h0(j,i)=h0(j,i)+tmp
!      write(6,*) j,i,tmp
       enddo
       enddo
       call h0aa(-1,h0,a1,a2,a3)
       amat(1,:)=a1(:) ; amat(2,:)=a2(:) ; amat(3,:)=a3(:)
       write(6,*)
!      write(6,*) amat,' after'
       tmt=amat
       tmp=maxval(abs(pp))
       write(6,*) tmp,' max component in stress tensor'
       write(6,*) maxval(abs(tmt-tms)),' max component in diff'
       end
!234567890
!      Written by In-Ho Lee, KRISS, December 30, 2015.
       subroutine ion_rlx(fxyz,xyz,lmove,prlxion,natot,nspecies1,nelements1)
       implicit none
       integer natot,nspecies1,nelements1(nspecies1)
       real*8 fxyz(natot,3),xyz(natot,3),prlxion
       logical lmove(natot,3)
       integer i,j,k,na
       real*8 tmp,rlxmax

       write(6,*) 'ionic relaxtion'
!      write(6,*) xyz, ' before'
       write(6,*)
!      note xyz is in cartesian
       rlxmax=0.10d0
       na=0
       do i=1,nspecies1
       do j=1,nelements1(i)
       na=na+1
!      write(6,*) na
       do k=1,3
       if(lmove(na,k))then
       tmp=prlxion*fxyz(na,k) ; if(tmp < -rlxmax) tmp=-rlxmax ; if(tmp > rlxmax) tmp=rlxmax
       xyz(na,k)=xyz(na,k)+tmp
!      write(6,*) k,tmp
                      endif
       enddo
       enddo
       enddo
!      write(6,*) xyz, ' after'
       tmp=maxval(abs(fxyz))
       write(6,*) tmp,' max component in fxyz'
       end
!234567890
!      Written by In-Ho Lee, KRISS, December 30, 2015.
      subroutine h0aa(jsign,h0,a1,a2,a3)
       implicit none
       integer jsign
       real*8 h0(3,3),a1(3),a2(3),a3(3)

       if(jsign ==  1)then
       h0(1,1)=a1(1) ; h0(2,1)=a1(2) ; h0(3,1)=a1(3)
       h0(1,2)=a2(1) ; h0(2,2)=a2(2) ; h0(3,2)=a2(3)
       h0(1,3)=a3(1) ; h0(2,3)=a3(2) ; h0(3,3)=a3(3)
                      endif
       if(jsign == -1)then
       a1(1)=h0(1,1) ; a1(2)=h0(2,1) ; a1(3)=h0(3,1)
       a2(1)=h0(1,2) ; a2(2)=h0(2,2) ; a2(3)=h0(3,2)
       a3(1)=h0(1,3) ; a3(2)=h0(2,3) ; a3(3)=h0(3,3)
                      endif
       end
!234567890
!      Written by In-Ho Lee, KRISS, December 30, 2015.
       subroutine cross3(a,b,c)
       implicit none
       real*8 a(3),b(3),c(3)

       c(1)=a(2)*b(3)-a(3)*b(2)
       c(2)=a(3)*b(1)-a(1)*b(3)
       c(3)=a(1)*b(2)-a(2)*b(1)
       end
!234567890
!      Written by In-Ho Lee, KRISS, December 30, 2015.
       subroutine diff_poscar(lfault,lconv)
       USE celldesign, ONLY : nspecies1,nelements1
       implicit none
       logical lfault,lconv
       real*8 tmp,tmq,vec1(3),vec2(3),vec3(3)
       integer j,k,natot
       logical lexist
       
       lfault=.false.
       lconv=.false.
       natot=0
       do j=1,nspecies1
       natot=natot+nelements1(j)
       enddo
       tmp=-1.d15
       tmq=-1.d15
       inquire(file='POSCAR_design',exist=lexist)
       if(.not. lexist)then
       return
                       endif
       inquire(file='POSCAR',exist=lexist)
       if(.not. lexist)then
       return
                       endif
       open(1,file='POSCAR',form='formatted') 
       open(2,file='POSCAR_design',form='formatted') 
       read(1,*,err=911,end=911)
       read(1,*,err=911,end=911)
       read(2,*,err=911,end=911)
       read(2,*,err=911,end=911)
       do k=1,3
       read(1,*,err=911,end=911) vec1(1),vec1(2),vec1(3)
       read(2,*,err=911,end=911) vec2(1),vec2(2),vec2(3)
       vec3=vec1-vec2 ; tmp=max(tmp,maxval(abs(vec3)))
       enddo
       read(1,*,err=911,end=911)
       read(1,*,err=911,end=911)
       read(1,*,err=911,end=911)
       read(2,*,err=911,end=911)
       read(2,*,err=911,end=911)
       read(2,*,err=911,end=911)
       do j=1,natot
       read(1,*,err=911) vec1(1),vec1(2),vec1(3)
       read(2,*,err=911) vec2(1),vec2(2),vec2(3)
       vec3=vec1-vec2 ; tmq=max(tmq,maxval(abs(vec3)))
       enddo
       goto 999
  911  continue
       lfault=.true.
  999  continue
       close(1)
       close(2)
       write(6,*) tmp,tmq
       if(tmp < 0.0001d0 .and. tmq < 0.0001d0) lconv=.true.
       end
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine read_stdout_design(lfault_stdout)
       implicit none
       character*280 stdname
       logical lfault_stdout
       character*6 ctest6,c6
       character*2 c2
       character*4 c4,cc4
       character*5 c5
       character*7 c7
       character*8 c8
       character*9 c9
       character*10 c10
       character*11 c11
       character*14 c14
       logical lfault
       logical lexist

       stdname='stdout.log'
       lfault=.false.
       inquire(file=trim(stdname),exist=lexist)
       if(.not. lexist)then
       lfault_stdout=lfault
       return
                       endif
       open(18,file=trim(stdname),form='formatted')
       do
       read(18,*,err=911,end=999) ctest6

       if(ctest6 == 'BRMIX:')then
       backspace(18)
       read(18,*,err=911,end=999) ctest6, c4,c7,c8
       if(trim(c4) == 'very' .and. trim(c7) == 'serious' .and. trim(c8) == 'problems')then
       print*, c4,' ',c7,' ',c8
       goto 911
                                                                                      endif
                             endif
       if(ctest6 == 'intern')then
       backspace(18)
       read(18,*,err=911,end=999) c8, c5, c14
       if(trim(c8) == 'internal' .and. trim(c5) == 'ERROR' .and. trim(c14) ==  'RSPHER:running' )then
       print*, c8,' ',c5,' ',c14
       goto 911
                                                                                                 endif
                             endif
       if(ctest6 == 'APPLIC')then
       backspace(18)
       read(18,*,err=911,end=999) c11, c10
       if(trim(c11) == 'APPLICATION' .and. trim(c10) =='TERMINATED')then
       print*, c11,' ',c10
       goto 911
                                                                    endif
                             endif
       if(ctest6 == 'Ctrl-C')then
       backspace(18)
       read(18,*,err=911,end=999) ctest6, c9
       if(trim(c9) == 'caught...' )then
       print*, ctest6,' ',c9
       goto 911
                                   endif
                             endif
!
       if(ctest6 == 'ERROR:')then
       backspace(18)
       read(18,*,err=911,end=999) ctest6, c5
       if(trim(c5) == 'while' )then
       print*, ctest6,' ',c5
       goto 911
                               endif
       backspace(18)
       read(18,*,err=911,end=999) ctest6, c6
       if(trim(c6) == 'charge')then
       print*, ctest6,' ',c6
       goto 911
                               endif
                             endif
       if(ctest6 == 'NKPT>N')then
       backspace(18)
       read(18,*,err=911,end=999) c10
       if(trim(c10) == 'NKPT>NKDIM')then
       print*, c10
       goto 911
                                    endif
                             endif
       if(trim(ctest6) == 'exit')then
       backspace(18)
       read(18,*,err=911,end=999) c4,ctest6,c2,cc4
       if(trim(ctest6) == 'status' .and. trim(c2) == 'of' .and. trim(cc4) == 'rank')then
       print*, c4,' ',ctest6,' ',c2,' ',cc4
       goto 911
                                                                                    endif
                                 endif
!
       enddo
  911  continue
       lfault=.true.
  999  continue
       close(18)
!      write(6,*) 'in stdout.log', lfault
!
!      if(lfault)then
!      write(6,*) 'there is a fault sign from stdout.log'
!      call sleep(3)
!                endif
       lfault_stdout=lfault
       call read_stdout_log10(lfault_stdout)
       end
!234567890
!      http://error.wiki/VASP
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine read_stdout_log10(lfault_stdout)
       USE strings, ONLY : parse,value
       implicit none
       character*280 stdname
       logical lfault_stdout
       integer nargs
       character*280 str1
       character*280 args(40)
       character*20 delims
       logical lfault
       logical lexist

       stdname='stdout.log'
       lfault=.false.
       inquire(file=trim(stdname),exist=lexist)
       if(.not. lexist)then
       lfault_stdout=lfault
       return
                       endif
       open(18,file=trim(stdname),form='formatted')
       do
       read(18,'(a280)',err=911,end=999) str1
!      write(6,*) len_trim(str1)
       delims=' '
       call parse(str1,delims,args,nargs)
       if(nargs > 0)then
       if(args(1) == 'No' .and. args(2) == 'initial' .and. args(3) == 'positions') goto 911
       if(args(1) == 'NKPT>NKDIM') goto 911
       if(args(1) == 'BRMIX:' .and. args(2) == 'very' .and. args(3) == 'serious' .and. args(4) == 'problems') goto 911
       if(args(1) == 'ERROR:' .and. args(2) == 'while') goto 911
       if(args(1) == 'ERROR:' .and. args(2) == 'charge') goto 911
       if(args(1) == 'ERROR:' .and. args(2) == 'there') goto 911
       if(args(1) == 'ERROR:' .and. args(2) == 'missing') goto 911
       if(args(1) == 'Error'  .and. args(2) == 'EDDDAV:' .and. args(5) == 'ZHEGV') goto 911
       if(args(1) == 'Error'  .and. args(2) == 'EDDRMM:' .and. args(5) == 'ZHEGV') goto 911
       if(args(5) == 'segmentation' .and. args(6) == 'fault') goto 911
       if(args(1) == 'ERROR'   .and. args(3) == 'supplied') goto 911
       if(args(1) == 'ERROR'   .and. args(2) == 'code' .and. args(3) == 'was') goto 911
       if(args(1) == 'ERROR'   .and. args(2) == 'in' .and. args(3) == 'subspace') goto 911
       if(args(1) == 'LAPACK:'  .and. args(2) == 'Routine' .and. args(3) == 'ZPOTRF') goto 911
       if(args(1) == 'internal' .and. args(2) == 'error') goto 911
       if(args(2) == 'internal' .and. args(3) == 'error') goto 911
       if(args(4) == 'internal' .and. args(5) == 'error') goto 911
       if(args(1) == 'internal' .and. args(2) == 'ERROR' .and. args(3) == 'RSPHER:running') goto 911
       if(args(1) == 'Hard' .and. args(2) == 'potentials') goto 911
       if(args(1) == 'Suspicious' .and. args(2) == 'behaviour') goto 911
       if(args(1) == 'Large' .and. args(2) == 'positive' .and. args(3) == 'energies') goto 911
       if(args(1) == 'APPLICATION' .and. args(2) == 'TERMINATED') goto 911
       if(args(1) == 'Ctrl-C' .and. args(2) == 'caught...') goto 911
       if(args(1) == 'exit' .and. args(2) == 'status' .and. args(3) == 'of' .and. args(4) == 'rank') goto 911
       if(args(1) == 'integer' .and. args(2) == 'divide' .and. args(3) == 'by') goto 911
       if(args(1) == 'Calculation' .and. args(2) == 'hangs' .and. args(3) == 'at') goto 911
       if(args(1) == 'Fatal' .and. args(2) == 'error' .and. args(3) == 'in') goto 911
       if(args(4) == 'accuracy' .and. args(5) == 'cannot' .and. args(6) == 'be') goto 911
       if(args(3) == 'accuracy' .and. args(4) == 'cannot' .and. args(5) == 'be') goto 911
                    endif
       enddo
  911  continue
       lfault=.true.
  999  continue
       close(18)
       if(lfault)then
       write(6,*) 'there is a fault sign from stdout.log'
       call sleep(3)
                 endif
       lfault_stdout=lfault
       if(lfault) call write_lexit()
       end
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine init_seed()
       implicit none
       integer n,ival(8),jv(3),i
       integer, allocatable :: kseed(:)

       call date_and_time(values=ival)
       jv(1) = ival(8) + 2048*ival(7)
       jv(2) = ival(6) + 64*ival(5) ! value(4) isn't really random
       jv(3) = ival(3) + 32*ival(2) + 32*8*ival(1)
       call random_seed(size=n)
       allocate(kseed(n))
       call random_seed()          ! Give the seed an implementation-dependent kick
       call random_seed(get=kseed)
       do i=1, n
       kseed(i) = kseed(i) + jv(mod(i-1, 3) + 1)
       enddo
       call random_seed(put=kseed)
       deallocate(kseed)
       end
!23456789
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine f90sleep(dtinsec)
       implicit none
       real*8 :: dtinsec         ! desired sleep interval [s]
       integer,dimension(8) :: t ! arguments for date_and_time
       integer :: ms1,ms2  ! start and end times [ms]
       real*8 :: dt              ! desired sleep interval [ms]

       dt=dtinsec*1.d3
       call date_and_time(values=t)
       ms1=(t(5)*3600+t(6)*60+t(7))*1000+t(8)

       do
         call date_and_time(values=t)
         ms2=(t(5)*3600+t(6)*60+t(7))*1000+t(8)
         if(ms2-ms1>=dt)exit
       enddo
       end
!23456789
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
       pi=4.0d0*atan(1.0d0) ; tmp=180.0d0/pi
       alpha=tmp*acos(cosinea) ; beta=tmp*acos(cosineb) ; gama=tmp*acos(cosinec)
       write(6,'(6f16.5)') ra,rb,rc,alpha,beta,gama
       a1(:)=cmat(1,:) ; a2(:)=cmat(2,:) ; a3(:)=cmat(3,:)
       tmp=(a1(2)*a2(3)-a1(3)*a2(2))*a3(1) &
          +(a1(3)*a2(1)-a1(1)*a2(3))*a3(2) &
          +(a1(1)*a2(2)-a1(2)*a2(1))*a3(3)
       tmp=abs(tmp) ; areap=tmp/abs(cmat(3,3)) 
       write(6,'(f20.6,1x,a5,2x,f20.6,1x,a5)') tmp, '(A^3)', areap, '(A^2)'
       end
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine csa_vasp_banner1()
       implicit none
       integer i
       character*93 banner(8)

       banner(1)='       `..     `.. ..        `.             `..         `..      `.         `.. ..  `....... '
       banner(2)='  `..   `..`..    `..     `. ..            `..       `..      `. ..     `..    `..`..    `.. '
       banner(3)=' `..        `..          `.  `..            `..     `..      `.  `..     `..      `..    `.. '
       banner(4)=' `..          `..       `..   `..   `.....   `..   `..      `..   `..      `..    `.......   '
       banner(5)=' `..             `..   `...... `..            `.. `..      `...... `..        `.. `..        '
       banner(6)='  `..   `..`..    `.. `..       `..            `....      `..       `.. `..    `..`..        '
       banner(7)='    `....    `.. ..  `..         `..            `..      `..         `..  `.. ..  `..        '
       banner(8)='                                                                                     20130911'
       do i=1,8
       write(6,'(15x,a93)') banner(i)
       enddo
       banner(1)='    #     #     #     #     ######   #######  #     #   #####  '
       banner(2)='   # #    ##   ##    # #    #     #  #        #     #  #     # '
       banner(3)='  #   #   # # # #   #   #   #     #  #        #     #  #       '
       banner(4)=' #     #  #  #  #  #     #  #     #  #####    #     #   #####  '
       banner(5)=' #######  #     #  #######  #     #  #        #     #        # '
       banner(6)=' #     #  #     #  #     #  #     #  #        #     #  #     # '
       banner(7)=' #     #  #     #  #     #  ######   #######   #####    #####  '
       banner(8)='             Ab initio MAterials DEsign Using cSa              '
       do i=1,8
       write(6,'(30x,a63)') banner(i)
       enddo
       end
!
