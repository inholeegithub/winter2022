!      Written by In-Ho Lee, KRISS, September 11, 2013.
!
!      Conformational Space Annealing (CSA) with First-Principles Electronic Structure Calculations
!      Atomic positions and six lattice parameters (a,b,c,alpha,beta,gamma) are dynamical variables for the CSA.
!      References: Phys. Rev. Lett. 91, 080201 (2003);  Phys. Rev. B 90, 115209 (2014).
!
!       ######   ######     ###            ##     ##    ###     ######  ########
!      ##    ## ##    ##   ## ##           ##     ##   ## ##   ##    ## ##     ##
!      ##       ##        ##   ##          ##     ##  ##   ##  ##       ##     ##
!      ##        ######  ##     ## ####### ##     ## ##     ##  ######  ########
!      ##             ## #########          ##   ##  #########       ## ##     
!      ##    ## ##    ## ##     ##           ## ##   ##     ## ##    ## ##    
!       ######   ######  ##     ##            ###    ##     ##  ######  ##   
!           
!               #     #     #     #     ######   #######  #     #   ##### 
!              # #    ##   ##    # #    #     #  #        #     #  #     #
!             #   #   # # # #   #   #   #     #  #        #     #  #  
!            #     #  #  #  #  #     #  #     #  #####    #     #   ##### 
!            #######  #     #  #######  #     #  #        #     #        #
!            #     #  #     #  #     #  #     #  #        #     #  #     #
!            #     #  #     #  #     #  ######   #######   #####    #####
!                        Ab initio MAterials DEsign Using cSa           
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       module csa_application
       implicit none
       private
       save
       integer natom,nft,nlopt
       integer iobj
       integer nspecies
       real*8 refa1(3),refa2(3),refa3(3),refvol,voltol,extpress,au2ev,au2mbar,au2ang
       real*8 rc1,rc2,shift
       logical lvcs,lpbc
       integer, allocatable :: nelements(:),itype(:),ncoord(:)
       character*2, allocatable :: symbl(:)
       real*8, allocatable :: sigmamatrix(:,:)
       real*8, allocatable :: wrk2(:),wrk4(:)
       integer, allocatable :: iwrk2(:),iwrk4(:)
       public :: iobj
       public :: natom,nft,nlopt,rc1,rc2,shift,ncoord
       public :: iwrk2,wrk2,iwrk4,wrk4,refvol,voltol
       public :: lvcs,lpbc,nspecies,symbl,nelements,sigmamatrix,itype
       public :: au2mbar,au2ev,au2ang,refa1,refa2,refa3,extpress
       end module csa_application
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       module csa
       implicit none
       private
       save
       integer ndeg,npop,npop1,nmate,npert,nevol,idiff,jdiff,nfrac,ndirectory
       integer iseed1,iseed2,ndeg_r,npop_r,npop1_r,mmdummy
       real*8 dcut,davg,amp,drate,energy_best
       real*8, allocatable :: posi(:,:),posi1(:,:),energy_sorted(:),energy_sorted1(:),posi_best(:)
       real*8, allocatable :: qosi(:,:),qosi1(:,:),wrk1(:),qosi0(:),qosi00(:),prev(:)
       real*8, allocatable :: posi_r(:,:),posi1_r(:,:),energy_sorted_r(:),energy_sorted1_r(:),posi_best_r(:)
       integer, allocatable :: iwrk1(:)
       real*8 energy0,energy_best_r,davg_r
       logical literative
       logical lquit
       character*280 cwd
       public :: csa_initial,csa_first_bank,csa_evolution,csa_final,lquit

       contains
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine csa_initial()
       USE csa_application, ONLY : iobj
       USE csa_application, ONLY : nft,nlopt,natom,rc1,rc2,shift,refvol,voltol
       USE csa_application, ONLY : iwrk2,wrk2,iwrk4,wrk4,ncoord,lpbc,lvcs
       USE csa_application, ONLY : nspecies,nelements,symbl,itype,sigmamatrix
       USE csa_application, ONLY : refa1,refa2,refa3,extpress,au2ev,au2mbar,au2ang
       implicit none
       integer idiff0,npop0,npop10,nmate0,npert0,nfrac0,nevol0,iseed10,iseed20
       real*8 amp0,drate0,factor
       real*8 cmatrix(3,3),s6(6),cellvol0,extpress0,vtest
       real*8 covlaentrr
       integer i,j,na,ish
       logical lexist,lnewjob

       read(5,*) nspecies
       allocate(symbl(nspecies),nelements(nspecies))
       allocate(sigmamatrix(nspecies,nspecies))
       read(5,*) (symbl(i),i=1,nspecies)
       read(5,*) (nelements(i),i=1,nspecies)
       read(5,*) cellvol0,extpress0,voltol
       read(5,*) cmatrix(1,1),cmatrix(1,2),cmatrix(1,3)
       read(5,*) cmatrix(2,1),cmatrix(2,2),cmatrix(2,3)
       read(5,*) cmatrix(3,1),cmatrix(3,2),cmatrix(3,3)
       do i=1,nspecies
       symbl(i)=adjustl(symbl(i))
       read(5,*) (sigmamatrix(i,j),j=1,nspecies)
       enddo
       read(5,*) lvcs,lpbc,iobj
       read(5,*) ndirectory
!      read(5,*) cwd
!      cwd=adjustl(cwd) 
       read(5,'(a280)') cwd
       cwd=adjustl(cwd)
       do i=1,280
       if(cwd(i:i) == ' ')then
       j=i
                          exit
                          endif
       enddo
       do i=j,280
       cwd(i:i)=' '
       enddo
       cwd=trim(cwd)
       i=len_trim(cwd) ; if(cwd(i:i) /= '/') cwd=trim(cwd)//'/' ; cwd=trim(cwd)
       read(5,*) idiff0,nevol0
       read(5,*) npop0,npop10
       npop10=npop0
       read(5,*) nmate0,npert0,nfrac0
       read(5,*) amp0,drate0
       read(5,*) iseed10,iseed20,lnewjob
       jdiff=1 
       if(idiff0 < 0)then
       idiff0=iabs(idiff0)
       jdiff=-1
                     endif
       if(ndirectory <= 0) ndirectory=npop0
       call init_seed()
       if(iseed10 <= 0 .or. iseed20 <= 0)then
       call random_number(vtest)
       iseed10=vtest*31328.d0+1.d0
       call random_number(vtest)
       iseed20=vtest*30081.d0+1.d0
                                         endif
!
       inquire(file='POTCAR',exist=lexist)
       if(.not. lexist)then
       write(6,*) 'POTCAR is not present.'
                       stop
                       endif
       inquire(file='INCAR_rlx',exist=lexist)
       if(.not. lexist)then
       write(6,*) 'INCAR_rlx is not present.'
                       stop
                       endif
       inquire(file='INCAR_rlxall',exist=lexist)
       if(.not. lexist)then
       write(6,*) 'INCAR_rlxall is not present.'
                       stop
                       endif
       inquire(file='INCAR_bs',exist=lexist)
       if(.not. lexist)then
       write(6,*) 'INCAR_bs is not present.'
                       stop
                       endif
       inquire(file='CSA_SOLDIER.pbs',exist=lexist)
       if(.not. lexist)then
       write(6,*) 'CSA_SOLDIER.pbs is not present.'
                       stop
                       endif
!
       factor=1.0d0
       na=0
       do i=1,nspecies
       do j=1,nspecies
       if(sigmamatrix(i,j) <= 0.0d0)then
       na=1
       factor=min(factor,abs(sigmamatrix(i,j)))
       sigmamatrix(i,j)=abs(sigmamatrix(i,j))
                                    endif
       enddo
       enddo
       do i=1,nspecies
       do j=1,nspecies
       if(abs(sigmamatrix(i,j)) < 1.d-1)then
       sigmamatrix(i,j)=covlaentrr(symbl(i))+covlaentrr(symbl(j))
       sigmamatrix(i,j)=sigmamatrix(i,j)*0.4d0
                                        endif
       enddo
       enddo
       if(na > 0)then
       if(factor > 1.0d0) factor=1.0d0
       if(factor < 0.1d0) factor=0.1d0
       do i=1,nspecies
       do j=1,nspecies
       sigmamatrix(i,j)=covlaentrr(symbl(i))+covlaentrr(symbl(j))
       enddo
       enddo
       sigmamatrix=abs(factor)*sigmamatrix
                 endif
       do i=1,nspecies
       do j=1,nspecies
       if(sigmamatrix(i,j) < 0.11d0) sigmamatrix(i,j)=0.11d0 
       enddo
       enddo
       sigmamatrix=(sigmamatrix+transpose(sigmamatrix))/2.0d0
!
       i=0
       i=1
       if(i==1)then
       write(6,*) nspecies
       write(6,'(20(2x,a2,1x))') (symbl(i),i=1,nspecies)
       write(6,'(20(i4,1x))') (nelements(i),i=1,nspecies)
       write(6,'(3f18.8)') cellvol0,extpress0,voltol
       write(6,'(3f22.12)') cmatrix(1,1),cmatrix(1,2),cmatrix(1,3)
       write(6,'(3f22.12)') cmatrix(2,1),cmatrix(2,2),cmatrix(2,3)
       write(6,'(3f22.12)') cmatrix(3,1),cmatrix(3,2),cmatrix(3,3)
       do i=1,nspecies
       write(6,'(20f10.4)') (sigmamatrix(i,j),j=1,nspecies)
       enddo
       write(6,*) 'sigmamatrix in Angstrom'
       do i=1,nspecies
       do j=1,nspecies
       write(6,'(2x,a2,2x,a2,f12.4)') symbl(i),symbl(j),sigmamatrix(i,j)
       enddo
       enddo
       if(iobj == 0) write(6,*) 'enthalpy minimization (See INCAR)'
       if(iobj == 1) write(6,*) 'direct band gap optimization (See objective funtion)'
       if(iobj == 2) write(6,*) 'electronic DOS at Fermi level maximization (See objective funtion)'
       if(iobj == 3) write(6,*) 'electronic DOS slope at Fermi level maximization (See objective funtion)'
       if(iobj == 4) write(6,*) 'electronic DOS derived effective mass maximization (See objective funtion)'
       write(6,*) lvcs,lpbc,iobj
       write(6,*) ndirectory
       write(6,*) trim(cwd)
       write(6,*) idiff0,nevol0
       write(6,*) npop0,npop10
       write(6,*) nmate0,npert0,nfrac0
       write(6,*) amp0,drate0
       write(6,*) iseed10,iseed20,lnewjob
               endif
!
       if(lnewjob)then
       inquire(file='fort.1',exist=lexist)
       if(lexist)then
       open(1,file='fort.1',form='formatted')
       close(1,status='delete')
       write(6,*) 'fort.1 is deleted.'
                 endif
                  else
       inquire(file='fort.1',exist=lexist)
       if(.not. lexist)then
       write(6,*) 'fort.1 is not present. Thus, this is a new job.'
       lnewjob=.true.
                       else
       write(6,*) 'fort.1 is present. This is an iterative job.'
                       endif
                  endif
!
       na=0
       do j=1,nspecies
       na=na+nelements(j)
       enddo
       allocate(itype(na))
       na=0
       do i=1,nspecies
       do j=1,nelements(i)
       na=na+1
       itype(na)=i
       enddo
       enddo
       natom=na
!
       if(.not. lpbc) lvcs=.false.
       if(lvcs) lpbc=.true.
       call latmat(s6,cmatrix,0)
       call latmatvol(s6,cmatrix,cellvol0)
       call latmat(s6,cmatrix,0)
       vtest=(cmatrix(1,2)*cmatrix(2,3)-cmatrix(1,3)*cmatrix(2,2))*cmatrix(3,1) &
            +(cmatrix(1,3)*cmatrix(2,1)-cmatrix(1,1)*cmatrix(2,3))*cmatrix(3,2) &
            +(cmatrix(1,1)*cmatrix(2,2)-cmatrix(1,2)*cmatrix(2,1))*cmatrix(3,3)
       vtest=abs(vtest)
       print*, vtest,' vtest,cellvol0',cellvol0
       print*, cmatrix(1,1),cmatrix(1,2),cmatrix(1,3)
       print*, cmatrix(2,1),cmatrix(2,2),cmatrix(2,3)
       print*, cmatrix(3,1),cmatrix(3,2),cmatrix(3,3)
       if(abs(vtest-cellvol0) > 1.d-8)then
       write(6,*) 'something went wrong'
                                      stop
                                      endif
!
       au2mbar=2.9421912d2 ; au2ev=13.6058d0*2.d0 ; au2ang=0.529177d0
       extpress=extpress0/au2mbar    ! extpress0 is in units of Mbar
!
       refa1(:)=cmatrix(1,:) ; refa2(:)=cmatrix(2,:) ; refa3(:)=cmatrix(3,:) ; refvol=cellvol0
       if(voltol < 0.d0)then
       voltol=abs(voltol)
                        endif
       print*, voltol,'voltol (%)',voltol*100.d0
!
       iseed1=iseed10 ; iseed2=iseed20
       iseed1=mod(iseed1,31328+1) ; iseed2=mod(iseed2,30081+1)
       call rmarin(iseed1,iseed2)
       rc1=0.d0 ; rc2=0.d0
       do i=1,nspecies
       do j=1,nspecies
       rc1=rc1+sigmamatrix(i,j)
       rc2=rc2+1.d0
       enddo
       enddo
       rc1=rc1/rc2
!!     rc1=sigmamatrix(1,1)*1.45d0
       rc1=rc1*1.45d0
       rc2=rc1*(1.70d0/1.35d0)
       shift=0.0d0
!
       ndeg=3*natom+6
       npop=npop0
       npop1=npop10
       if(npop  <= 0) npop=20
       if(npop1 <= 0) npop1=npop
       nmate=nmate0*npop
       npert=npert0*npop
       if(nmate <= 0) nmate=80*npop
       if(npert <= 0) npert=20*npop
       nevol=nevol0
       if(nevol <= 0) nevol=30
       amp=amp0
       nfrac=nfrac0
       if(nfrac <= 0) nfrac=4
       if(natom <= 4)then
       nfrac=2
       write(6,*) 'natom,nfrac',natom,nfrac
                     endif
       drate=drate0
!
       idiff=idiff0
!      write(6,*) nspecies,'nspecies,na',na
       write(6,'(20(2x,a2,1x))') (symbl(j),j=1,nspecies)
       write(6,'(20(i4))') (nelements(j),j=1,nspecies)
       write(6,*) lvcs,'lvcs,lpbc',lpbc
       write(6,'(2i5,2x,a24,2x,2i5)') npop,npop1,'npop,npop1,iseed1,iseed2',iseed1,iseed2
       write(6,'(2x,a5,2x,i5)') 'idiff',idiff
!      write(6,'(2f10.5,2x,3i7,1x,a33,1x,i7)') amp,drate,nmate,npert,nfrac,'amp,drate,nmate,npert,nfrac,nevol',nevol
       write(6,'(f6.3,f15.9,1x,3i7,1x,a33,1x,i7)') amp,drate,nmate,npert,nfrac,'amp,drate,nmate,npert,nfrac,nevol',nevol
       write(6,'(3f13.6,2x,a13)') rc1,rc2,shift,'rc1,rc2,shift'
!
       allocate(posi(ndeg,npop)) ; allocate(posi1(ndeg,npop1))
       allocate(energy_sorted(npop)) ; allocate(energy_sorted1(npop1))
       allocate(qosi(ndeg,npop),qosi1(ndeg,npop1))
       allocate(posi_best(ndeg))
       allocate(qosi0(ndeg),qosi00(ndeg))
       allocate(iwrk1(max(npop,npop1,3))) ; allocate(wrk1(max(npop,npop1,3)))
       allocate(prev(npop))
       allocate(ncoord(natom))
       allocate(wrk2(natom)) ; allocate(iwrk2(natom))
       allocate(wrk4(natom*natom)) ; allocate(iwrk4(natom*natom))
       energy_sorted=1.d20 ; energy_sorted1=1.d20
       ish=ndeg-6
       do i=1,6
       posi(ish+i,:)=s6(i)
       posi1(ish+i,:)=s6(i)
       enddo
       literative=.false.
       inquire(file='fort.1',exist=lexist)
       if(lexist)then
       call csa_bank_dump(1)
       if(ndeg_r == ndeg) literative=.true.
       write(6,*) 'iterative, npop1_r,npop1',npop1_r,npop1
       write(6,*) 'iterative, npop_r,npop',npop_r,npop
       write(6,*) 'iterative, energy_best_r',energy_best_r
       j=min(npop1,npop1_r)
       posi1(:,1:j)=posi1_r(:,1:j)
       energy_sorted1(1:j)=energy_sorted1_r(1:j)
       j=min(npop,npop_r)
       posi(:,1:j)=posi_r(:,1:j)
       energy_sorted(1:j)=energy_sorted_r(1:j)
       posi_best=posi_best_r
       energy_best=energy_best_r
                 endif
       nft=0 ; nlopt=0 ; energy_best=1.d23
       open(7,file='fort.7',form='formatted')
       open(8,file='fort.8',form='formatted')
       call flush(6)
       mmdummy=0
       inquire(file='fort.2',exist=lexist)
       if(lexist)then
       open(2,file='fort.2',form='formatted',status='old',position='append')
                 else
       open(2,file='fort.2',form='formatted')
                 endif
       end subroutine csa_initial
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine csa_final()
       USE csa_application, ONLY : iwrk2,wrk2,iwrk4,wrk4
       USE csa_application, ONLY : nft,nlopt,nelements,symbl,natom,sigmamatrix,itype,ncoord
       implicit none

       call onedprint10(energy_sorted,npop)
       write(6,*) nft,nlopt,' nft,nlopt',energy_best,natom
       write(7,*) nft,nlopt,' nft,nlopt',energy_best,natom
       write(8,'(i5,2x,f19.9,2x,2i9)') natom,energy_best,nft,nlopt
       call flush(8)
       call csa_bank_dump(0)
!
       if(allocated(posi_r)) deallocate(posi_r)
       if(allocated(posi1_r)) deallocate(posi1_r)
       if(allocated(energy_sorted_r)) deallocate(energy_sorted_r)
       if(allocated(energy_sorted1_r)) deallocate(energy_sorted1_r)
       if(allocated(posi_best_r)) deallocate(posi_best_r)
       deallocate(posi,posi1) ; deallocate(qosi,qosi1)
       deallocate(posi_best)
       deallocate(qosi0,qosi00)
       deallocate(energy_sorted) ; deallocate(energy_sorted1)
       deallocate(prev)
       deallocate(ncoord)
       deallocate(sigmamatrix)
       deallocate(symbl,nelements,itype)
       deallocate(iwrk1,iwrk2,iwrk4)
       deallocate(wrk1,wrk2,wrk4)
       close(7)
       close(8)
       close(2)
       end subroutine csa_final
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine csa_first_bank()
       USE csa_application, ONLY : natom,lvcs,lpbc,refa1,refa2,refa3,voltol
       USE csa_application, ONLY : refvol,nspecies,nelements,symbl,sigmamatrix
       implicit none
       integer i,j,ish,ispgrp
       real*8 tmq,tmr,tmp,amatrix(3,3),s6(6),cellvol0,vtest,s60(6)
       logical lflagls
       real*8, allocatable :: sigmamatrix0(:,:)
       real ranmar
 
       ish=ndeg-6
       if(lpbc)then
       amatrix(1,:)=refa1(:) ; amatrix(2,:)=refa2(:) ; amatrix(3,:)=refa3(:)
       cellvol0=(amatrix(1,2)*amatrix(2,3)-amatrix(1,3)*amatrix(2,2))*amatrix(3,1) &
               +(amatrix(1,3)*amatrix(2,1)-amatrix(1,1)*amatrix(2,3))*amatrix(3,2) &
               +(amatrix(1,1)*amatrix(2,2)-amatrix(1,2)*amatrix(2,1))*amatrix(3,3)
       cellvol0=abs(cellvol0)
       call latmat(s6,amatrix,0)
       s60=s6
       do i=1,6
       posi(ish+i,:)=s6(i)
       posi1(ish+i,:)=s6(i)
       enddo
               endif

!      dynamic variables in a random initialization mode
       if(lpbc)      then
       do j=1,npop
       if(lvcs)then
       vtest=cellvol0*(1.d0+voltol*(ranmar()-0.5)*2.d0)
       call gen_lattice_matrix(amatrix,s6,vtest)
       do i=1,6
       qosi0(ish+i)=s6(i)
       enddo
               else
       do i=1,6
       qosi0(ish+i)=s60(i)
       enddo
               endif
       do i=1,ndeg-6
       qosi0(i)=ranmar()
       enddo
       call danglingbond_care()
!
       if(ranmar() < 1.10)then
       if(lpbc)then
       qosi00(:)=qosi0(:)
       ispgrp=0
       allocate(sigmamatrix0(nspecies,nspecies))
       sigmamatrix0=sigmamatrix*0.5d0
       call gen_latt_site(ispgrp,ndeg,nspecies,nelements,symbl,sigmamatrix0,voltol,refvol,qosi0,lpbc,lvcs,lflagls)
       deallocate(sigmamatrix0)
       if(.not. lflagls) qosi0(:)=qosi00(:)
               endif
                          endif
!
       posi(:,j)=qosi0(:)
       enddo
       do j=1,npop1
       if(lvcs)then
       vtest=cellvol0*(1.d0+voltol*(ranmar()-0.5)*2.d0)
       call gen_lattice_matrix(amatrix,s6,vtest)
       do i=1,6
       qosi0(ish+i)=s6(i)
       enddo
               else
       do i=1,6
       qosi0(ish+i)=s60(i)
       enddo
               endif
       do i=1,ndeg-6
       qosi0(i)=ranmar()
       enddo
       call danglingbond_care()
!
       if(ranmar() < 1.10)then
       if(lpbc)then
       qosi00(:)=qosi0(:)
       ispgrp=0
       allocate(sigmamatrix0(nspecies,nspecies))
       sigmamatrix0=sigmamatrix*0.5d0
       call gen_latt_site(ispgrp,ndeg,nspecies,nelements,symbl,sigmamatrix0,voltol,refvol,qosi0,lpbc,lvcs,lflagls)
       deallocate(sigmamatrix0)
       if(.not. lflagls) qosi0(:)=qosi00(:)
               endif
                          endif
!
       posi1(:,j)=qosi0(:)
       enddo
                     endif
       if(.not. lpbc)then
       tmp=amp
       do j=1,npop
       do i=1,natom*3
       qosi0(i)=(ranmar()-0.5)*tmp
       enddo
       call danglingbond_care()
       posi(:,j)=qosi0(:)
       enddo
       do j=1,npop1
       do i=1,natom*3
       qosi0(i)=(ranmar()-0.5)*tmp
       enddo
       call danglingbond_care()
       posi1(:,j)=qosi0(:)
       enddo
                     endif
!
       if(.not. lvcs)then
       amatrix(1,:)=refa1(:) ; amatrix(2,:)=refa2(:) ; amatrix(3,:)=refa3(:)
       call latmat(s6,amatrix,0)
       do i=1,6
       posi(ish+i,:)=s6(i)
       posi1(ish+i,:)=s6(i)
       enddo
       if(lpbc) write(6,*) 'reference lattice vectors are used, lvcs',lvcs
                     endif
       if(.not. lpbc)then
       amatrix(1,:)=refa1(:) ; amatrix(2,:)=refa2(:) ; amatrix(3,:)=refa3(:)
       call latmat(s6,amatrix,0)
       do i=1,6
       posi(ish+i,:)=s6(i)
       posi1(ish+i,:)=s6(i)
       enddo
       write(6,*) 'reference lattice vectors are never used, lpbc',lpbc
       write(6,*) 'it is a nominal one'
                     endif
!
       if(npop <= npop1)then
       if(.not. literative)then
       call master_slave(npop1,ndirectory,-1)
       j=min(npop,npop1)
       posi(:,1:j)=posi1(:,1:j) ; energy_sorted(1:j)=energy_sorted1(1:j)
                           else
       call master_slave(npop1-npop1_r,ndirectory,-1)
       do j=1,npop1_r
       if(npop1-j+1 <1) exit
       posi1(:,npop1-j+1)=posi1_r(:,j)
       energy_sorted1(npop1-j+1)=energy_sorted1_r(j)
       enddo
       do j=1,npop_r
       if(npop-j+1 <1) exit
       posi(:,npop-j+1)=posi_r(:,j)
       energy_sorted(npop-j+1)=energy_sorted_r(j)
       enddo
       do j=1,npop1-npop1_r
       if(j <= npop)then
       posi(:,j)=posi1(:,j)
       energy_sorted(j)=energy_sorted1(j)
                    endif
       enddo
                           endif
       call csa_bank_sort(1)
       call csa_bank_sort(0)
                        endif
!
       if(npop > npop1)then
       if(.not. literative)then
       call master_slave(npop,ndirectory,1)
       j=min(npop,npop1)
       posi1(:,1:j)=posi(:,1:j) ; energy_sorted1(1:j)=energy_sorted(1:j)
                           else
       call master_slave(npop-npop_r,ndirectory,1)
       do j=1,npop1_r
       if(npop1-j+1 <1) exit
       posi1(:,npop1-j+1)=posi1_r(:,j)
       energy_sorted1(npop1-j+1)=energy_sorted1_r(j)
       enddo
       do j=1,npop_r
       if(npop-j+1 <1) exit
       posi(:,npop-j+1)=posi_r(:,j)
       energy_sorted(npop-j+1)=energy_sorted_r(j)
       enddo
       do j=1,npop-npop_r
       if(j <= npop1)then
       posi1(:,j)=posi(:,j)
       energy_sorted1(j)=energy_sorted(j)
                     endif
       enddo
                           endif
       call csa_bank_sort(0)
       call csa_bank_sort(1)
                       endif
!
       call csa_difference(1,0,davg,tmq,tmr)
       write(6,'(f18.8,2x,a4)') davg,'davg'
       dcut=davg/2.0d0
       if(literative)then
       if(davg_r > 0.0d0) dcut=davg_r/2.0d0
       write(6,*) 'iterative,davg_r,davg',davg_r,davg
                     endif
!
       call csa_update_best(0)
       write(6,'(1x,e22.12,2x,a37)') energy_best,'objective functions in the first bank'
       call onedprint10(energy_sorted,npop)
       end subroutine csa_first_bank
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine csa_difference(imode,jctl,davg,dsig,dista)
       USE csa_application, ONLY : nspecies,nelements,symbl,sigmamatrix,rc1,rc2,lpbc
       USE prdf, ONLY : prdf_init,prdf_final,get_prdf,prdf_cmp
       USE bldist, ONLY : bldist_init,bldist_final,get_blsrtd,bldist_cmp
       USE qlabmod, ONLY : qlab_init,qlab_final,get_qlab,qlab_cmp
       USE elemdist, ONLY : elemdist_init,elemdist_final,get_hist,elemdist_cmp
       implicit none
       integer imode,jctl
       real*8 davg,dsig,dista
       real*8 rmax0,tmp,tmq
       integer j,i222,ktmp

       rmax0=6.d0
       if(imode == 1 .or. imode == 2)then
       do j=1,npop
       if(idiff == 1)then
       if(j == 1) call elemdist_init(rmax0,nspecies,nelements,rc1,rc2,npop,lpbc)
       call get_hist(j,posi(1,j))
                     endif
       if(idiff == 2)then
       i222=1
       if(j == 1) call bldist_init(rmax0,jdiff,i222,nelements,nspecies,npop,lpbc)
       call get_blsrtd(j,posi(1,j))
                     endif
       if(idiff == 3)then
       i222=2
       if(j == 1) call bldist_init(rmax0,jdiff,i222,nelements,nspecies,npop,lpbc)
       call get_blsrtd(j,posi(1,j))
                     endif
       if(idiff == 4)then
       if(j == 1) call qlab_init(rmax0,nspecies,nelements,sigmamatrix,npop,lpbc)
       call get_qlab(j,posi(1,j))
                     endif
       if(idiff == 5)then
       if(j == 1) call prdf_init(rmax0,nspecies,nelements,symbl,npop,lpbc)
       call get_prdf(j,posi(1,j))
                     endif
       if(idiff == 6)then
       ktmp=-npop
       if(j == 1) call prdf_init(rmax0,nspecies,nelements,symbl,ktmp,lpbc)
       call get_prdf(j,posi(1,j))
                     endif
       enddo
                                     endif
       if(imode == 1)then
       if(idiff == 1) call elemdist_final(1,davg,dsig)
       if(idiff == 2) call bldist_final(1,davg,dsig)
       if(idiff == 3) call bldist_final(1,davg,dsig)
       if(idiff == 4) call qlab_final(1,davg,dsig)
       if(idiff == 5) call prdf_final(1,davg,dsig)
       if(idiff == 6) call prdf_final(1,davg,dsig)
                     endif
       if(imode == 2)then
       if(jctl == 0)then
       if(idiff == 1) call get_hist(0,qosi0)
       if(idiff == 2) call get_blsrtd(0,qosi0)
       if(idiff == 3) call get_blsrtd(0,qosi0)
       if(idiff == 4) call get_qlab(0,qosi0)
       if(idiff == 5) call get_prdf(0,qosi0)
       if(idiff == 6) call get_prdf(0,qosi0)
                    endif
                     endif
       if(imode == 3)then
       if(jctl /= 0)then
       if(idiff == 1) call elemdist_cmp(0,jctl,dista)
       if(idiff == 2) call bldist_cmp(0,jctl,dista)
       if(idiff == 3) call bldist_cmp(0,jctl,dista)
       if(idiff == 4) call qlab_cmp(0,jctl,dista)
       if(idiff == 5) call prdf_cmp(0,jctl,dista)
       if(idiff == 6) call prdf_cmp(0,jctl,dista)
                    endif
                     endif
       if(imode == -3)then
       if(idiff == 1) call elemdist_final(0,tmp,tmq)
       if(idiff == 2) call bldist_final(0,tmp,tmq)
       if(idiff == 3) call bldist_final(0,tmp,tmq)
       if(idiff == 4) call qlab_final(0,tmp,tmq)
       if(idiff == 5) call prdf_final(0,tmp,tmq)
       if(idiff == 6) call prdf_final(0,tmp,tmq)
                      endif
       end subroutine csa_difference
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine csa_evolution()
       USE csa_application, ONLY : nft,nlopt,natom
       implicit none
       integer ievol,iok
       real*8 test1

       iok=0
       prev=energy_sorted
       ievol=0
       write(6,*) dcut,' dcut, ievol',ievol
       do ievol=1,nevol

       call csa_perturbation_mate()
       wrk1(1:npop)=energy_sorted(1:npop)-prev(1:npop)
       test1=maxval(abs(wrk1(1:npop)))
       write(6,'(i8,f13.6,1x,f18.9,2x,2i10,2x, e11.3)') ievol,dcut,energy_best,nft,nlopt,test1
       write(7,'(i8,f13.6,1x,f18.9,2x,2i10,2x, e11.3)') ievol,dcut,energy_best,nft,nlopt,test1
       prev=energy_sorted
!      if(iok == 1)then
!      write(6,'(i5,1x,a66,1x,i5)') natom,'OK +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',ievol
!      write(8,'(i5,1x,a66,1x,i5)') natom,'OK +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',ievol
!      exit
!                  endif

       dcut=dcut*abs(drate)
       if(dcut < davg*1.d-3) dcut=davg*1.d-3
       write(6,*) dcut,' dcut, ievol',ievol
       call flush(6)
       call flush(7)
       enddo
       if(iok == 0)then
       write(6,'(i5,1x,a66,1x,i5)') natom,'************ ************ ***************** ******** *** *********',ievol
       write(8,'(i5,1x,a66,1x,i5)') natom,'************ ************ ***************** ******** *** *********',ievol
                   endif
       end subroutine csa_evolution
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine csa_update_conformations()
       implicit none

       call csa_keep_diversity()
       call csa_update_best(0)
       end subroutine csa_update_conformations
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine csa_keep_diversity()
       implicit none
       integer j,j0,kase
       real*8 dista,dista0,sstt,tmq,tmr,tms
       logical lsimilar

       call csa_difference(2,0,tmq,tmr,tms)
       j0=1 ; dista0=2.0d22
       lsimilar=.false.
       do j=1,npop
       call csa_difference(3,j,tmq,tmr,dista)
       if(dista < dcut)then
       lsimilar=.true.
       if(dista0 > dista)then
       dista0=dista ; j0=j
                         endif
                       endif
       enddo
       call csa_difference(-3,0,tmq,tmr,tms)
!
       kase=0
       if(.not. lsimilar)then
       kase=-1
       if(energy_sorted(npop) > energy0)then
       sstt=energy_sorted(npop) 
       posi(:,npop)=qosi0(:)
       energy_sorted(npop)=energy0
       write(6,'(1x,a20,1x,3e20.8)') 'new type, introduced',energy0,sstt,energy0-sstt
       kase=1
                                        endif
                         else
       kase=-2
       if(energy_sorted(j0) > energy0)then
       sstt=energy_sorted(j0)
       posi(:,j0)=qosi0(:)
       energy_sorted(j0)=energy0
       write(6,'(1x,a18,3x,3e20.8)') 'old type, replaced',energy0,sstt,energy0-sstt
       kase=2
                                      endif
                         endif
       call csa_bank_sort(0)
       mmdummy=mmdummy+1
       call contcarseries(mmdummy,energy0,kase)
       end subroutine csa_keep_diversity
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine csa_update_best(kprint)
       implicit none
       integer kprint
       integer j
!      integer iok

       if(energy_best > energy_sorted(1))then
       energy_best=energy_sorted(1)
       posi_best(:)=posi(:,1)  
                                         endif
!      if(iok == 1)then
!      lquit=.true.
!                  endif
       if(kprint > 0)then
       write(6,'(5f18.8,2x,a4)') (energy_sorted(j),j=1,5), 'best'
       call onedprint10(energy_sorted,npop)
                     endif
       end subroutine csa_update_best
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine csa_perturbation_mate()
       implicit none
       integer nwork

       lquit=.false.
       nwork=npert+nmate
       call master_slave(nwork,ndirectory,2)
       end subroutine csa_perturbation_mate
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine csa_rnd_lattice_basis(icm,nwork,iseq)
       USE csa_application, ONLY : iobj
       USE csa_application, ONLY : natom,shift,wrk2,iwrk2,lvcs,lpbc,nspecies,nelements,itype,refvol,voltol
       USE csa_application, ONLY : symbl,sigmamatrix
       implicit none
       integer icm,nwork,iseq(nwork)
       real*8 vec(3),uec(3),rx(3,3),ry(3,3),rz(3,3),ss(3,3),tt(3,3),ddv1,ddv2,ddv3,shift00,dv1,dv2,dv3
       real*8 a1(3),a2(3),a3(3),cmatrix(3,3),amatrix(3,3),r6(6),s6(6),t6(6),tmp,tmq,xmix
       real*8 sig0,xl0,xl,cellvol0,xxr,xxi,vtest,wec(3,2)
       integer ish,i,i1,i2,nb,jfrac,j,ii,jj,kk,ll,islice,jslice,kslice,lslice,jtr,ity
       integer icase,nptsize,jcase,jmutation,jparents(2)
       integer ispgrp
       logical llattice,lcate,lflagls
       real*8, allocatable :: sigmamatrix0(:,:)
       real ranmar
 
       icase=1
       icase=2
       tmp=amp*2.d0
       nptsize=5
       jcase=3
       jcase=2
       jcase=1   
!
       ish=ndeg-6 ; lcate=.false.
!
       if(iseq(icm)==1)then
       i1=dble(ranmar())*nptsize+1 ; tmq=tmp/dble(i1)
       xxi=ranmar()+0.3
!      either from first bank or from  bank
       if(ranmar() < 0.10)then
       j=dble(ranmar())*npop1+1 ; qosi(:,1)=posi1(:,j)
       jmutation=-j
!
!
                          else
       j=dble(ranmar())*npop+1 ; qosi(:,1)=posi(:,j)
       jmutation=j
                          endif
       qosi0(:)=qosi(:,1)
       if(ranmar() < 0.05)then
       qosi0(:)=qosi(:,1)
       do i=1,natom
       qosi0(3*(i-1)+1)=qosi(3*(i-1)+1,1)+(ranmar()-0.5)
       qosi0(3*(i-1)+2)=qosi(3*(i-1)+2,1)+(ranmar()-0.5)
       qosi0(3*(i-1)+3)=qosi(3*(i-1)+3,1)+(ranmar()-0.5)
       enddo
                          else
       i1=dble(ranmar())*nptsize+1 ; tmq=tmp/dble(i1)
       if(.not. lpbc)then
       do i=1,natom
       qosi0(3*(i-1)+1)=qosi(3*(i-1)+1,1)+(ranmar()-0.5)*tmq
       qosi0(3*(i-1)+2)=qosi(3*(i-1)+2,1)+(ranmar()-0.5)*tmq
       qosi0(3*(i-1)+3)=qosi(3*(i-1)+3,1)+(ranmar()-0.5)*tmq
       enddo
       if(ranmar() < 0.05)then
       j=dble(ranmar())*3+1
       vtest=1.d19
       do i=1,natom
       if(vtest > qosi0(3*(i-1)+j)) vtest=qosi0(3*(i-1)+j)
       enddo
       xxr=vtest
       vtest=-1.d19
       do i=1,natom
       if(vtest < qosi0(3*(i-1)+j)) vtest=qosi0(3*(i-1)+j)
       enddo
       xxr=vtest-xxr
!      xxr=xxr*(1.d0+ranmar()-0.5)
       do i=1,natom
       qosi0(3*(i-1)+j)=qosi(3*(i-1)+j,1)+(tmq)*sin(xxi*3.141592d0*qosi(3*(i-1)+j,1)/xxr)
       enddo
                          endif
                     endif
       if(lpbc)then
       do i=1,6
       t6(i)=qosi0(ish+i)
       enddo
       call latmat(t6,cmatrix,1)
       a1(:)=cmatrix(1,:) ; a2(:)=cmatrix(2,:) ; a3(:)=cmatrix(3,:)
       cellvol0=(cmatrix(1,2)*cmatrix(2,3)-cmatrix(1,3)*cmatrix(2,2))*cmatrix(3,1) &
               +(cmatrix(1,3)*cmatrix(2,1)-cmatrix(1,1)*cmatrix(2,3))*cmatrix(3,2) &
               +(cmatrix(1,1)*cmatrix(2,2)-cmatrix(1,2)*cmatrix(2,1))*cmatrix(3,3)
       cellvol0=abs(cellvol0)
               endif
       if(lpbc)then
       dv1=sqrt(dot_product(a1,a1)) ; dv2=sqrt(dot_product(a2,a2)) ; dv3=sqrt(dot_product(a3,a3))
       ddv1=dv1/(1.d-2) ; ddv2=dv2/(1.d-2) ; ddv3=dv3/(1.d-2)
       shift00=1.d-3/(max(dv1,dv2,dv3))
               else
       shift00=shift
       dv1=1.d0 ; dv2=1.d0 ; dv3=1.d0
       ddv1=1.d0/(1.d-2) ; ddv2=1.d0/(1.d-2) ; ddv3=1.d0/(1.d-2)
               endif
       if(lpbc)then
       do i=1,natom
       qosi0(3*(i-1)+1)=qosi(3*(i-1)+1,1)+(ranmar()-0.5)*(tmq/dv1)
       qosi0(3*(i-1)+2)=qosi(3*(i-1)+2,1)+(ranmar()-0.5)*(tmq/dv2)
       qosi0(3*(i-1)+3)=qosi(3*(i-1)+3,1)+(ranmar()-0.5)*(tmq/dv3)
       enddo
       if(ranmar() < 0.05)then
       j=dble(ranmar())*3+1
       if(j ==1) xxr=dv1
       if(j ==2) xxr=dv2
       if(j ==3) xxr=dv3
       do i=1,natom
       qosi0(3*(i-1)+j)=qosi(3*(i-1)+j,1)+(tmq/xxr)*sin(xxi*3.141592653589793238d0*qosi(3*(i-1)+j,1))
       enddo
                          endif
               endif
                          endif
!
       if(nspecies >= 2)then
       if(natom > 2)then
       do i2=1,1+3*dble(ranmar())+int(natom/3)
  410  continue
       i=dble(ranmar())*natom+1
       j=dble(ranmar())*natom+1
       if(itype(i) == itype(j)) goto 410
       wec(1,1)=qosi0(3*(i-1)+1)
       wec(2,1)=qosi0(3*(i-1)+2)
       wec(3,1)=qosi0(3*(i-1)+3)
       wec(1,2)=qosi0(3*(j-1)+1)
       wec(2,2)=qosi0(3*(j-1)+2)
       wec(3,2)=qosi0(3*(j-1)+3)
       qosi0(3*(i-1)+1)=wec(1,2)
       qosi0(3*(i-1)+2)=wec(2,2)
       qosi0(3*(i-1)+3)=wec(3,2)
       qosi0(3*(j-1)+1)=wec(1,1)
       qosi0(3*(j-1)+2)=wec(2,1)
       qosi0(3*(j-1)+3)=wec(3,1)
       enddo
                    endif
                        endif
!
       if(ranmar() <  0.10)then
       if(nspecies > 1    )then
       do i=1,natom
       wrk2(i)=ranmar()
       enddo
       call sortnr(natom,wrk2,iwrk2)
       qosi(:,1)=qosi0(:)
       do i=1,natom
       j=iwrk2(i)
       qosi0(3*(i-1)+1)=qosi(3*(j-1)+1,1)
       qosi0(3*(i-1)+2)=qosi(3*(j-1)+2,1)
       qosi0(3*(i-1)+3)=qosi(3*(j-1)+3,1)
       enddo
                           endif
                           endif
!
       if(lvcs)then
       do i=1,6
       r6(i)=qosi0(ish+i)
       enddo
       jtr=0
       do 
       jtr=jtr+1
       if(jcase ==1)then
       sig0=1.d0 ; xl0=0.d0
       do i=1,6
       call gauss(sig0,xl0,xl)
       t6(i)=r6(i)+(0.04d0*r6(i))*xl
       enddo
!
       if(ranmar() > 0.7)then
       i=dble(ranmar())*3+1
  331  continue
       vtest=(0.80d0+ranmar()*0.40d0)
       if(vtest < 0.95 .or. vtest > 1.05)then
       t6(i)=r6(i)*vtest
                                         else
       goto 331
                                         endif
                         endif
!
       if(jtr >200)then
       vtest=refvol*(1.d0+voltol*(ranmar()-0.5)*2.d0)
       call latmatvol(t6,cmatrix,vtest)
       call latmat(t6,cmatrix,0)
                   endif
!
                    endif
       if(jcase ==2)then
       do i=1,6
       t6(i)=r6(i)+(0.25d0*r6(i))*(ranmar()-0.5)
       enddo
                    endif
       if(jcase ==3)then
!      vtest=cellvol0*(1.d0+voltol*(ranmar()-0.5)*2.d0)
       vtest=refvol*(1.d0+voltol*(ranmar()-0.5)*2.d0)
       call gen_lattice_matrix(amatrix,s6,vtest)
       islice=dble(ranmar())*5+1
       xmix=0.1d0*dble(islice)
       do i=1,6
       t6(i)=r6(i)*(1.d0-xmix)+s6(i)*xmix
       enddo
                    endif
       if(jcase ==4)then
       call latmat(r6,cmatrix,1)
       call lat_mutation(cmatrix)
       call latmat(t6,cmatrix,0)
       vtest=refvol*(1.d0+voltol*(ranmar()-0.5)*2.d0)
       call latmatvol(t6,cmatrix,vtest)
       call latmat(t6,cmatrix,0)
                    endif
       call latmat(t6,cmatrix,1)
!
       vtest=(cmatrix(1,2)*cmatrix(2,3)-cmatrix(1,3)*cmatrix(2,2))*cmatrix(3,1) &
            +(cmatrix(1,3)*cmatrix(2,1)-cmatrix(1,1)*cmatrix(2,3))*cmatrix(3,2) &
            +(cmatrix(1,1)*cmatrix(2,2)-cmatrix(1,2)*cmatrix(2,1))*cmatrix(3,3)
       vtest=abs(vtest)
       if(vtest > refvol*(1.d0-voltol) .and. vtest < refvol*(1.d0+voltol))then
       call check_lat(llattice,lcate,cmatrix)
                                                                          else
       llattice=.false.
                                                                          endif
       if(llattice) goto 555
       if(jtr > 300)then
       do i=1,6
       t6(i)=r6(i)
       enddo
       write(6,*) 'we have trouble with lattice variations 300'
       goto 555
                    endif
       enddo
 555   continue
       if(jtr > 100) write(6,'(i5,2x,a13)')  jtr,'jtr-   code=1'
       do i=1,6
       qosi0(ish+i)=t6(i)
       enddo
               endif
!      differential evolution type : mutation and crossover for position vectors
       if(ranmar() < -0.1               )then
       if(icm  > npop)                   then
       call onedffvltn(ndeg,npop,posi,qosi0)
       do i=1,6
       qosi0(ish+i)=t6(i)
       enddo
       if(lpbc)then
       do i=1,natom*3
       qosi0(i)=qosi0(i)-anint(qosi0(i))
       if(qosi0(i) <= 0.d0) qosi0(i)=qosi0(i)+1.d0
       enddo
               endif
                                         endif
                                         endif
!      soft mutation : assupmtion : already optimized fixed lattice parameters
       if(ranmar() <  0.7 .and. iobj == 0)then
       if(icm  > npop)                    then
       j=dble(ranmar())*npop+1
       do i=1,1
       j=min(dble(j),dble(ranmar())*npop+1)
       enddo
!
       jmutation=j
       qosi(:,1)=posi(:,j)
       qosi0(:)=qosi(:,1)
       if(lpbc)then
       call softmutation(ndeg,qosi0,j,amp,t6)
               else
       call softmutation1(ndeg,qosi0,j,amp)
               endif
       if(ranmar() < 0.10)then
       if(lpbc)then
!      if(lvcs)then
       qosi(:,1)=qosi0(:)
       ispgrp=0
       allocate(sigmamatrix0(nspecies,nspecies))
       sigmamatrix0=sigmamatrix*0.5d0
       call gen_latt_site(ispgrp,ndeg,nspecies,nelements,symbl,sigmamatrix0,voltol,refvol,qosi0,lpbc,lvcs,lflagls)
       deallocate(sigmamatrix0)
       if(.not. lflagls) qosi0(:)=qosi(:,1)
!              endif
               endif
                          endif
                                          endif
                                          endif
       if(jmutation <0)then
       xmix=energy_sorted1(-jmutation)
                       else
       xmix=energy_sorted(jmutation)
                       endif
       write(6,'(i5,1x,e22.12,1x,a18,1x,i6)') jmutation,xmix,'jmutation,jparents',icm
                       endif
!
       if(iseq(icm)==2)then
       if(icase == 1)then
       ii=dble(ranmar())*npop+1
  11   continue
       jj=dble(ranmar())*npop+1
       if( ii == jj) goto 11
                     endif
       if(icase == 2)then
       ii=npop/2 ; ii=-dble(ii)*log(ranmar()) ; if(ii <= 0 .or. ii > npop) ii=npop*dble(ranmar())+1
  22   continue
       jj=npop/2 ; jj=-dble(jj)*log(ranmar()) ; if(jj <= 0 .or. jj > npop) jj=npop*dble(ranmar())+1
       if( ii == jj) goto 22
                     endif
       if(icase == 3)then
       ii=min(dble(ranmar())*npop+1,dble(ranmar())*npop+1)
  33   continue
       jj=min(dble(ranmar())*npop+1,dble(ranmar())*npop+1)
       if( ii == jj) goto 33
                     endif
       qosi(:,ii)=posi(:,ii)
       qosi(:,jj)=posi(:,jj)
       jparents(1)=ii
       jparents(2)=jj
!      either from first bank or from  bank
       if(ranmar() < 0.10)then
       j=dble(ranmar())*npop1+1
       qosi(:,jj)=posi1(:,j)
!      jparents(2)=-jj
       jparents(2)=-j
                          endif
       if(lpbc)then
       call tocarx(qosi(1,ii))
       call tocarx(qosi(1,jj))
               endif
       call centering(qosi(1,ii))
       call centering(qosi(1,jj))
!
       if(natom > 20)then
       call gen_randrot(rx,ry,rz)
       tt=matmul(ry,rz) ; ss=matmul(rx,tt)
       do i=1,natom
       vec(1)=qosi(3*(i-1)+1,ii)
       vec(2)=qosi(3*(i-1)+2,ii)
       vec(3)=qosi(3*(i-1)+3,ii)
       uec=matmul(ss,vec)
       qosi(3*(i-1)+1,ii)=uec(1)
       qosi(3*(i-1)+2,ii)=uec(2)
       qosi(3*(i-1)+3,ii)=uec(3)
       enddo
                     endif
!
       call gen_randrot(rx,ry,rz)
       tt=matmul(ry,rz) ; ss=matmul(rx,tt)
       do i=1,natom
       vec(1)=qosi(3*(i-1)+1,jj)
       vec(2)=qosi(3*(i-1)+2,jj)
       vec(3)=qosi(3*(i-1)+3,jj)
       uec=matmul(ss,vec)
       qosi(3*(i-1)+1,jj)=uec(1)
       qosi(3*(i-1)+2,jj)=uec(2)
       qosi(3*(i-1)+3,jj)=uec(3)
       enddo
!
       qosi0(:)=qosi(:,ii)
       do i=1,natom
       wrk2(i)=qosi0(3*(i-1)+1)
       enddo
       call sortnr(natom,wrk2,iwrk2)
       nb=0
       do i1=1,nspecies
       do i=1,natom
       i2=iwrk2(i)
       if(itype(i2) == i1)then
       nb=nb+1
       qosi(3*(nb-1)+1,ii)=qosi0(3*(i2-1)+1)
       qosi(3*(nb-1)+2,ii)=qosi0(3*(i2-1)+2)
       qosi(3*(nb-1)+3,ii)=qosi0(3*(i2-1)+3)
                          endif
       enddo
       enddo
       if(lpbc)then
       do i=1,6
       t6(i)=qosi0(ish+i)
       enddo
       call latmat(t6,cmatrix,1)
       a1(:)=cmatrix(1,:) ; a2(:)=cmatrix(2,:) ; a3(:)=cmatrix(3,:)
       cellvol0=(cmatrix(1,2)*cmatrix(2,3)-cmatrix(1,3)*cmatrix(2,2))*cmatrix(3,1) &
               +(cmatrix(1,3)*cmatrix(2,1)-cmatrix(1,1)*cmatrix(2,3))*cmatrix(3,2) &
               +(cmatrix(1,1)*cmatrix(2,2)-cmatrix(1,2)*cmatrix(2,1))*cmatrix(3,3)
       cellvol0=abs(cellvol0)
               endif
       if(lpbc)then
       dv1=sqrt(dot_product(a1,a1)) ; dv2=sqrt(dot_product(a2,a2)) ; dv3=sqrt(dot_product(a3,a3))
       ddv1=dv1/(1.d-2) ; ddv2=dv2/(1.d-2) ; ddv3=dv3/(1.d-2)
       shift00=1.d-3/(max(dv1,dv2,dv3))
               else
       shift00=shift
       dv1=1.d0 ; dv2=1.d0 ; dv3=1.d0
       ddv1=1.d0/(1.d-2) ; ddv2=1.d0/(1.d-2) ; ddv3=1.d0/(1.d-2)
               endif
!
       qosi0(:)=qosi(:,jj)
       do i=1,natom
       wrk2(i)=qosi0(3*(i-1)+1)
       enddo
       call sortnr(natom,wrk2,iwrk2)
       nb=0
       do i1=1,nspecies
       do i=1,natom
       i2=iwrk2(i)
       if(itype(i2) == i1)then
       nb=nb+1
       qosi(3*(nb-1)+1,jj)=qosi0(3*(i2-1)+1)
       qosi(3*(nb-1)+2,jj)=qosi0(3*(i2-1)+2)
       qosi(3*(nb-1)+3,jj)=qosi0(3*(i2-1)+3)
                          endif
       enddo
       enddo
!      additional rotation around x axis
       call gen_randrot(rx,ry,rz)
       do i=1,natom
       vec(1)=qosi(3*(i-1)+1,jj)
       vec(2)=qosi(3*(i-1)+2,jj)
       vec(3)=qosi(3*(i-1)+3,jj)
       uec=matmul(rx,vec)
       qosi(3*(i-1)+1,jj)=uec(1)
       qosi(3*(i-1)+2,jj)=uec(2)
       qosi(3*(i-1)+3,jj)=uec(3)
       enddo
!
       nb=0
       do ity=1,nspecies
       jfrac=dble(ranmar())*(nfrac-1)+1
       islice=(dble(jfrac)/dble(nfrac)) *dble(nelements(ity))
       if(islice <= 1) islice=1 ; if(islice >= nelements(ity)-1) islice=nelements(ity)-1
       i2=0
       do i=1,islice
       nb=nb+1
       i2=i2+1
       qosi0(3*(nb-1)+1)=qosi(3*(nb-1)+1,ii)+(ranmar()-0.5)/ddv1-shift00*0.5d0
       qosi0(3*(nb-1)+2)=qosi(3*(nb-1)+2,ii)+(ranmar()-0.5)/ddv2
       qosi0(3*(nb-1)+3)=qosi(3*(nb-1)+3,ii)+(ranmar()-0.5)/ddv3
       enddo
       do i=islice+1,nelements(ity)
       nb=nb+1
       i2=i2+1
       qosi0(3*(nb-1)+1)=qosi(3*(nb-1)+1,jj)+(ranmar()-0.5)/ddv1+shift00*0.5d0
       qosi0(3*(nb-1)+2)=qosi(3*(nb-1)+2,jj)+(ranmar()-0.5)/ddv2
       qosi0(3*(nb-1)+3)=qosi(3*(nb-1)+3,jj)+(ranmar()-0.5)/ddv3
       enddo
       if(i2 /= nelements(ity))then
       write(6,*) 'something went wrong 2-1'
                               stop
                               endif
       enddo
       if(nb /= natom)then
       write(6,*) 'something went wrong 2'
                      stop
                      endif
!
       if(lvcs)then
       jtr=0
       do 
       jtr=jtr+1
       islice=dble(ranmar())*5+1
       if(islice <= 1) islice=1 ; if(islice >= 5) islice=5
       nb=0
       do i=1,islice
       nb=nb+1
       qosi0(ish+i)=qosi(ish+i,ii)
       enddo
       do i=islice+1,6
       nb=nb+1
       qosi0(ish+i)=qosi(ish+i,jj)
       enddo
       if(nb /= 6)then
       write(6,*) 'something went wrong 2-vcs'
                  stop
                  endif
       if(jtr > 100)then
       sig0=1.d0 ; xl0=0.d0
       do i=1,6
       call gauss(sig0,xl0,xl)
       if(ranmar() > 0.5)then
       qosi0(ish+i)=qosi(ish+i,ii)+(0.04d0*qosi(ish+i,ii))*xl
                         else
       qosi0(ish+i)=qosi(ish+i,jj)+(0.04d0*qosi(ish+i,jj))*xl
                         endif
       enddo
                    endif
       if(jtr > 200)then
       sig0=1.d0 ; xl0=0.d0
       do i=1,6
       call gauss(sig0,xl0,xl)
       qosi0(ish+i)=qosi(ish+i,ii)+(0.04d0*qosi(ish+i,ii))*xl
       enddo
                    endif
       if(jtr > 300)then
       sig0=1.d0 ; xl0=0.d0
       do i=1,6
       call gauss(sig0,xl0,xl)
       qosi0(ish+i)=qosi(ish+i,jj)+(0.01d0*qosi(ish+i,jj))*xl
       if(ranmar() > 0.5) qosi0(ish+i)=qosi(ish+i,ii)+(0.01d0*qosi(ish+i,ii))*xl
       enddo
                    endif
       do i=1,6
       t6(i)=qosi0(ish+i)
       enddo
       call latmat(t6,cmatrix,1)
       vtest=(cmatrix(1,2)*cmatrix(2,3)-cmatrix(1,3)*cmatrix(2,2))*cmatrix(3,1) &
            +(cmatrix(1,3)*cmatrix(2,1)-cmatrix(1,1)*cmatrix(2,3))*cmatrix(3,2) &
            +(cmatrix(1,1)*cmatrix(2,2)-cmatrix(1,2)*cmatrix(2,1))*cmatrix(3,3)
       vtest=abs(vtest)
       if(vtest > refvol*(1.d0-voltol) .and. vtest < refvol*(1.d0+voltol))then
       call check_lat(llattice,lcate,cmatrix)
                                                                          else
       llattice=.false.
                                                                          endif
       if(llattice) goto 666
       if(jtr > 400)then
       do i=1,6
       t6(i)=qosi(ish+i,ii)
       enddo
       vtest=refvol*(1.d0+voltol*(ranmar()-0.5)*2.d0)
       call gen_lattice_matrix(cmatrix,t6,vtest)
       goto 666
                    endif
       enddo
 666   continue
       if(jtr > 400) write(6,*) 'we have trouble with ii,jj 400 : new lattice vector'
       if(jtr > 100) write(6,'(i5,2x,a13)')  jtr,'jtr    code=2'
       do i=1,6
       qosi0(ish+i)=t6(i)
       enddo
               endif
       if(lpbc) call tolatx(qosi0)
!
       tmp=energy_sorted(jparents(1))
       if(jparents(2) <0)then
       xmix=energy_sorted1(-jparents(2))
                         else
       xmix=energy_sorted(jparents(2))
                         endif
       write(6,'(2i5,1x,2e22.12,1x,a8,1x,i6)') jparents(1),jparents(2),tmp,xmix,'jparents',icm
                       endif
!
       if(iseq(icm)==3)then
       ii=dble(ranmar())*npop+1
       jj=dble(ranmar())*npop+1
       kk=dble(ranmar())*npop+1
       qosi(:,ii)=posi(:,ii)
       qosi(:,jj)=posi(:,jj)
       qosi(:,kk)=posi(:,kk)
!      either from first bank or from  bank
       if(ranmar() < 0.10)then
       j=dble(ranmar())*npop1+1
       qosi(:,jj)=posi1(:,j)
                          endif
!      either from first bank or from  bank
       if(ranmar() < 0.10)then
       j=dble(ranmar())*npop1+1
       qosi(:,kk)=posi1(:,j)
                          endif
       if(lpbc)then
       call tocarx(qosi(1,ii))
       call tocarx(qosi(1,jj))
       call tocarx(qosi(1,kk))
               endif
       call centering(qosi(1,ii))
       call centering(qosi(1,jj))
       call centering(qosi(1,kk))
!
       call gen_randrot(rx,ry,rz)
       do i=1,natom
       vec(1)=qosi(3*(i-1)+1,jj)
       vec(2)=qosi(3*(i-1)+2,jj)
       vec(3)=qosi(3*(i-1)+3,jj)
       uec=matmul(rz,vec)
       vec=matmul(ry,uec)
       uec=matmul(rx,vec)
       qosi(3*(i-1)+1,jj)=uec(1)
       qosi(3*(i-1)+2,jj)=uec(2)
       qosi(3*(i-1)+3,jj)=uec(3)
       enddo
!
       call gen_randrot(rx,ry,rz)
       do i=1,natom
       vec(1)=qosi(3*(i-1)+1,kk)
       vec(2)=qosi(3*(i-1)+2,kk)
       vec(3)=qosi(3*(i-1)+3,kk)
       uec=matmul(rz,vec)
       vec=matmul(ry,uec)
       uec=matmul(rx,vec)
       qosi(3*(i-1)+1,kk)=uec(1)
       qosi(3*(i-1)+2,kk)=uec(2)
       qosi(3*(i-1)+3,kk)=uec(3)
       enddo
!
       qosi0(:)=qosi(:,ii)
       do i=1,natom
       wrk2(i)=qosi0(3*(i-1)+1)
       enddo
       call sortnr(natom,wrk2,iwrk2)
       nb=0
       do i1=1,nspecies
       do i=1,natom
       i2=iwrk2(i)
       if(itype(i2) == i1)then
       nb=nb+1
       qosi(3*(nb-1)+1,ii)=qosi0(3*(i2-1)+1)
       qosi(3*(nb-1)+2,ii)=qosi0(3*(i2-1)+2)
       qosi(3*(nb-1)+3,ii)=qosi0(3*(i2-1)+3)
                          endif
       enddo
       enddo
       if(lpbc)then
       do i=1,6
       t6(i)=qosi0(ish+i)
       enddo
       call latmat(t6,cmatrix,1)
       a1(:)=cmatrix(1,:) ; a2(:)=cmatrix(2,:) ; a3(:)=cmatrix(3,:)
       cellvol0=(cmatrix(1,2)*cmatrix(2,3)-cmatrix(1,3)*cmatrix(2,2))*cmatrix(3,1) &
               +(cmatrix(1,3)*cmatrix(2,1)-cmatrix(1,1)*cmatrix(2,3))*cmatrix(3,2) &
               +(cmatrix(1,1)*cmatrix(2,2)-cmatrix(1,2)*cmatrix(2,1))*cmatrix(3,3)
       cellvol0=abs(cellvol0)
               endif
       if(lpbc)then
       dv1=sqrt(dot_product(a1,a1)) ; dv2=sqrt(dot_product(a2,a2)) ; dv3=sqrt(dot_product(a3,a3))
       ddv1=dv1/(1.d-2) ; ddv2=dv2/(1.d-2) ; ddv3=dv3/(1.d-2)
       shift00=1.d-3/(max(dv1,dv2,dv3))
               else
       shift00=shift
       dv1=1.d0 ; dv2=1.d0 ; dv3=1.d0
       ddv1=1.d0/(1.d-2) ; ddv2=1.d0/(1.d-2) ; ddv3=1.d0/(1.d-2)
               endif
!
       qosi0(:)=qosi(:,jj)
       do i=1,natom
       wrk2(i)=qosi0(3*(i-1)+1)
       enddo
       call sortnr(natom,wrk2,iwrk2)
       nb=0
       do i1=1,nspecies
       do i=1,natom
       i2=iwrk2(i)
       if(itype(i2) == i1)then
       nb=nb+1
       qosi(3*(nb-1)+1,jj)=qosi0(3*(i2-1)+1)
       qosi(3*(nb-1)+2,jj)=qosi0(3*(i2-1)+2)
       qosi(3*(nb-1)+3,jj)=qosi0(3*(i2-1)+3)
                          endif
       enddo
       enddo
!
       qosi0(:)=qosi(:,kk)
       do i=1,natom
       wrk2(i)=qosi0(3*(i-1)+1)
       enddo
       call sortnr(natom,wrk2,iwrk2)
!      do i=1,natom
!      i1=iwrk2(i)
!      qosi(3*(i-1)+1,kk)=qosi0(3*(i1-1)+1)
!      qosi(3*(i-1)+2,kk)=qosi0(3*(i1-1)+2)
!      qosi(3*(i-1)+3,kk)=qosi0(3*(i1-1)+3)
!      enddo
       nb=0
       do i1=1,nspecies
       do i=1,natom
       i2=iwrk2(i)
       if(itype(i2) == i1)then
       nb=nb+1
       qosi(3*(nb-1)+1,kk)=qosi0(3*(i2-1)+1)
       qosi(3*(nb-1)+2,kk)=qosi0(3*(i2-1)+2)
       qosi(3*(nb-1)+3,kk)=qosi0(3*(i2-1)+3)
                          endif
       enddo
       enddo
!      additional rotation around x axis
       call gen_randrot(rx,ry,rz)
       do i=1,natom
       vec(1)=qosi(3*(i-1)+1,jj)
       vec(2)=qosi(3*(i-1)+2,jj)
       vec(3)=qosi(3*(i-1)+3,jj)
       uec=matmul(rx,vec)
       qosi(3*(i-1)+1,jj)=uec(1)
       qosi(3*(i-1)+2,jj)=uec(2)
       qosi(3*(i-1)+3,jj)=uec(3)
       enddo
!
       nb=0
       do ity=1,nspecies
       i1=dble(ranmar())*(nfrac-1)+1
       i=dble(ranmar())*(nfrac-1)+1
       jslice=(dble(i1)/dble(nfrac)) *dble(nelements(ity))
       kslice=(dble(i)/dble(nfrac)) *dble(nelements(ity))
       if(jslice <= 1) jslice=1 ; if(jslice >= nelements(ity)-1) jslice=nelements(ity)-1
       if(kslice <= 1) kslice=1 ; if(kslice >= nelements(ity)-1) kslice=nelements(ity)-1
       i=jslice
       i1=kslice
       jslice=min(i1,i)
       kslice=max(i1,i)
       i2=0
       do i=1,jslice
       nb=nb+1
       i2=i2+1
       qosi0(3*(nb-1)+1)=qosi(3*(nb-1)+1,ii)+(ranmar()-0.5)/ddv1-shift00*0.5d0
       qosi0(3*(nb-1)+2)=qosi(3*(nb-1)+2,ii)+(ranmar()-0.5)/ddv2
       qosi0(3*(nb-1)+3)=qosi(3*(nb-1)+3,ii)+(ranmar()-0.5)/ddv3
       enddo
       do i=jslice+1,kslice
       nb=nb+1
       i2=i2+1
       qosi0(3*(nb-1)+1)=qosi(3*(nb-1)+1,jj)+(ranmar()-0.5)/ddv1
       qosi0(3*(nb-1)+2)=qosi(3*(nb-1)+2,jj)+(ranmar()-0.5)/ddv2
       qosi0(3*(nb-1)+3)=qosi(3*(nb-1)+3,jj)+(ranmar()-0.5)/ddv3
       enddo
       do i=kslice+1,nelements(ity)
       nb=nb+1
       i2=i2+1
       qosi0(3*(nb-1)+1)=qosi(3*(nb-1)+1,kk)+(ranmar()-0.5)/ddv1+shift00*0.5d0
       qosi0(3*(nb-1)+2)=qosi(3*(nb-1)+2,kk)+(ranmar()-0.5)/ddv2
       qosi0(3*(nb-1)+3)=qosi(3*(nb-1)+3,kk)+(ranmar()-0.5)/ddv3
       enddo
       if(i2 /= nelements(ity))then
       write(6,*) 'something went wrong 3-1'
                               stop
                               endif
       enddo
       if(nb /= natom)then
       write(6,*) 'something went wrong 3'
                      stop
                      endif
!
       if(lvcs)then
       jtr=0
       do 
       jtr=jtr+1
       jslice=dble(ranmar())*5+1
       kslice=dble(ranmar())*5+1
       if(jslice <= 1) jslice=1 ; if(jslice >= 5) jslice=5
       if(kslice <= 1) kslice=1 ; if(kslice >= 5) kslice=5
       i=jslice
       i1=kslice
       jslice=min(i1,i)
       kslice=max(i1,i)
       nb=0
       do i=1,jslice
       nb=nb+1
       qosi0(ish+i)=qosi(ish+i,ii)
       enddo
       do i=jslice+1,kslice
       nb=nb+1
       qosi0(ish+i)=qosi(ish+i,jj)
       enddo
       do i=kslice+1,6
       nb=nb+1
       qosi0(ish+i)=qosi(ish+i,kk)
       enddo
       if(nb /= 6)then
       write(6,*) 'something went wrong 3-vcs'
                  stop
                  endif
!
       if(jtr > 100)then
       sig0=1.d0 ; xl0=0.d0
       do i=1,6
       call gauss(sig0,xl0,xl)
       vtest=ranmar()
       if(vtest < 0.33333)                       qosi0(ish+i)=qosi(ish+i,ii)+(0.04d0*qosi(ish+i,ii))*xl
       if(vtest > 0.33333 .and. vtest < 0.66666) qosi0(ish+i)=qosi(ish+i,jj)+(0.04d0*qosi(ish+i,jj))*xl
       if(vtest > 0.66666)                       qosi0(ish+i)=qosi(ish+i,kk)+(0.04d0*qosi(ish+i,kk))*xl
       enddo
                    endif
!
       if(jtr > 200)then
       sig0=1.d0 ; xl0=0.d0
       do i=1,6
       call gauss(sig0,xl0,xl)
       qosi0(ish+i)=qosi(ish+i,ii)+(0.04d0*qosi(ish+i,ii))*xl
       enddo
                    endif
       do i=1,6
       t6(i)=qosi0(ish+i)
       enddo
       call latmat(t6,cmatrix,1)
       vtest=(cmatrix(1,2)*cmatrix(2,3)-cmatrix(1,3)*cmatrix(2,2))*cmatrix(3,1) &
            +(cmatrix(1,3)*cmatrix(2,1)-cmatrix(1,1)*cmatrix(2,3))*cmatrix(3,2) &
            +(cmatrix(1,1)*cmatrix(2,2)-cmatrix(1,2)*cmatrix(2,1))*cmatrix(3,3)
       vtest=abs(vtest)
       if(vtest > refvol*(1.d0-voltol) .and. vtest < refvol*(1.d0+voltol))then
       call check_lat(llattice,lcate,cmatrix)
                                                                          else
       llattice=.false.
                                                                          endif
       if(llattice) goto 777
       if(jtr > 400)then
       do i=1,6
       t6(i)=qosi(ish+i,ii)
       enddo
!      write(6,*) 'we have trouble with ii,jj,kk'
       goto 777
                    endif
       enddo
 777   continue
       if(jtr > 400) write(6,*) 'we have trouble with ii,jj,kk, 400'
       if(jtr > 100) write(6,'(i5,2x,a13)')  jtr,'jtr    code=3'
       do i=1,6
       qosi0(ish+i)=t6(i)
       enddo
               endif
       if(lpbc) call tolatx(qosi0)
                       endif
!
       if(iseq(icm)==4)then
       ii=dble(ranmar())*npop+1
       jj=dble(ranmar())*npop+1
       kk=dble(ranmar())*npop+1
       ll=dble(ranmar())*npop+1
       qosi(:,ii)=posi(:,ii)
       qosi(:,jj)=posi(:,jj)
       qosi(:,kk)=posi(:,kk)
       qosi(:,ll)=posi(:,ll)
!      either from first bank or from  bank
       if(ranmar() < 0.10)then
       j=dble(ranmar())*npop1+1
       qosi(:,jj)=posi1(:,j)
                          endif
!      either from first bank or from  bank
       if(ranmar() < 0.10)then
       j=dble(ranmar())*npop1+1
       qosi(:,kk)=posi1(:,j)
                          endif
!      either from first bank or from  bank
       if(ranmar() < 0.10)then
       j=dble(ranmar())*npop1+1
       qosi(:,ll)=posi1(:,j)
                          endif
       if(lpbc)then
       call tocarx(qosi(1,ii))
       call tocarx(qosi(1,jj))
       call tocarx(qosi(1,kk))
       call tocarx(qosi(1,ll))
               endif
       call centering(qosi(1,ii))
       call centering(qosi(1,jj))
       call centering(qosi(1,kk))
       call centering(qosi(1,ll))
!
       call gen_randrot(rx,ry,rz)
       do i=1,natom
       vec(1)=qosi(3*(i-1)+1,jj)
       vec(2)=qosi(3*(i-1)+2,jj)
       vec(3)=qosi(3*(i-1)+3,jj)
       uec=matmul(rz,vec)
       vec=matmul(ry,uec)
       uec=matmul(rx,vec)
       qosi(3*(i-1)+1,jj)=uec(1)
       qosi(3*(i-1)+2,jj)=uec(2)
       qosi(3*(i-1)+3,jj)=uec(3)
       enddo
!
       call gen_randrot(rx,ry,rz)
       do i=1,natom
       vec(1)=qosi(3*(i-1)+1,kk)
       vec(2)=qosi(3*(i-1)+2,kk)
       vec(3)=qosi(3*(i-1)+3,kk)
       uec=matmul(rz,vec)
       vec=matmul(ry,uec)
       uec=matmul(rx,vec)
       qosi(3*(i-1)+1,kk)=uec(1)
       qosi(3*(i-1)+2,kk)=uec(2)
       qosi(3*(i-1)+3,kk)=uec(3)
       enddo
!
       call gen_randrot(rx,ry,rz)
       do i=1,natom
       vec(1)=qosi(3*(i-1)+1,ll)
       vec(2)=qosi(3*(i-1)+2,ll)
       vec(3)=qosi(3*(i-1)+3,ll)
       uec=matmul(rz,vec)
       vec=matmul(ry,uec)
       uec=matmul(rx,vec)
       qosi(3*(i-1)+1,ll)=uec(1)
       qosi(3*(i-1)+2,ll)=uec(2)
       qosi(3*(i-1)+3,ll)=uec(3)
       enddo
!
       qosi0(:)=qosi(:,ii)
       do i=1,natom
       wrk2(i)=qosi0(3*(i-1)+1)
       enddo
       call sortnr(natom,wrk2,iwrk2)
       nb=0
       do i1=1,nspecies
       do i=1,natom
       i2=iwrk2(i)
       if(itype(i2) == i1)then
       nb=nb+1
       qosi(3*(nb-1)+1,ii)=qosi0(3*(i2-1)+1)
       qosi(3*(nb-1)+2,ii)=qosi0(3*(i2-1)+2)
       qosi(3*(nb-1)+3,ii)=qosi0(3*(i2-1)+3)
                          endif
       enddo
       enddo
       if(lpbc)then
       do i=1,6
       t6(i)=qosi0(ish+i)
       enddo
       call latmat(t6,cmatrix,1)
       a1(:)=cmatrix(1,:) ; a2(:)=cmatrix(2,:) ; a3(:)=cmatrix(3,:)
       cellvol0=(cmatrix(1,2)*cmatrix(2,3)-cmatrix(1,3)*cmatrix(2,2))*cmatrix(3,1) &
               +(cmatrix(1,3)*cmatrix(2,1)-cmatrix(1,1)*cmatrix(2,3))*cmatrix(3,2) &
               +(cmatrix(1,1)*cmatrix(2,2)-cmatrix(1,2)*cmatrix(2,1))*cmatrix(3,3)
       cellvol0=abs(cellvol0)
               endif
       if(lpbc)then
       dv1=sqrt(dot_product(a1,a1)) ; dv2=sqrt(dot_product(a2,a2)) ; dv3=sqrt(dot_product(a3,a3))
       ddv1=dv1/(1.d-2) ; ddv2=dv2/(1.d-2) ; ddv3=dv3/(1.d-2)
       shift00=1.d-3/(max(dv1,dv2,dv3))
               else
       shift00=shift
       dv1=1.d0 ; dv2=1.d0 ; dv3=1.d0
       ddv1=1.d0/(1.d-2) ; ddv2=1.d0/(1.d-2) ; ddv3=1.d0/(1.d-2)
               endif
!
       qosi0(:)=qosi(:,jj)
       do i=1,natom
       wrk2(i)=qosi0(3*(i-1)+1)
       enddo
       call sortnr(natom,wrk2,iwrk2)
       nb=0
       do i1=1,nspecies
       do i=1,natom
       i2=iwrk2(i)
       if(itype(i2) == i1)then
       nb=nb+1
       qosi(3*(nb-1)+1,jj)=qosi0(3*(i2-1)+1)
       qosi(3*(nb-1)+2,jj)=qosi0(3*(i2-1)+2)
       qosi(3*(nb-1)+3,jj)=qosi0(3*(i2-1)+3)
                          endif
       enddo
       enddo
!
       qosi0(:)=qosi(:,kk)
       do i=1,natom
       wrk2(i)=qosi0(3*(i-1)+1)
       enddo
       call sortnr(natom,wrk2,iwrk2)
       nb=0
       do i1=1,nspecies
       do i=1,natom
       i2=iwrk2(i)
       if(itype(i2) == i1)then
       nb=nb+1
       qosi(3*(nb-1)+1,kk)=qosi0(3*(i2-1)+1)
       qosi(3*(nb-1)+2,kk)=qosi0(3*(i2-1)+2)
       qosi(3*(nb-1)+3,kk)=qosi0(3*(i2-1)+3)
                          endif
       enddo
       enddo
!
       qosi0(:)=qosi(:,ll)
       do i=1,natom
       wrk2(i)=qosi0(3*(i-1)+1)
       enddo
       call sortnr(natom,wrk2,iwrk2)
       nb=0
       do i1=1,nspecies
       do i=1,natom
       i2=iwrk2(i)
       if(itype(i2) == i1)then
       nb=nb+1
       qosi(3*(nb-1)+1,ll)=qosi0(3*(i2-1)+1)
       qosi(3*(nb-1)+2,ll)=qosi0(3*(i2-1)+2)
       qosi(3*(nb-1)+3,ll)=qosi0(3*(i2-1)+3)
                          endif
       enddo
       enddo
!      additional rotation around x axis
       call gen_randrot(rx,ry,rz)
       do i=1,natom
       vec(1)=qosi(3*(i-1)+1,jj)
       vec(2)=qosi(3*(i-1)+2,jj)
       vec(3)=qosi(3*(i-1)+3,jj)
       uec=matmul(rx,vec)
       qosi(3*(i-1)+1,jj)=uec(1)
       qosi(3*(i-1)+2,jj)=uec(2)
       qosi(3*(i-1)+3,jj)=uec(3)
       enddo
!
       nb=0
       do ity=1,nspecies
       i=dble(ranmar())*(nfrac-1)+1
       jslice=(dble(i)/dble(nfrac)) *dble(nelements(ity))
       i=dble(ranmar())*(nfrac-1)+1
       kslice=(dble(i)/dble(nfrac)) *dble(nelements(ity))
       i=dble(ranmar())*(nfrac-1)+1
       lslice=(dble(i)/dble(nfrac)) *dble(nelements(ity))
       if(jslice <= 1) jslice=1 ; if(jslice >= nelements(ity)-1) jslice=nelements(ity)-1
       if(kslice <= 1) kslice=1 ; if(kslice >= nelements(ity)-1) kslice=nelements(ity)-1
       if(lslice <= 1) lslice=1 ; if(lslice >= nelements(ity)-1) lslice=nelements(ity)-1
       wrk1(1)=dble(jslice)
       wrk1(2)=dble(kslice)
       wrk1(3)=dble(lslice)
       call sortnr(3,wrk1,iwrk1)
       jslice=wrk1(iwrk1(1))
       kslice=wrk1(iwrk1(2))
       lslice=wrk1(iwrk1(3))
       i2=0
       do i=1,jslice
       nb=nb+1
       i2=i2+1
       qosi0(3*(nb-1)+1)=qosi(3*(nb-1)+1,ii)+(ranmar()-0.5)/ddv1-shift00*0.5d0
       qosi0(3*(nb-1)+2)=qosi(3*(nb-1)+2,ii)+(ranmar()-0.5)/ddv2
       qosi0(3*(nb-1)+3)=qosi(3*(nb-1)+3,ii)+(ranmar()-0.5)/ddv3
       enddo
       do i=jslice+1,kslice
       nb=nb+1
       i2=i2+1
       qosi0(3*(nb-1)+1)=qosi(3*(nb-1)+1,jj)+(ranmar()-0.5)/ddv1
       qosi0(3*(nb-1)+2)=qosi(3*(nb-1)+2,jj)+(ranmar()-0.5)/ddv2
       qosi0(3*(nb-1)+3)=qosi(3*(nb-1)+3,jj)+(ranmar()-0.5)/ddv3
       enddo
       do i=kslice+1,lslice
       nb=nb+1
       i2=i2+1
       qosi0(3*(nb-1)+1)=qosi(3*(nb-1)+1,kk)+(ranmar()-0.5)/ddv1
       qosi0(3*(nb-1)+2)=qosi(3*(nb-1)+2,kk)+(ranmar()-0.5)/ddv2
       qosi0(3*(nb-1)+3)=qosi(3*(nb-1)+3,kk)+(ranmar()-0.5)/ddv3
       enddo
       do i=lslice+1,nelements(ity)
       nb=nb+1
       i2=i2+1
       qosi0(3*(nb-1)+1)=qosi(3*(nb-1)+1,ll)+(ranmar()-0.5)/ddv1+shift00*0.5d0
       qosi0(3*(nb-1)+2)=qosi(3*(nb-1)+2,ll)+(ranmar()-0.5)/ddv2
       qosi0(3*(nb-1)+3)=qosi(3*(nb-1)+3,ll)+(ranmar()-0.5)/ddv3
       enddo
       if(i2 /= nelements(ity))then
       write(6,*) 'something went wrong 4-1'
                               stop
                               endif
       enddo
       if(nb /= natom)then
       write(6,*) 'something went wrong 4'
                      stop
                      endif
!
       if(lvcs)then
       jtr=0
       do 
       jtr=jtr+1
       jslice=dble(ranmar())*5+1
       kslice=dble(ranmar())*5+1
       lslice=dble(ranmar())*5+1
       if(jslice <= 1) jslice=1 ; if(jslice >= 5) jslice=5
       if(kslice <= 1) kslice=1 ; if(kslice >= 5) kslice=5
       if(lslice <= 1) lslice=1 ; if(lslice >= 5) lslice=5
       wrk1(1)=dble(jslice)
       wrk1(2)=dble(kslice)
       wrk1(3)=dble(lslice)
       call sortnr(3,wrk1,iwrk1)
       jslice=wrk1(iwrk1(1))
       kslice=wrk1(iwrk1(2))
       lslice=wrk1(iwrk1(3))
       nb=0
       do i=1,jslice
       nb=nb+1
       qosi0(ish+i)=qosi(ish+i,ii)
       enddo
       do i=jslice+1,kslice
       nb=nb+1
       qosi0(ish+i)=qosi(ish+i,jj)
       enddo
       do i=kslice+1,lslice
       nb=nb+1
       qosi0(ish+i)=qosi(ish+i,kk)
       enddo
       do i=lslice+1,6
       nb=nb+1
       qosi0(ish+i)=qosi(ish+i,ll)
       enddo
       if(nb /= 6)then
       write(6,*) 'something went wrong 4-vcs'
                  stop
                  endif
!
       if(jtr > 100)then
       sig0=1.d0 ; xl0=0.d0
       do i=1,6
       call gauss(sig0,xl0,xl)
       vtest=ranmar()
       if(vtest < 0.2500)                      qosi0(ish+i)=qosi(ish+i,ii)+(0.04d0*qosi(ish+i,ii))*xl
       if(vtest > 0.2500 .and. vtest < 0.5000) qosi0(ish+i)=qosi(ish+i,jj)+(0.04d0*qosi(ish+i,jj))*xl
       if(vtest > 0.5000 .and. vtest < 0.7500) qosi0(ish+i)=qosi(ish+i,kk)+(0.04d0*qosi(ish+i,kk))*xl
       if(vtest > 0.7500)                      qosi0(ish+i)=qosi(ish+i,ll)+(0.04d0*qosi(ish+i,ll))*xl
       enddo
                    endif
!
       if(jtr > 200)then
       sig0=1.d0 ; xl0=0.d0
       do i=1,6
       call gauss(sig0,xl0,xl)
       qosi0(ish+i)=qosi(ish+i,ii)+(0.04d0*qosi(ish+i,ii))*xl
       enddo
                    endif
       do i=1,6
       t6(i)=qosi0(ish+i)
       enddo
       call latmat(t6,cmatrix,1)
       vtest=(cmatrix(1,2)*cmatrix(2,3)-cmatrix(1,3)*cmatrix(2,2))*cmatrix(3,1) &
            +(cmatrix(1,3)*cmatrix(2,1)-cmatrix(1,1)*cmatrix(2,3))*cmatrix(3,2) &
            +(cmatrix(1,1)*cmatrix(2,2)-cmatrix(1,2)*cmatrix(2,1))*cmatrix(3,3)
       vtest=abs(vtest)
       if(vtest > refvol*(1.d0-voltol) .and. vtest < refvol*(1.d0+voltol))then
       call check_lat(llattice,lcate,cmatrix)
                                                                          else
       llattice=.false.
                                                                          endif
       if(llattice) goto 888
       if(jtr > 400)then
       do i=1,6
       t6(i)=qosi(ish+i,ii)
       enddo
!      write(6,*) 'we have trouble with ii,jj,kk,ll'
       goto 888
                    endif
       enddo
 888   continue
       if(jtr > 400) write(6,*) 'we have trouble with ii,jj,kk,ll 400'
       if(jtr > 100) write(6,'(i5,2x,a13)')  jtr,'jtr    code=4'
       do i=1,6
       qosi0(ish+i)=t6(i)
       enddo
               endif
       if(lpbc) call tolatx(qosi0)
                       endif
!      special process for low coordinated atoms and contacts
       call danglingbond_care()
       
       end subroutine csa_rnd_lattice_basis
!
!      Written by In-Ho Lee, KRISS, January 28, 2013.
       subroutine tocarx(qqq)
       USE csa_application, ONLY : natom
       implicit none
       real*8 qqq(ndeg)
       integer j,i,ish
       real*8 t6(6),cmatrix(3,3),a1(3),a2(3),a3(3),x,y,z

       ish=ndeg-6
       do i=1,6
       t6(i)=qqq(ish+i)
       enddo
       call latmat(t6,cmatrix,1)
       a1(:)=cmatrix(1,:) ; a2(:)=cmatrix(2,:) ; a3(:)=cmatrix(3,:)
!
       do j=1,natom
       x=a1(1)*qqq(3*(j-1)+1)+a2(1)*qqq(3*(j-1)+2)+a3(1)*qqq(3*(j-1)+3)
       y=a1(2)*qqq(3*(j-1)+1)+a2(2)*qqq(3*(j-1)+2)+a3(2)*qqq(3*(j-1)+3)
       z=a1(3)*qqq(3*(j-1)+1)+a2(3)*qqq(3*(j-1)+2)+a3(3)*qqq(3*(j-1)+3)
       qqq(3*(j-1)+1)=x
       qqq(3*(j-1)+2)=y
       qqq(3*(j-1)+3)=z
       enddo
       end subroutine tocarx
!
!      Written by In-Ho Lee, KRISS, January 28, 2013.
       subroutine tolatx(qqq)
       USE csa_application, ONLY : natom
       implicit none
       real*8 qqq(ndeg)
       real*8 b(3,3),devid
       integer j,i,ish
       real*8 t6(6),cmatrix(3,3),a1(3),a2(3),a3(3),d1,d2,d3

       ish=ndeg-6
       do i=1,6
       t6(i)=qqq(ish+i)
       enddo
       call latmat(t6,cmatrix,1)
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
       end subroutine tolatx

!      Written by In-Ho Lee, KRISS, January 28, 2013.
       subroutine lat_mutation(cmatrix)
       implicit none
       real*8 cmatrix(3,3)
       real*8 unitm(3,3),tmpmat(3,3),gaussm(3,3),strainedm(3,3),x,y
       integer i,j
       real ranmar

       tmpmat=cmatrix
       unitm=0.d0
       do i=1,3
       unitm(i,i)=1.d0
       enddo
       do i=1,3
       do j=1,3
       x=ranmar()
       y=ranmar()
       if(y > 0.5) gaussm(i,j)=exp(-x**2)*0.25
       if(y <=0.5) gaussm(i,j)=-exp(-x**2)*0.25
       if(i == j)  gaussm(i,j)=abs(gaussm(i,j))*3.0
       enddo
       enddo
       strainedm=unitm+gaussm
       cmatrix=matmul(tmpmat,strainedm)
       end subroutine lat_mutation
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine csa_bank_dump(idirection)
       implicit none
       integer idirection

       if(idirection == 0)then
       open(1,file='fort.1',form='formatted')
       write(1,*) ndeg,npop,npop1
       write(1,*) posi1,energy_sorted1
       write(1,*) posi,energy_sorted
       write(1,*) posi_best,energy_best,davg
       close(1)
                          endif
       if(idirection == 1)then
       open(1,file='fort.1',form='formatted')
       read(1,*) ndeg_r,npop_r,npop1_r
       if(ndeg_r /= ndeg)then
       write(6,*) 'system size mismatch'
       write(6,*) 'fort.1 and csa.in are different from each other'
       close(1)
       return
                         endif
       allocate(posi_r(ndeg,npop_r))
       allocate(posi1_r(ndeg,npop1_r))
       allocate(energy_sorted_r(npop_r))
       allocate(energy_sorted1_r(npop1_r))
       allocate(posi_best_r(ndeg))
       read(1,*) posi1_r,energy_sorted1_r
       read(1,*) posi_r,energy_sorted_r
       read(1,*) posi_best_r,energy_best_r,davg_r
       close(1)
                          endif
       end subroutine csa_bank_dump
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine csa_bank_sort(ibank)
       implicit none
       integer ibank
       integer j,j1

       if(ibank ==1)then
       wrk1(1:npop1)=energy_sorted1(1:npop1) ; qosi1=posi1 ; call sortnr(npop1,wrk1,iwrk1)
       do j=1,npop1
       j1=iwrk1(j)
       posi1(:,j)=qosi1(:,j1)
       energy_sorted1(j)=wrk1(j1)
       enddo
                    else
       wrk1(1:npop)=energy_sorted(1:npop) ; qosi=posi ; call sortnr(npop,wrk1,iwrk1)
       do j=1,npop
       j1=iwrk1(j)
       posi(:,j)=qosi(:,j1)
       energy_sorted(j)=wrk1(j1)
       enddo
                    endif
       end subroutine csa_bank_sort
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine danglingbond_care()
       USE csa_application, ONLY : natom,rc1,ncoord,wrk2,iwrk2,lpbc
       implicit none
       integer i,iord
       integer itgt,imvg,i3,ish
       real*8 dv1,dv2,dv3,cmatrix(3,3),t6(6),a1(3),a2(3),a3(3)
       real ranmar

       ish=ndeg-6
       if(lpbc)then
       do i=1,6
       t6(i)=qosi0(ish+i)
       enddo
       call latmat(t6,cmatrix,1)
       a1(:)=cmatrix(1,:) ; a2(:)=cmatrix(2,:) ; a3(:)=cmatrix(3,:)
       dv1=sqrt(dot_product(a1,a1))/1.1d0
       dv2=sqrt(dot_product(a2,a2))/1.1d0
       dv3=sqrt(dot_product(a3,a3))/1.1d0
               else
       dv1=0.5d0/rc1
       dv2=0.5d0/rc1
       dv3=0.5d0/rc1
               endif
       if(.not. lpbc)then
       call centering(qosi0)
       do iord=1,1
       call gen_coordination()
       do i=1,natom
       wrk2(i)=dble(ncoord(i))
       enddo
       call sortnr(natom,wrk2,iwrk2)
       do i=1,natom
       if(wrk2(iwrk2(i)) <  1)then
       qosi0(3*(iwrk2(i)-1)+1)=qosi0(3*(iwrk2(i)-1)+1)*0.1
       qosi0(3*(iwrk2(i)-1)+2)=qosi0(3*(iwrk2(i)-1)+2)*0.1
       qosi0(3*(iwrk2(i)-1)+3)=qosi0(3*(iwrk2(i)-1)+3)*0.1
                              endif
       enddo
       if(natom >2)then
       imvg=iwrk2(1) ; itgt=iwrk2(2) ; i3=iwrk2(3)
       qosi0(3*(imvg-1)+1)=qosi0(3*(i3-1)+1)+(ranmar()-0.5)/dv1
       qosi0(3*(imvg-1)+2)=qosi0(3*(i3-1)+2)+(ranmar()-0.5)/dv2
       qosi0(3*(imvg-1)+3)=qosi0(3*(i3-1)+3)+(ranmar()-0.5)/dv3
       qosi0(3*(itgt-1)+1)=qosi0(3*(i3-1)+1)+(ranmar()-0.5)/dv1
       qosi0(3*(itgt-1)+2)=qosi0(3*(i3-1)+2)+(ranmar()-0.5)/dv2
       qosi0(3*(itgt-1)+3)=qosi0(3*(i3-1)+3)+(ranmar()-0.5)/dv3
                   endif
       call repulsion_care()
       enddo
                     endif
       if(lpbc)then
       do iord=1,1
       call gen_coordination()
       do i=1,natom
       wrk2(i)=dble(ncoord(i))
       enddo
       call sortnr(natom,wrk2,iwrk2)
       imvg=iwrk2(1) ; itgt=iwrk2(2)
       qosi0(3*(imvg-1)+1)=qosi0(3*(itgt-1)+1)+(ranmar()-0.5)/dv1
       qosi0(3*(imvg-1)+2)=qosi0(3*(itgt-1)+2)+(ranmar()-0.5)/dv2
       qosi0(3*(imvg-1)+3)=qosi0(3*(itgt-1)+3)+(ranmar()-0.5)/dv3
       call repulsion_care()
       enddo
               endif
       end subroutine danglingbond_care
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine repulsion_care()
       USE csa_application, ONLY : natom,lpbc,sigmamatrix,nspecies,itype
       USE csa_application, ONLY : iwrk2,iwrk4
       implicit none
       real*8 dv1,dv2,dv3,repul0,repul1,repul,ddxx,cmatrix(3,3),t6(6),a1(3),a2(3),a3(3),cellvol0
       integer i,j,itr,ktr,ljitramax,ish
       integer iswp,jswp,itrx
       real*8 wec(3,2)
       integer nprint
       integer kcase
       real ranmar

       do i=1,natom*3
       call tonormal(qosi0(i))
       enddo
       nprint=0
       nprint=1
       repul0=0.d0
       repul1=0.d0
       do i=1,nspecies
       do j=i,nspecies
       repul0=repul0+1.d0
       repul1=repul1+sigmamatrix(i,j)
       enddo
       enddo
       repul1=repul1/repul0

       kcase=2
       kcase=1
       kcase=3

       ddxx=1.5d0 *repul1
       if(kcase ==1) ddxx=1.5d0 *repul1
       if(kcase ==2) ddxx=1.5d0 *sigmamatrix(1,1)


       ljitramax=100000*natom

!      ljitramax=0

       if(ljitramax <= 0) return
       iwrk2=1

       ish=ndeg-6
       if(lpbc)then
       do i=1,6
       t6(i)=qosi0(ish+i)
       enddo
       call latmat(t6,cmatrix,1)
       a1(:)=cmatrix(1,:) ; a2(:)=cmatrix(2,:) ; a3(:)=cmatrix(3,:)
       cellvol0=(cmatrix(1,2)*cmatrix(2,3)-cmatrix(1,3)*cmatrix(2,2))*cmatrix(3,1) &
               +(cmatrix(1,3)*cmatrix(2,1)-cmatrix(1,1)*cmatrix(2,3))*cmatrix(3,2) &
               +(cmatrix(1,1)*cmatrix(2,2)-cmatrix(1,2)*cmatrix(2,1))*cmatrix(3,3)
       cellvol0=abs(cellvol0)
       repul=(cellvol0/dble(natom))**(1.d0/3.d0)
       ddxx=repul/2.d0

       dv1=sqrt(dot_product(a1,a1))/ddxx
       dv2=sqrt(dot_product(a2,a2))/ddxx
       dv3=sqrt(dot_product(a3,a3))/ddxx
       call direct_pbc(qosi0)
       qosi00=qosi0
       call cal_repulsion(qosi00,repul0)
       repul=0.d0 ; itr=0
       repul1=repul0
       if(repul1 < 1.d-8) goto 101
       do itr=1,ljitramax
!
       if(nspecies >= 2)then
       if(natom > 2)then
       do itrx=1,1+3*dble(ranmar())+int(natom/3)
  510  continue
       iswp=dble(ranmar())*natom+1 ; jswp=dble(ranmar())*natom+1
       if(itype(iswp) == itype(jswp)) goto 510
       wec(1,1)=qosi0(3*(iswp-1)+1)
       wec(2,1)=qosi0(3*(iswp-1)+2)
       wec(3,1)=qosi0(3*(iswp-1)+3)
       wec(1,2)=qosi0(3*(jswp-1)+1)
       wec(2,2)=qosi0(3*(jswp-1)+2)
       wec(3,2)=qosi0(3*(jswp-1)+3)
       qosi0(3*(iswp-1)+1)=wec(1,2)
       qosi0(3*(iswp-1)+2)=wec(2,2)
       qosi0(3*(iswp-1)+3)=wec(3,2)
       qosi0(3*(jswp-1)+1)=wec(1,1)
       qosi0(3*(jswp-1)+2)=wec(2,1)
       qosi0(3*(jswp-1)+3)=wec(3,1)
       enddo
                    endif
                        endif
!
       do ktr=1,1000
!
       do i=1,natom
       if(iwrk2(i) ==1)then
       qosi0(3*(i-1)+1)=qosi00(3*(i-1)+1)+(ranmar()-0.5)/dv1
       qosi0(3*(i-1)+2)=qosi00(3*(i-1)+2)+(ranmar()-0.5)/dv2
       qosi0(3*(i-1)+3)=qosi00(3*(i-1)+3)+(ranmar()-0.5)/dv3
       j=iwrk4(i)
       qosi0(3*(j-1)+1)=qosi00(3*(j-1)+1)+(ranmar()-0.5)/dv1
       qosi0(3*(j-1)+2)=qosi00(3*(j-1)+2)+(ranmar()-0.5)/dv2
       qosi0(3*(j-1)+3)=qosi00(3*(j-1)+3)+(ranmar()-0.5)/dv3
                       endif
       enddo
       if(ranmar() < 0.05)then
       do i=1,natom
       qosi0(3*(i-1)+1)=qosi00(3*(i-1)+1)+(ranmar()-0.5)/dv1
       qosi0(3*(i-1)+2)=qosi00(3*(i-1)+2)+(ranmar()-0.5)/dv2
       qosi0(3*(i-1)+3)=qosi00(3*(i-1)+3)+(ranmar()-0.5)/dv3
       enddo
                          endif
       call direct_pbc(qosi0)
       call cal_repulsion(qosi0,repul)
       if(repul0 > repul)then
!      write(6,'(a7,2x,2e18.8)') 'updated', repul,repul0
       qosi00(:)=qosi0(:)
       repul0=repul
!      if(repul0 <1.d-8) exit
       if(repul0 <1.d-8) goto 101
                         endif
!
       enddo
       enddo
  101  continue
       if(nprint ==1)then
       write(6,'(i8,1x,a9,1x,2f16.1,2x,f18.8)') itr,'itr,repul',repul0,repul1,cellvol0
       if(repul > 1.d-8) write(6,'(i8,1x,a9,1x,4f18.8)') itr,'itr,repul',repul,repul0,repul1,cellvol0
                     endif
!      do i=1,natom
!      write(6,'(3f18.8)') qosi00(3*(i-1)+1),qosi00(3*(i-1)+2),qosi00(3*(i-1)+3)
!      enddo
       qosi0=qosi00
       call direct_pbc(qosi0)
               else
       call centering(qosi0)
       qosi00=qosi0
       call cal_repulsion(qosi00,repul0)
       repul=0.d0 ; itr=0
       repul1=repul0
       if(repul1 < 1.d-8) goto 102
       do itr=1,ljitramax
!
       if(nspecies >= 2)then
       if(natom > 2)then
       do itrx=1,1+3*dble(ranmar())+int(natom/3)
  610  continue
       iswp=dble(ranmar())*natom+1 ; jswp=dble(ranmar())*natom+1
       if(itype(iswp) == itype(jswp)) goto 610
       wec(1,1)=qosi0(3*(iswp-1)+1)
       wec(2,1)=qosi0(3*(iswp-1)+2)
       wec(3,1)=qosi0(3*(iswp-1)+3)
       wec(1,2)=qosi0(3*(jswp-1)+1)
       wec(2,2)=qosi0(3*(jswp-1)+2)
       wec(3,2)=qosi0(3*(jswp-1)+3)
       qosi0(3*(iswp-1)+1)=wec(1,2)
       qosi0(3*(iswp-1)+2)=wec(2,2)
       qosi0(3*(iswp-1)+3)=wec(3,2)
       qosi0(3*(jswp-1)+1)=wec(1,1)
       qosi0(3*(jswp-1)+2)=wec(2,1)
       qosi0(3*(jswp-1)+3)=wec(3,1)
       enddo
                    endif
                        endif
!
       do ktr=1,1000
!
       do i=1,natom
       if(iwrk2(i) == 1) then
       qosi0(3*(i-1)+1)=qosi00(3*(i-1)+1)+(ranmar()-0.5)*ddxx*(1.d0+ranmar()-0.5)
       qosi0(3*(i-1)+2)=qosi00(3*(i-1)+2)+(ranmar()-0.5)*ddxx*(1.d0+ranmar()-0.5)
       qosi0(3*(i-1)+3)=qosi00(3*(i-1)+3)+(ranmar()-0.5)*ddxx*(1.d0+ranmar()-0.5)
       j=iwrk4(i)
       qosi0(3*(j-1)+1)=qosi00(3*(j-1)+1)+(ranmar()-0.5)*ddxx*(1.d0+ranmar()-0.5)
       qosi0(3*(j-1)+2)=qosi00(3*(j-1)+2)+(ranmar()-0.5)*ddxx*(1.d0+ranmar()-0.5)
       qosi0(3*(j-1)+3)=qosi00(3*(j-1)+3)+(ranmar()-0.5)*ddxx*(1.d0+ranmar()-0.5)
                         endif
       enddo
       if(ranmar() < 0.05)then
       do i=1,natom
       qosi0(3*(i-1)+1)=qosi00(3*(i-1)+1)+(ranmar()-0.5)*ddxx*(1.d0+ranmar()-0.5)
       qosi0(3*(i-1)+2)=qosi00(3*(i-1)+2)+(ranmar()-0.5)*ddxx*(1.d0+ranmar()-0.5)
       qosi0(3*(i-1)+3)=qosi00(3*(i-1)+3)+(ranmar()-0.5)*ddxx*(1.d0+ranmar()-0.5)
       enddo
                          endif
       call centering(qosi0)
       call cal_repulsion(qosi0,repul)
       if(repul0 > repul)then
!      write(6,'(a7,2x,2e18.8)') 'updated', repul,repul0
       qosi00(:)=qosi0(:)
       repul0=repul
!      if(repul0 <1.d-8) exit
       if(repul0 <1.d-8) goto 102
                         endif
!
       enddo
       enddo
  102  continue
       if(nprint ==1)then
       write(6,'(i8,1x,a9,1x,2f16.1)') itr,'itr,repul',repul0,repul1
       if(repul > 1.d-8) write(6,'(i8,1x,a9,1x,3f18.8)') itr,'itr,repul',repul,repul0,repul1
                     endif
!      do i=1,natom
!      write(6,'(3f18.8)') qosi00(3*(i-1)+1),qosi00(3*(i-1)+2),qosi00(3*(i-1)+3)
!      enddo
       qosi0=qosi00
       call centering(qosi0)
               endif
       end subroutine repulsion_care
!       
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine cal_repulsion(qqq,repul)
       USE csa_application, ONLY : natom,itype,sigmamatrix,lpbc
       USE csa_application, ONLY : iwrk2,iwrk4
       implicit none
       real*8 repul,qqq(ndeg)
       integer i,j,ish
       real*8 x,y,z,r,d1,d2,d3,ctest,sig,rho,rs,pi,a1(3),a2(3),a3(3),cmatrix(3,3),t6(6)

       iwrk2=0
       if(lpbc)then
       ish=ndeg-6
       do i=1,6
       t6(i)=qqq(ish+i)
       enddo
       call latmat(t6,cmatrix,1)
       a1(:)=cmatrix(1,:) ; a2(:)=cmatrix(2,:) ; a3(:)=cmatrix(3,:)
               endif
       sig=1.5d0
       if(lpbc)then
       ctest=a1(1)*a2(2)*a3(3)-a1(2)*a2(1)*a3(3)-a1(1)*a2(3)*a3(2)   &
            +a1(3)*a2(1)*a3(2)+a1(2)*a2(3)*a3(1)-a1(3)*a2(2)*a3(1)
       ctest=abs(ctest)
       pi=4.d0*atan(1.d0)
       rho=dble(natom)/ctest
       rs=(3.d0/(rho*4.d0*pi))**(1.d0/3.d0)
       sig=(rs*2.d0)*(2.d0)**(-1.d0/6.d0)
               endif
       repul=0.d0
       do i=1,natom-1
       do j=i+1,natom
       if(lpbc)then
       d1=qqq(3*(i-1)+1)-qqq(3*(j-1)+1)
       d2=qqq(3*(i-1)+2)-qqq(3*(j-1)+2)
       d3=qqq(3*(i-1)+3)-qqq(3*(j-1)+3)
       d1=d1-anint(d1)
       d2=d2-anint(d2)
       d3=d3-anint(d3)
       x=d1*a1(1)+d2*a2(1)+d3*a3(1)
       y=d1*a1(2)+d2*a2(2)+d3*a3(2)
       z=d1*a1(3)+d2*a2(3)+d3*a3(3)
               else
       x=qqq(3*(j-1)+1)-qqq(3*(i-1)+1)
       y=qqq(3*(j-1)+2)-qqq(3*(i-1)+2)
       z=qqq(3*(j-1)+3)-qqq(3*(i-1)+3)
               endif
       r=sqrt(x*x+y*y+z*z)
       if(r < 1.d-8) r=1.d-8
       sig=sigmamatrix(itype(i),itype(j))
       if( r <  sig)then
       repul=repul+1.d0
       iwrk2(i)=1
       iwrk2(j)=1
       iwrk4(i)=j
       iwrk4(j)=i
                    endif
       enddo
       enddo
       end subroutine cal_repulsion
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine gen_coordination()
       USE csa_application, ONLY : natom,rc1,ncoord,lpbc
       implicit none
       integer i,j,kk1,ish
       real*8 x,y,z,r,d1,d2,d3,t6(6),cmatrix(3,3),a1(3),a2(3),a3(3)

       if(lpbc)then
       ish=ndeg-6
       do i=1,6
       t6(i)=qosi0(ish+i)
       enddo
       call latmat(t6,cmatrix,1)
       a1(:)=cmatrix(1,:) ; a2(:)=cmatrix(2,:) ; a3(:)=cmatrix(3,:)
               endif
       do i=1,natom
       kk1=0
       do j=1,natom
       if(j == i) cycle
       if(lpbc)then
       d1=qosi0(3*(i-1)+1)-qosi0(3*(j-1)+1)
       d2=qosi0(3*(i-1)+2)-qosi0(3*(j-1)+2)
       d3=qosi0(3*(i-1)+3)-qosi0(3*(j-1)+3)
       d1=d1-anint(d1)
       d2=d2-anint(d2)
       d3=d3-anint(d3)
       x=d1*a1(1)+d2*a2(1)+d3*a3(1)
       y=d1*a1(2)+d2*a2(2)+d3*a3(2)
       z=d1*a1(3)+d2*a2(3)+d3*a3(3)
               else
       x=qosi0(3*(i-1)+1)-qosi0(3*(j-1)+1)
       y=qosi0(3*(i-1)+2)-qosi0(3*(j-1)+2)
       z=qosi0(3*(i-1)+3)-qosi0(3*(j-1)+3)
               endif
       r=sqrt(x*x+y*y+z*z)
       if( r <= rc1)then
       kk1=kk1+1
                    endif
       enddo
       ncoord(i)=kk1
       enddo
       end subroutine gen_coordination
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine gen_lattice_matrix(amatrix,s6,cellvol0)
       implicit none
       real*8 amatrix(3,3),s6(6),cellvol0
       integer isgindex,jtr
       logical llattice,lcate
       real ranmar

       lcate=.false.
       jtr=0
       do 
  111  continue
       jtr=jtr+1
       isgindex=dble(ranmar())*230+1
       call gen_sg_lat(isgindex,cellvol0,amatrix)
       if(jtr > 400)then
       write(6,*) ' problem in gen_lattice_matrix '
       call gen_latnosym(amatrix,cellvol0)
                    endif
       call check_lat(llattice,lcate,amatrix)
       if(.not. llattice)then
!      write(6,*) 'not a lattice'
       goto  111
                         else
       goto  222
                         endif
       enddo
 222   continue
       call latmat(s6,amatrix,0)
       end subroutine gen_lattice_matrix
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine gen_latnosym(wmat,volume)
       implicit none
       real*8 volume,wmat(3,3)
       real*8 randomlat(6),rmat(3,3),tmq,tmr,slat(6),pi
       integer i
       real ranmar
 
        wmat=0.0d0
        do i=1,6
        randomlat(i)=ranmar()
        enddo
        pi=4.0d0*atan(1.0d0)
        do i=4,6
        randomlat(i)=randomlat(i)*pi/2.d0
        enddo
        call latmat(randomlat,rmat,1)
        tmr=(rmat(1,2)*rmat(2,3)-rmat(1,3)*rmat(2,2))*rmat(3,1)  &
           +(rmat(1,3)*rmat(2,1)-rmat(1,1)*rmat(2,3))*rmat(3,2)  &
           +(rmat(1,1)*rmat(2,2)-rmat(1,2)*rmat(2,1))*rmat(3,3)
        tmq=volume/tmr ; tmq=tmq**(1.0d0/3.0d0)
        slat(1)=randomlat(1)*tmq
        slat(2)=randomlat(2)*tmq
        slat(3)=randomlat(3)*tmq
        slat(4)=randomlat(4) ; slat(5)=randomlat(5) ; slat(6)=randomlat(6)
        call latmat(slat,wmat,1)
        end subroutine gen_latnosym
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine check_lat_cyc(t6,cellvol0,llattice)
       USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY : IEEE_IS_FINITE,IEEE_IS_NAN
       implicit none
       logical llattice
       real*8 t6(6),cellvol0
       real*8 cmatrix(3,3),bmatrix(3,3),s6(6)
       integer i,j

       llattice=.false.
       call latmatvol(t6,cmatrix,cellvol0)
!---{
!
       do i=1,3
       do j=1,3
       if(ieee_is_nan(cmatrix(i,j)))then
       llattice=.false.
                                    return
                                    endif
       if(.not. ieee_is_finite(cmatrix(i,j)))then
       llattice=.false.
                                             return
                                             endif
       enddo
       enddo
!
!---}
       call latmat(t6,bmatrix,1)
       call latmat(s6,bmatrix,0)
!      write(6,'(6f20.10)') ((s6(i)-t6(i)),i=1,6)
       if(sum(abs(s6-t6)) < 1.d-10) llattice=.true.
!      write(6,'(6f20.10,2x,l1)') ((s6(i)-t6(i)),i=1,6), llattice
       end subroutine check_lat_cyc
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine check_lat(lflag,l2d,cmatrix)
       USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY : IEEE_IS_FINITE,IEEE_IS_NAN
       implicit none
       logical lflag,l2d
       real*8 cmatrix(3,3)
       real*8 altm(3,3),ra,rb,rc,alpha,beta,gama,cosinea,cosineb,cosinec,pi,tmp
       real*8 uec(3),vec(3),wec(3),uu,vv,ww
       real*8 cosine1,cosine2,cosine3
       integer i,j
  
       altm=cmatrix
       pi=4.0d0*atan(1.0d0)
!---{
!
       do i=1,3
       do j=1,3
       if(ieee_is_nan(altm(i,j)))then
       lflag=.false.     
                                 return
                                 endif
       if(.not. ieee_is_finite(altm(i,j)))then
       lflag=.false.     
                                          return
                                          endif
       enddo
       enddo
!
!---}
       ra=sqrt(altm(1,1)**2+altm(1,2)**2+altm(1,3)**2)
       rb=sqrt(altm(2,1)**2+altm(2,2)**2+altm(2,3)**2)
       rc=sqrt(altm(3,1)**2+altm(3,2)**2+altm(3,3)**2)
       cosinea=(altm(2,1)*altm(3,1)+altm(2,2)*altm(3,2)+altm(2,3)*altm(3,3))/rb/rc
       cosineb=(altm(1,1)*altm(3,1)+altm(1,2)*altm(3,2)+altm(1,3)*altm(3,3))/rc/ra
       cosinec=(altm(1,1)*altm(2,1)+altm(1,2)*altm(2,2)+altm(1,3)*altm(2,3))/ra/rb  
       tmp=180.0d0/pi
       alpha=tmp*acos(cosinea) ; beta=tmp*acos(cosineb) ; gama=tmp*acos(cosinec)
       lflag=.true.
       if(.not. l2d)then
       if(ra    < 1.2d0  .or. rb    < 1.2d0 .or. rc < 1.2d0) lflag=.false.
       if(alpha < 20.0d0 .or. alpha > 160.0d0 ) lflag=.false.
       if(beta  < 20.0d0 .or. beta  > 160.0d0 ) lflag=.false.
       if(gama  < 20.0d0 .or. gama  > 160.0d0 ) lflag=.false.
       if(ra/rb > 6.0d0  .or. ra/rb < 0.3d0 ) lflag=.false.
       if(ra/rc > 6.0d0  .or. ra/rc < 0.3d0 ) lflag=.false.
       if(rb/rc > 6.0d0  .or. rb/rc < 0.3d0 ) lflag=.false.
                    else
       if(ra    < 1.2d0 .or. rb     < 1.2d0 ) lflag=.false.
       if(alpha < 20.0d0 .or. alpha > 160.0d0 ) lflag=.false.
       if(beta  < 20.0d0 .or. beta  > 160.0d0 ) lflag=.false.
       if(gama  < 20.0d0 .or. gama  > 160.0d0 ) lflag=.false.
       if(ra/rb > 6.0d0 .or. ra/rb  < 0.3d0 ) lflag=.false.
                    endif
!      October 5, 2017
       uec(:)=altm(1,:)+altm(2,:)
       vec(:)=altm(2,:)+altm(3,:)
       wec(:)=altm(3,:)+altm(1,:)
       uu=sqrt(uec(1)**2+uec(2)**2+uec(3)**2)
       vv=sqrt(vec(1)**2+vec(2)**2+vec(3)**2)
       ww=sqrt(wec(1)**2+wec(2)**2+wec(3)**2)
       cosine1=(uec(1)*altm(3,1)+uec(2)*altm(3,2)+uec(3)*altm(3,3))/uu/rc
       cosine2=(vec(1)*altm(1,1)+vec(2)*altm(1,2)+vec(3)*altm(1,3))/vv/ra
       cosine3=(wec(1)*altm(2,1)+wec(2)*altm(2,2)+wec(3)*altm(2,3))/ww/rb
       tmp=180.0d0/pi
       cosine1=tmp*acos(cosine1) 
       cosine2=tmp*acos(cosine2) 
       cosine3=tmp*acos(cosine3)
       if(l2d)then
       if(cosine1 < 20.d0 .or. cosine1 > 160.d0) lflag=.false.
              else
       if(cosine1 < 20.d0 .or. cosine1 > 160.d0) lflag=.false.
       if(cosine2 < 20.d0 .or. cosine2 > 160.d0) lflag=.false.
       if(cosine3 < 20.d0 .or. cosine3 > 160.d0) lflag=.false.
              endif
       end subroutine check_lat
!       
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine master_slave(nwork,ndir,kcmd)
       implicit none
       integer nwork,ndir,kcmd
       integer iw,mm
       integer, allocatable :: iseq(:)
       logical, allocatable :: loccupied(:)

       if(nwork <=0) return
!
       if(kcmd == 2)then
       allocate(iseq(nwork))
       call perturbation_seq(npert,nmate,nwork,iseq)
                    else
       if(nwork > 0)then
       allocate(iseq(nwork))
       iseq=1
                    endif
                    endif
!
       iw=ndir ;  call gen_directories(iw)
       allocate(loccupied(ndir))
       loccupied=.false.
       mm=0
       do iw=1,min(nwork,ndir)
       mm=mm+1
       call send_exe(mm,ndir,loccupied,kcmd,nwork,iseq)
       enddo
!
       do iw=1,nwork
       call receive(ndir,loccupied,kcmd)
       if(drate < 0.0d0 .and. iw > npop) dcut=dcut*abs(drate)
       if(dcut < davg*1.d-3) dcut=davg*1.d-3
!-----[
       call csa_bank_dump(0)
!-----]
       if(mm < nwork)then
       mm=mm+1
       call send_exe(mm,ndir,loccupied,kcmd,nwork,iseq)
                     endif
       enddo
       iw=-ndir ; call gen_directories(iw)
       if(allocated(iseq)) deallocate(iseq)
       if(allocated(loccupied)) deallocate(loccupied)
       end subroutine master_slave
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine receive(ndir,loccupied,kcmd)
       USE csa_application, ONLY : nlopt,lpbc,iobj
       implicit none
       integer ndir,kcmd
       logical loccupied(ndir)
       integer jd,man,ifile
       integer lcmd
       real*8 egp1,egp2,test,efermi,tstm1,tstm2,gapsize
       character*280 file_names(20)
       character*280 pname1
       character*280 cmd
       logical lfault1,lfault2,lfault3,lfault_stdout
       logical lexist,lexist1,lexist2,lexist3,lexist20
       real ranmar

       lcmd=iobj
       do 
       lfault1=.false.
       lfault2=.false.
       lfault3=.false.
       jd=dble(ranmar())*ndir+1
!      if(loccupied(jd))then
       call iofilearray(jd,file_names)
       inquire(file=trim(file_names(6)),exist=lexist)
       inquire(file=trim(file_names(7)),exist=lexist1)
       inquire(file=trim(file_names(10)),exist=lexist3)
       if(ranmar() < 0.01) &
       call genstopcar(file_names(2),file_names(6),file_names(18),file_names(17))
       lexist2=.false. ; call jobstatus(file_names(19),ifile) ; if(ifile == 1) lexist2=.true.
       if(lexist2)then
       if(lexist )then
       if(lexist1)then
       lfault_stdout=.false.
       inquire(file=trim(file_names(20)),exist=lexist20)
       if(.not. lexist20) lfault_stdout=.true.
       if(lexist20) call read_stdout_log(file_names(20),lfault_stdout)
       if(lfault_stdout)then
!      write(6,'(a82)') 'there is a fault ; kill the job-->stdout.log (Ctrl) -> touch OUTCAR, CONTCAR, STOP'
       write(6,*) 'there is a falut, we discard this conformation'
!
       energy0=1.d21 ; egp1=0.d0 ; egp2=0.d0
       tstm1=1.d6 ; tstm2=1.d6 ; gapsize=-1.d21
       man=npop ; pname1=trim(file_names(3))//'_' 
       pname1=trim(pname1) ; call read_poscar_bac(man,pname1)
       test=1.d21
       goto 444
                        endif
       lfault1=.false.
       call read_outcar(man,file_names(3),file_names(6),efermi,lfault1)
       lfault2=.false.
       if(.not. lfault1) call read_contcar(file_names(7),lfault2)
       lfault3=.false.
       egp1=0.d0 ; egp2=0.d0
       tstm1=1.d6 ; tstm2=1.d6 ; gapsize=-1.d21
       if(lcmd == 1) call egp_test(file_names(6),file_names(12),egp1,egp2,lfault3)
       if(lcmd == 2) call eds_test(file_names(6),file_names(12),test,lfault3)
       if(lcmd == 3) call eds_test1(file_names(6),file_names(12),test,lfault3)
       if(lcmd == 4) call eds_test2(file_names(6),file_names(12),test,lfault3)
       if(lcmd == 5) call emass_test(file_names(6),file_names(12),tstm1,tstm2,gapsize,lfault3)
       if(lcmd == 6) call get_spbd(lpbc,ndeg,qosi0,energy0,lfault3)
       if(lcmd == 7) call get_intensity(file_names(13),energy0,lfault3)
  444  continue
!      enthalpy minimization 
       if(lcmd == 0) call object(egp1,egp2,energy0,lcmd,lfault_stdout)
!      direct band gap optimization 
       if(lcmd == 1) call object(egp1,egp2,energy0,lcmd,lfault_stdout)
!      electronic DOS at Fermi level maximization 
       if(lcmd == 2) energy0=test
       if(lcmd == 2) call object(egp1,egp2,energy0,lcmd,lfault_stdout)
!      electronic DOS slope at Fermi level maximization 
       if(lcmd == 3) energy0=test
       if(lcmd == 3) call object(egp1,egp2,energy0,lcmd,lfault_stdout)
!      electronic DOS derived effetive mass maximization 
       if(lcmd == 4) energy0=test
       if(lcmd == 4) call object(egp1,egp2,energy0,lcmd,lfault_stdout)
!      effetive mass minimization
       if(lcmd == 5) energy0=min(tstm1,tstm2)-gapsize
       if(lcmd == 5) call object(tstm1,tstm2,energy0,lcmd,lfault_stdout)
!      special bond length preference
       if(lcmd == 6) call object(tstm1,tstm2,energy0,lcmd,lfault_stdout)
!      xrd simulation vs experimental data
       if(lcmd == 7) call object(tstm1,tstm2,energy0,lcmd,lfault_stdout)
       if(lfault1) energy0=1.d21
       if(lfault2) energy0=1.d21
       if(lcmd == 1 .and. lfault3) energy0=1.d21
       if(lcmd == 2 .and. lfault3) energy0=1.d21
       if(lcmd == 3 .and. lfault3) energy0=1.d21
       if(lcmd == 4 .and. lfault3) energy0=1.d21
       if(lcmd == 5 .and. lfault3) energy0=1.d21
       if(lcmd == 6 .and. lfault3) energy0=1.d21
       if(lcmd == 7 .and. lfault3) energy0=1.d21
!
       nlopt=nlopt+1
       if(.not. lpbc)then
       call centering(qosi0)
                     endif
!
       loccupied(jd)=.false.
       call jobstatus0(file_names(19))
       cmd='rm -f '//trim(cwd)//trim(file_names(10)) ; cmd=trim(cmd)
       call system(cmd)
       call system('sleep 0.5')
       goto 100
                  endif
                  endif
                  endif
!                       endif
!      
       call system('sleep 0.1')
       enddo
  100  continue
       if(kcmd == 1)then
       if(man > npop) man=npop
       if(man <    1) man=npop
       posi(:,man)=qosi0(:)
       energy_sorted(man)=energy0
                    endif
       if(kcmd == -1)then
       if(man > npop1) man=npop1
       if(man <     1) man=npop1
       posi1(:,man)=qosi0(:)
       energy_sorted1(man)=energy0
                     endif
       if(kcmd == 2)then
       call csa_update_conformations()
                    endif
       call flush(6)
       call system('sleep 1')
       end subroutine receive
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine send_exe(mm,ndir,loccupied,kcmd,nwork,iseq)
       USE csa_application, ONLY : lpbc
       implicit none
       integer mm,ndir,kcmd,nwork,iseq(nwork)
       logical loccupied(ndir)
       integer jd,i,ish
       real*8 r6(6),cmatrix(3,3),a1(3),a2(3),a3(3),ddg,pi
       character*280 file_names(20),tmpname
       character*280 cmd
       logical lexist18

       do jd=1,ndir
       if(.not. loccupied(jd))then
       call iofilearray(jd,file_names)
       call csa_rnd_lattice_basis(mm,nwork,iseq)
       ish=ndeg-6
       if(lpbc)then
       do i=1,6
       r6(i)=qosi0(ish+i)
       enddo
       call latmat(r6,cmatrix,1)
       a1(:)=cmatrix(1,:) ; a2(:)=cmatrix(2,:) ; a3(:)=cmatrix(3,:)
               endif
       if(lpbc)then
       call direct_pbc(qosi0)
               else
       call centering(qosi0)
               endif
       inquire(file=trim(file_names(5)),exist=lexist18)
       if(lexist18)then
       open(44,file=trim(file_names(5)),form='formatted')
       close(44,status='delete')
                   endif
       inquire(file=trim(file_names(18)),exist=lexist18)
       if(lexist18)then
       open(44,file=trim(file_names(18)),form='formatted')
       close(44,status='delete')
                   endif
       call write_poscar(mm,file_names(3))
       pi=4.0d0*atan(1.0d0)
       ddg=(2.0d0*pi)*0.12d0 ; tmpname=trim(file_names(5))//'_012'
       call write_kpoints(ddg,a1,a2,a3,tmpname)
       ddg=(2.0d0*pi)*0.06d0 ; tmpname=trim(file_names(5))//'_006'
       call write_kpoints(ddg,a1,a2,a3,tmpname)
       ddg=(2.0d0*pi)*0.03d0 ; tmpname=trim(file_names(5))//'_003'
       call write_kpoints(ddg,a1,a2,a3,tmpname)
       ddg=(2.0d0*pi)*0.02d0 ; tmpname=trim(file_names(5))//'_002'
       call write_kpoints(ddg,a1,a2,a3,tmpname)
       ddg=(2.0d0*pi)*0.00d0 ; tmpname=trim(file_names(5))//'_000'
       call write_kpoints(ddg,a1,a2,a3,tmpname)
       call system('sleep 0.1')
       cmd='cp '//trim(file_names(3))//' '//trim(file_names(3))//'_' 
       cmd=trim(cmd) ; call system(cmd)
!
!      cmd='cd ./'//trim(file_names(1))//' ; '//'qsub ./CSA_SOLDIER.pbs'
       cmd='cd ./'//trim(file_names(1))//' ; '//'sbatch ./CSA_SOLDIER.pbs'
       cmd=trim(cmd) ; call system(cmd)
       loccupied(jd)=.true.
       call system('sleep 0.1')
       exit
                              endif
       enddo
       end subroutine send_exe
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine read_stdout_log(stdname,lfault_stdout)
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

       lfault=.false.
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
!                endif
       lfault_stdout=lfault
       call read_stdout_log1(stdname,lfault_stdout)
       end subroutine read_stdout_log
!
!234567890
!      http://error.wiki/VASP
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine read_stdout_log1(stdname,lfault_stdout)
       USE strings, ONLY : parse,value
       implicit none
       character*280 stdname
       logical lfault_stdout
       integer nargs
       character*200 str1
       character*200 args(40)
       character*20 delims
       logical lfault

       lfault=.false.
       open(18,file=trim(stdname),form='formatted')
       do
       read(18,'(a200)',err=911,end=999) str1
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
       if(args(1) == 'internal' .and. args(2) == 'ERROR' .and. args(3) == 'SETYLM_AUG:') goto 911
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
!      October 1, 2017
       if(args(2) == 'BAD' .and. args(3) == 'TERMINATION' .and. args(4) == 'OF') goto 911
                    endif
       enddo
  911  continue
       lfault=.true.
  999  continue
       close(18)
       if(lfault)then
       write(6,*) 'there is a fault sign from stdout.log'
       call system('sleep 3')
                 endif
       lfault_stdout=lfault
       end subroutine read_stdout_log1
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine read_poscar_bac(man1,pname1)
       implicit none
       integer man1
       character*280 pname1
       logical lexist
       integer ios
  
       man1=npop
       inquire(file=trim(pname1),exist=lexist)
       if(.not. lexist)then
       write(6,*) 'somehow POSCAR_ is not present'
                       endif
       open(71,file=trim(pname1),form='formatted')
       read(71,*,iostat=ios) man1
       if(ios /=0) man1=npop
       close(71)
       end subroutine read_poscar_bac
!
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine read_outcar(man0,pname,otname,efermi,lfault1)
!      USE csa_application, ONLY : au2ang,au2ev,extpress
       USE csa_application, ONLY : nspecies,nelements,symbl
       implicit none
       integer man0
       logical lfault1
       character*280 pname,otname
       real*8 efermi,ccmat(3,3),tmp
       real*8 energy00_read,cellvol0,vtest,cmat(3,3),scale1,pi,s6(6)
       logical lfault,lexist
       integer ios
       integer i
       real*8, external :: atomicmass
  
       lfault=.false.
       man0=npop
       energy00_read=1.d19
       inquire(file=trim(pname),exist=lexist)
       if(.not. lexist)then
       write(6,*) 'somehow POSCAR is not present'
       cmat=0.d0 ; cmat(1,1)=1.d0 ; cmat(2,2)=1.d0 ; cmat(3,3)=1.d0
       lfault=.true.
       goto 888
                       else
       open(71,file=trim(pname),form='formatted')
       read(71,*,iostat=ios) man0 
       if(ios /=0) goto 911
       read(71,*,iostat=ios) scale1
       if(ios /=0) goto 911
       read(71,*,iostat=ios) cmat(1,1),cmat(1,2),cmat(1,3)
       if(ios /=0) goto 911
       read(71,*,iostat=ios) cmat(2,1),cmat(2,2),cmat(2,3)
       if(ios /=0) goto 911
       read(71,*,iostat=ios) cmat(3,1),cmat(3,2),cmat(3,3)
       if(ios /=0) goto 911
       goto 999
  911  continue
       lfault=.true.
  999  continue
       close(71)
       if(lfault) goto 888
!
       if(scale1 > 0.d0)then
       cmat=cmat*scale1
                        endif
       if(scale1 < 0.d0)then
       vtest=(cmat(1,2)*cmat(2,3)-cmat(1,3)*cmat(2,2))*cmat(3,1) &
            +(cmat(1,3)*cmat(2,1)-cmat(1,1)*cmat(2,3))*cmat(3,2) &
            +(cmat(1,1)*cmat(2,2)-cmat(1,2)*cmat(2,1))*cmat(3,3)
       vtest=(abs(scale1)/abs(vtest))**(1.d0/3.d0)
       cmat=cmat*vtest
                        endif
!
                       endif
       lfault=.false.
       call get_etot_enth1(otname,energy00_read,efermi,ccmat,lfault)
       if(.not. lfault) cmat=ccmat
  888  continue
       if(.not. lfault)then
       cellvol0=(cmat(1,2)*cmat(2,3)-cmat(1,3)*cmat(2,2))*cmat(3,1) &
               +(cmat(1,3)*cmat(2,1)-cmat(1,1)*cmat(2,3))*cmat(3,2) &
               +(cmat(1,1)*cmat(2,2)-cmat(1,2)*cmat(2,1))*cmat(3,3)
       cellvol0=abs(cellvol0)
!
       tmp=0.d0
       do i=1,nspecies
       tmp=tmp+atomicmass(symbl(i))*nelements(i)
       enddo
       tmp=tmp/cellvol0*(1.660539d-24)/(1.d-24)
!
       pi=4.d0*atan(1.d0)
       call latmat(s6,cmat,0)
       write(6,'(i6,1x,f10.3,1x,f16.7,1x,a6,1x,3f11.4,1x,4f10.4)') man0,cellvol0,energy00_read, 'outcar', &
       s6(1),s6(2),s6(3),s6(4)*180.d0/pi,s6(5)*180.d0/pi,s6(6)*180.d0/pi,tmp
                       endif
       if(lfault)then
       energy00_read=1.d19
       write(6,*) 'there is a fault with OUTCAR'
       write(6,*) 'trial conformation is discarded with high energy'
                 endif
       energy0=energy00_read
       lfault1=lfault
       end subroutine read_outcar
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine get_etot_enth1(otname,energy00_read,efermi,ccmat,lfault)
       USE strings, ONLY : parse,value
       implicit none
       character*280 otname
       logical lfault
       real*8 energy00_read,efermi,ccmat(3,3)
       character*200 str1
       integer ios,nargs
       character*200 args(40)
       character*20 delims
       real*8 enthalpy,etot
       logical lenth
 
       efermi=1.d8
       etot=1.d8
       enthalpy=1.d8
       lfault=.false.
       lenth=.false.
       open(81,file=trim(otname),form='formatted')
       do 
       read(81,'(a200)',err=911,end=999) str1
       delims=' '
       call parse(str1,delims,args,nargs)
       if(nargs == 7)then
       if(args(1) == 'energy'  )then
       if(args(2) == 'without' )then
       if(args(3) == 'entropy=')then
!      call value(args(4),etot,ios)
       call value(args(7),etot,ios)
       if(ios /= 0) etot=1.d8
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
       call value(args(5),enthalpy,ios)
       if(ios /= 0) enthalpy=1.d8
!      print*, str1
!      print*, enthalpy
                                endif
                                endif
                                endif
                     endif
       if(nargs >= 7)then
       if(args(1) == 'E-fermi')then
       if(args(2) == ':'      )then
       call value(args(3),efermi,ios)
       if(ios /= 0) efermi=1.d8
!      print*, str1
!      print*, efermi
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
       read(81,*,err=911,end=999) ccmat(1,1),ccmat(1,2),ccmat(1,3)
       read(81,*,err=911,end=999) ccmat(2,1),ccmat(2,2),ccmat(2,3)
       read(81,*,err=911,end=999) ccmat(3,1),ccmat(3,2),ccmat(3,3)
                                  endif
                                  endif
                                  endif
                                  endif
                                  endif
                                  endif
                     endif
!      if(nargs == 8)then
!      if(args(1) == 'General')then
!      if(args(2) == 'timing' )then
!      if(args(3) == 'and'    )then
!      if(args(4) == 'accounting')then
!      print*, str1
!      print*, 'exit',enthalpy
!      goto 999
!                                 endif
!                              endif
!                              endif
!                              endif
!                    endif
       enddo
  911  continue
       lfault=.true.
  999  continue
       close(81)
       energy00_read=etot
       if(lenth) energy00_read=enthalpy
       if(efermi >= 1.d8 .or. energy00_read >= 1.d8) lfault=.true.
       end subroutine get_etot_enth1
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine read_contcar(cname,lfault2)
       USE csa_application, ONLY : nspecies,nelements
       USE csa_application, ONLY : lpbc
       implicit none
       logical lfault2
       character*280 cname
       integer i,j,na,ish
       real*8 scale1,aa1(3),aa2(3),aa3(3),x,y,z,d1,d2,d3,cmatrix(3,3),r6(6),vtest
       character*200 string0
       character*9 ch9_1,ch9_2
       logical lexist
  
       lfault2=.false.
       inquire(file=trim(cname),exist=lexist)
       if(.not. lexist)then
       write(6,*) 'somehow CONTCAR is not present'
       lfault2=.true.
                       return
                       endif
       open(81,file=trim(cname),form='formatted')
       read(81,*,err=911,end=999) string0
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
       call latmat(r6,cmatrix,0)
       ish=ndeg-6
       do i=1,6
       qosi0(ish+i)=r6(i)
       enddo
!
       read(81,*,err=911,end=999) 
       read(81,*,err=911,end=999) (nelements(j),j=1,nspecies) 
       na=sum(nelements)
       read(81,*,err=911,end=999) ch9_1
       if(ch9_1(1:1) == 'S') ch9_1='Selective'
       if(ch9_1(1:1) == 's') ch9_1='Selective'
       if(ch9_1(1:1) == 'C') ch9_1='Cartesian'
       if(ch9_1(1:1) == 'c') ch9_1='Cartesian'
       if(ch9_1(1:1) == 'K') ch9_1='Cartesian'
       if(ch9_1(1:1) == 'k') ch9_1='Cartesian'
       if(ch9_1(1:1) == 'D') ch9_1='Direct'
       if(ch9_1(1:1) == 'd') ch9_1='Direct'
       if(ch9_1 == 'Selective' .or. ch9_1 == 'selective')then
       read(81,*,err=911,end=999) ch9_2
       if(ch9_2(1:1) == 'C') ch9_2='Cartesian'
       if(ch9_2(1:1) == 'c') ch9_2='Cartesian'
       if(ch9_2(1:1) == 'K') ch9_2='Cartesian'
       if(ch9_2(1:1) == 'k') ch9_2='Cartesian'
       if(ch9_2(1:1) == 'D') ch9_2='Direct'
       if(ch9_2(1:1) == 'd') ch9_2='Direct'
       if(ch9_2 == 'Cartesian' .or. ch9_2 == 'cartesian')then
       do i=1,na
       read(81,*,err=911,end=999) x,y,z
       if(scale1 > 0.d0)then
       qosi0(3*(i-1)+1)=x*scale1
       qosi0(3*(i-1)+2)=y*scale1
       qosi0(3*(i-1)+3)=z*scale1
                        endif
       enddo
       if(lpbc) call tolatx(qosi0)
       goto 900
                                                         else
       do i=1,na
       read(81,*,err=911,end=999) d1,d2,d3
       qosi0(3*(i-1)+1)=d1
       qosi0(3*(i-1)+2)=d2
       qosi0(3*(i-1)+3)=d3
       enddo
       goto 900
                                                         endif
                                                         else
       if(ch9_1 == 'Cartesian' .or. ch9_1 == 'cartesian')then
       do i=1,na
       read(81,*,err=911,end=999) x,y,z
       if(scale1 > 0.d0)then
       qosi0(3*(i-1)+1)=x*scale1
       qosi0(3*(i-1)+2)=y*scale1
       qosi0(3*(i-1)+3)=z*scale1
                        endif
       enddo
       if(lpbc) call tolatx(qosi0)
       goto 900
                                                         else
       do i=1,na
       read(81,*,err=911,end=999) d1,d2,d3
       qosi0(3*(i-1)+1)=d1
       qosi0(3*(i-1)+2)=d2
       qosi0(3*(i-1)+3)=d3
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
       call latmat(r6,cmatrix,0)
       ish=ndeg-6
       do i=1,6
       qosi0(ish+i)=r6(i)
       enddo
                        endif
       end subroutine read_contcar
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine contcarseries(mmdummy0,test0,kase0)
       USE csa_application, ONLY : nspecies,nelements,symbl
       implicit none
       integer mmdummy0,kase0
       real*8 test0
       integer i,j,na,ish,jchoice
       real*8 r6(6),cmat(3,3),a1(3),a2(3),a3(3),vec(3)

       jchoice=0
       ish=ndeg-6
       do i=1,6
       r6(i)=qosi0(ish+i)
       enddo
       call latmat(r6,cmat,1)
       a1(:)=cmat(1,:) ; a2(:)=cmat(2,:) ; a3(:)=cmat(3,:)
       write(2,'(i7,1x,e24.12,1x,a13,1x,i3)') mmdummy0,test0,'contcarseries',kase0
       write(2,'(a3)') '1.0'
       write(2,'(3f23.16)') a1(1),a1(2),a1(3)
       write(2,'(3f23.16)') a2(1),a2(2),a2(3)
       write(2,'(3f23.16)') a3(1),a3(2),a3(3)
       write(2,'(20(2x,a2,1x))') (symbl(i),i=1,nspecies)
       write(2,'(20(i4,1x))') (nelements(i),i=1,nspecies)
       if(jchoice == 0) write(2,'(a6)') "Direct"
       if(jchoice /= 0) write(2,'(a9)') "Cartesian"
       na=0
       do i=1,nspecies
       do j=1,nelements(i)
       na=na+1
       enddo
       enddo
       do i=1,na*3
       call tonormal(qosi0(i))
       enddo
       na=0
       do i=1,nspecies
       do j=1,nelements(i)
       na=na+1
       if(jchoice /= 0)then
       vec(1)=qosi0(3*(na-1)+1)*a1(1)+qosi0(3*(na-1)+2)*a2(1)+qosi0(3*(na-1)+3)*a3(1)
       vec(2)=qosi0(3*(na-1)+1)*a1(2)+qosi0(3*(na-1)+2)*a2(2)+qosi0(3*(na-1)+3)*a3(2)
       vec(3)=qosi0(3*(na-1)+1)*a1(3)+qosi0(3*(na-1)+2)*a2(3)+qosi0(3*(na-1)+3)*a3(3)
       write(2,'(3f25.16)') vec(1),vec(2),vec(3)
                       else
       write(2,'(3f20.16)') qosi0(3*(na-1)+1),qosi0(3*(na-1)+2),qosi0(3*(na-1)+3)
                       endif
       enddo
       enddo
       call flush(2)
       end subroutine contcarseries
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine write_poscar(mm0,pname)
       USE csa_application, ONLY : nspecies,nelements,symbl
       implicit none
       integer mm0
       character*280 pname
       character*280 sname
       integer i,j,na,ish,jchoice
       real*8 r6(6),cmat(3,3),a1(3),a2(3),a3(3),vec(3)
       logical lexist

       jchoice=0
       ish=ndeg-6
       do i=1,6
       r6(i)=qosi0(ish+i)
       enddo
       call latmat(r6,cmat,1)
       a1(:)=cmat(1,:) ; a2(:)=cmat(2,:) ; a3(:)=cmat(3,:)
       open(71,file=trim(pname),form='formatted')
       write(71,*) mm0 
       write(71,'(a3)') '1.0'
       write(71,'(3f23.16)') a1(1),a1(2),a1(3)
       write(71,'(3f23.16)') a2(1),a2(2),a2(3)
       write(71,'(3f23.16)') a3(1),a3(2),a3(3)
       write(71,'(20(2x,a2,1x))') (symbl(i),i=1,nspecies)
       write(71,'(20(i4,1x))') (nelements(i),i=1,nspecies)
       if(jchoice == 0) write(71,'(a6)') "Direct"
       if(jchoice /= 0) write(71,'(a9)') "Cartesian"
       na=0
       do i=1,nspecies
       do j=1,nelements(i)
       na=na+1
       enddo
       enddo
       do i=1,na*3
       call tonormal(qosi0(i))
       enddo
       na=0
       do i=1,nspecies
       do j=1,nelements(i)
       na=na+1
       if(jchoice /= 0)then
       vec(1)=qosi0(3*(na-1)+1)*a1(1)+qosi0(3*(na-1)+2)*a2(1)+qosi0(3*(na-1)+3)*a3(1)
       vec(2)=qosi0(3*(na-1)+1)*a1(2)+qosi0(3*(na-1)+2)*a2(2)+qosi0(3*(na-1)+3)*a3(2)
       vec(3)=qosi0(3*(na-1)+1)*a1(3)+qosi0(3*(na-1)+2)*a2(3)+qosi0(3*(na-1)+3)*a3(3)
       write(71,'(3f25.16)') vec(1),vec(2),vec(3)
                       else
       write(71,'(3f20.16)') qosi0(3*(na-1)+1),qosi0(3*(na-1)+2),qosi0(3*(na-1)+3)
                       endif
       enddo
       enddo
       close(71)
       sname=trim(pname) ; j=len_trim(sname)-6
       sname=trim(sname(1:j))//'STOPCAR'
       inquire(file=trim(sname),exist=lexist)
       if(lexist)then
       open(90,file=trim(sname),form='formatted')
       close(90,status='delete')
                 endif
       end subroutine write_poscar
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine write_kpoints(dk,a1,a2,a3,kname)
       implicit none
       character*280 kname
       real*8 dk,a1(3),a2(3),a3(3)
       real*8 b1(3),b2(3),b3(3),omega
       real*8 dum,pi,ga,gb,gc
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
       call meshijk(dk,ga,dum,ka) ; call meshijk(dk,gb,dum,kb) ; call meshijk(dk,gc,dum,kc)
       if(ka > 15) ka=15
       if(kb > 15) kb=15
       if(kc > 15) kc=15
       if(iswitch == 1)then
       ka=ka+1+ka/2
       kb=kb+1+kb/2
       kc=kc+1+kc/2
                       endif
       open(71,file=trim(kname),form='formatted')
       write(71,"(a1)") 'A'
       write(71,"(i1)")  0
       write(71,"(a1)") 'G'
       write(71,"(1x,i3,2x,i3,2x,i3)") ka,kb,kc
       write(71,"(1x,i2,2x,i2,2x,i2)") 0,0,0
!      write(71,"(1x,3f4.2)") 0.,0.,0.
       close(71)
       return
       end subroutine write_kpoints
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine meshijk(dk,ga,dum,ka)
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
       return
       end subroutine meshijk
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine centering(qqq)
       USE csa_application, ONLY : natom
       implicit none
       real*8 qqq(ndeg)
       real*8 vec(3)
       integer i
 
       vec=0.d0
       do i=1,natom
       vec(1)=vec(1)+qqq(3*(i-1)+1)
       vec(2)=vec(2)+qqq(3*(i-1)+2)
       vec(3)=vec(3)+qqq(3*(i-1)+3)
       enddo
       vec=vec/dble(natom)
       do i=1,natom
       qqq(3*(i-1)+1)=qqq(3*(i-1)+1)-vec(1)
       qqq(3*(i-1)+2)=qqq(3*(i-1)+2)-vec(2)
       qqq(3*(i-1)+3)=qqq(3*(i-1)+3)-vec(3)
       enddo
       end subroutine centering
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine gen_randrot(rx,ry,rz)
       implicit none
       real*8 rx(3,3),ry(3,3),rz(3,3)
       real*8 pi,theta,costheta,sintheta
       real ranmar

       pi=4.d0*atan(1.d0)
       theta=ranmar()*2.d0*pi ; costheta=cos(theta) ; sintheta=sin(theta)
       rx(1,1)=1.d0
       rx(1,2)=0.d0
       rx(1,3)=0.d0
       rx(2,1)=0.d0
       rx(2,2)=costheta
       rx(2,3)=-sintheta
       rx(3,1)=0.d0
       rx(3,2)=-sintheta
       rx(3,3)=costheta
       theta=ranmar()*2.d0*pi ; costheta=cos(theta) ; sintheta=sin(theta)
       ry(1,1)=costheta
       ry(1,2)=0.d0
       ry(1,3)=sintheta
       ry(2,1)=0.d0
       ry(2,2)=1.d0
       ry(2,3)=0.d0
       ry(3,1)=-sintheta
       ry(3,2)=0.d0
       ry(3,3)=costheta
       theta=ranmar()*2.d0*pi ; costheta=cos(theta) ; sintheta=sin(theta)
       rz(1,1)=costheta
       rz(1,2)=-sintheta
       rz(1,3)=0.d0
       rz(2,1)=sintheta
       rz(2,2)=costheta
       rz(2,3)=0.d0
       rz(3,1)=0.d0
       rz(3,2)=0.d0
       rz(3,3)=1.d0
       end subroutine gen_randrot
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine direct_pbc(qqq)
       USE csa_application, ONLY : natom
       implicit none
       real*8 qqq(ndeg)
       real*8 vec(3)
       integer i

       do i=1,natom
       vec(1)=qqq(3*(i-1)+1)
       vec(2)=qqq(3*(i-1)+2)
       vec(3)=qqq(3*(i-1)+3)
       vec(1)=vec(1)-anint(vec(1)) 
       vec(2)=vec(2)-anint(vec(2)) 
       vec(3)=vec(3)-anint(vec(3)) 
       if(vec(1) < 0.d0) vec(1)=vec(1)+1.d0
       if(vec(2) < 0.d0) vec(2)=vec(2)+1.d0
       if(vec(3) < 0.d0) vec(3)=vec(3)+1.d0
       qqq(3*(i-1)+1)=vec(1)
       qqq(3*(i-1)+2)=vec(2)
       qqq(3*(i-1)+3)=vec(3)
       enddo
       end subroutine direct_pbc
!
       end module csa
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine gen_directories(ndir)
       implicit none
       integer ndir
       character*280 string
       character*280 cmd
       integer isize,i
       logical lexist

       isize=4
       if(ndir > 0)then
       do i=1,ndir
       call xnumeral(i,string,isize) ; string=trim(string)
       cmd='mkdir '//trim(string)                     ; cmd=trim(cmd) ; call system(cmd)
       cmd='cp ./CSA_SOLDIER.pbs '//trim(string)//'/' ; cmd=trim(cmd) ; call system(cmd)
       cmd='cp INCAR_rlx '//trim(string)//'/'         ; cmd=trim(cmd) ; call system(cmd)
       cmd='cp INCAR_rlxall '//trim(string)//'/'      ; cmd=trim(cmd) ; call system(cmd)
       cmd='cp INCAR_bs '//trim(string)//'/'          ; cmd=trim(cmd) ; call system(cmd)
       cmd='cp POTCAR  '//trim(string)//'/'           ; cmd=trim(cmd) ; call system(cmd)
       cmd=trim(string)//'/STOPCAR'
       inquire(file=trim(cmd),exist=lexist)
       if(lexist)then
       cmd='rm -rf '//trim(string)//'/STOPCAR'        ; cmd=trim(cmd) ; call system(cmd)
                 endif
       cmd=trim(string)//'/STOP'
       inquire(file=trim(cmd),exist=lexist)
       if(lexist)then
       cmd='rm -rf '//trim(string)//'/STOP'           ; cmd=trim(cmd) ; call system(cmd)
                 endif
       cmd=trim(string)//'/STATUS'
       inquire(file=trim(cmd),exist=lexist)
       if(lexist)then
       cmd='rm -rf '//trim(string)//'/STATUS'         ; cmd=trim(cmd) ; call system(cmd)
                 endif
       enddo
                   endif
       if(ndir < 0)then
       ndir=iabs(ndir)
       do i=1,ndir
       call xnumeral(i,string,isize) ; string=trim(string)
       cmd='rm -rf '//trim(string)//'/' ; cmd=trim(cmd)
       call system(cmd)
       call system('sleep 1')
       enddo
                   endif
       call system('sleep 1')
       return
       end
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine iofilearray(id,file_names)
       implicit none
       integer id
       character*280 file_names(20)
       integer isize
       character*280 string

       isize=4
       call xnumeral(id,string,isize)
       string=trim(string)
       file_names(1)=trim(string)//'/'
       file_names(2)=trim(string)//'/INCAR'
       file_names(3)=trim(string)//'/POSCAR'
       file_names(4)=trim(string)//'/POTCAR'
       file_names(5)=trim(string)//'/KPOINTS'
       file_names(6)=trim(string)//'/OUTCAR'
       file_names(7)=trim(string)//'/CONTCAR'
       file_names(10)=trim(string)//'/STOP'
       file_names(17)=trim(string)//'/OSZICAR'
       file_names(18)=trim(string)//'/STOPCAR'
       file_names(19)=trim(string)//'/STATUS'
       file_names(11)=trim(string)//'/DOSCAR'
       file_names(12)=trim(string)//'/EIGENVAL'
       file_names(13)=trim(string)//'/xrd.txt'
       file_names(20)=trim(string)//'/stdout.log'
       return
       end
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine jobstatus(fname,i)
       implicit none
       character*280 fname
       integer i
       logical lexist
       character*20 ch

       i=0
       inquire(file=trim(fname),exist=lexist)
       if(lexist)then
       open(77,file=trim(fname),form='formatted')
       do 
       read(77,*,end=999) ch 
       enddo 
  999  continue
       close(77)
       if(trim(ch) == 'DONE' .or. trim(ch) == 'done') i=1
       if(trim(ch) == 'Done' .or. trim(ch) == 'DOne') i=1
       if(trim(ch) == 'DONe'                        ) i=1
                 endif
       return
       end
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine jobstatus0(fname)
       implicit none
       character*280 fname
       character*280 cmd

!      open(77,file=trim(fname),form='formatted')
!      write(77,*) 'ING' 
!      close(77)
       cmd='echo "ING" >> '//trim(fname) ; cmd=trim(cmd)
       call system(cmd)
       call system('sleep 0.1')
       return
       end
!234567890
!      Written by In-Ho Lee, KRISS, November 7, 2018.
       subroutine genstopcar(incarname,outcarname,stopcarname,oszicarname)
       USE strings, ONLY : parse,value,lowercase
       implicit none
       character*280 incarname,outcarname,stopcarname,oszicarname
       integer i,j,kount,nsw,nequal
       real*8 test,arr(10000)
       logical lfault9,lexist9
       integer ios,nargs
       character*280 str1
       character*280 args(40)
       character*20 delims

       inquire(file=trim(incarname),exist=lexist9)
       if(.not. lexist9) goto 333
       inquire(file=trim(outcarname),exist=lexist9)
       if(.not. lexist9) goto 333
       inquire(file=trim(oszicarname),exist=lexist9)
       if(.not. lexist9) goto 333
       lfault9=.false.
       nsw=0
       open(38,file=trim(incarname),form='formatted')
       do
       read(38,'(a280)',err=311,end=499) str1
       delims=' '
       call parse(str1,delims,args,nargs)
       if(nargs >= 3)then
       if(lowercase(args(1)) == 'nsw')then
       call value(args(3),nsw,ios)
!      write(6,*) nsw,' nsw'
                                      endif
                     endif
       enddo
  311  continue
       lfault9=.true.
  499  continue
       close(38)
       if(nsw == 0)then
                   goto 333
                   endif
       do i=1,10000
       arr(i)=dble(i)
       enddo
       nequal=10
       kount=0
       i=0
       if(i == 1)then
       lfault9=.false.
       open(39,file=trim(outcarname),form='formatted')
       do
       read(39,'(a280)',err=111,end=299) str1
       delims=' '
       call parse(str1,delims,args,nargs)
       if(nargs == 7)then
       if(args(3) == 'entropy=')then
       if(args(1) == 'energy')then
       call value(args(7),test,ios)
       kount=kount+1 ; arr(kount)=test
       if(kount == 10000) kount=0
                              endif
                                endif
                     endif
       enddo
  111  continue
       lfault9=.true.
  299  continue
       close(39)
                 endif
       lfault9=.false.
       open(39,file=trim(oszicarname),form='formatted')
       do
       read(39,'(a280)',err=511,end=599) str1
       delims=' '
       call parse(str1,delims,args,nargs)
       if(nargs >= 7)then
       if(args(2) == 'F=')then
       if(args(4) == 'E0=')then
       call value(args(5),test,ios)
       kount=kount+1 ; arr(kount)=test
       if(kount == 10000) kount=0
                           endif
                          endif
                     endif
       enddo
  511  continue
       lfault9=.true.
  599  continue
       close(39)
!      write(6,*) kount,' steps'
       if(kount > 0)then
       j=0 ; test=arr(kount)
       if(kount >= nequal)then
       do i=kount,kount-nequal+1,-1
!      write(6,*) test,arr(i)
       if(abs(arr(i)-test) < 1.d-7)then
       j=j+1
                                   endif
       enddo
                          endif
       if(j >= nequal)then
!      write(6,*) 'condition for stopcar generation'
       i=0
       if(iabs(nsw-kount) > 20) i=1
       if(i == 1)then
       open(44,file=trim(stopcarname),form='formatted')
       write(44,*) 'LSTOP = .TRUE.'
       close(44)
                 endif
       goto 333
                      endif
                    endif 
  333  continue
       end
!234567890
!      Written by In-Ho Lee, KRISS, September 18, 2007.
       subroutine gauss(sigma,xl0,xl)
       IMPLICIT NONE
       real*8 sigma,xl0,xl
       real*8 r,v1,v2
       real ranmar

       r=2.0d0
       do while (r >= 1.d0)
       v1=2.0*ranmar()-1.0
       v2=2.0*ranmar()-1.0
       r=v1**2+v2**2
       enddo
       xl=v1*sqrt(-2.d0*log(r)/r)
       xl=xl0+sigma*xl
       end
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine perturbation_seq(npert1,nmate1,nn1,iseq1)
       implicit none
       integer npert1,nmate1,nn1,iseq1(nn1)
       integer, allocatable :: iwrk(:),jwrk(:)
       real*8, allocatable :: wrk(:)
       integer n,i
       real ranmar

       n=npert1+nmate1
       if(n /= nn1)then
       write(6,*) 'n /= nn1',n,nn1
                   stop
                   endif
       allocate(iwrk(n),jwrk(n)) ; allocate(wrk(n))
       do i=1,n
       wrk(i)=ranmar()
       enddo
       do i=1,npert1
       jwrk(i)=1
       enddo
       do i=1,nmate1
       jwrk(i+npert1)=2
       enddo
       call sortnr(n,wrk,iwrk)
       do i=1,n
       iseq1(i)=jwrk(iwrk(i))
       enddo
       deallocate(iwrk,jwrk) ; deallocate(wrk)
       end
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine onedffvltn(ndim,npop,x,y)
       implicit none
       integer ndim,npop
       real*8 x(ndim,npop),y(ndim)
       real*8, allocatable :: xnew(:,:)
       integer j

       allocate(xnew(ndim,npop))
       call dffvltn(ndim,npop,x,xnew)
       do j=1,ndim
       y(j)=xnew(j,1)
       enddo
       deallocate(xnew)
       end
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine onedffvltn1(ndeg,npop,x,y)
       implicit none
       integer ndeg,npop
       real*8 x(ndeg,npop),y(ndeg)
       real*8, allocatable :: xnew(:,:)
       real*8, allocatable :: x1(:,:),x1new(:,:),x2(:,:),x2new(:,:)
       integer j,n1,n2,i

       allocate(xnew(ndeg,npop))
       call dffvltn(ndeg,npop,x,xnew)
       do j=1,ndeg
       y(j)=xnew(j,1)
       enddo
       deallocate(xnew)
       
       n1=6 ; n2=npop
       allocate(x1(n1,n2),x1new(n1,n2))
       do i=1,n2
       do j=1,n1
       x1(j,i)=x(ndeg-6+j,i)
       enddo
       enddo
       call dffvltn(n1,n2,x1,x1new)
       do j=1,n1
       y(ndeg-6+j)=x1new(j,1)
       enddo
       deallocate(x1,x1new)

       n1=ndeg-6 ; n2=npop
       allocate(x2(n1,n2),x2new(n1,n2))
       do i=1,n2
       do j=1,n1
       x2(j,i)=x(j,i)
       enddo
       enddo
       call dffvltn(n1,n2,x2,x2new)
       do j=1,n1
       y(j)=x2new(j,1)
       enddo
       deallocate(x2,x2new)
       end
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine dffvltn(ndim,npop,x,xnew)
       implicit none
       integer ndim,npop
       real*8 x(ndim,npop),xnew(ndim,npop)
       integer iarr(4),jbrr,i,j
       real*8 tmt
       real*8 cr,ff
       real ranmar
!      ff :  differential weight
!      cr :  crossover probability
       xnew=x
       do i=1,npop
       if(npop >=4 )then
       iarr(1)=dble(ranmar())*npop+1
  111  continue
       iarr(2)=dble(ranmar())*npop+1
       if(iarr(2) == iarr(1)) goto 111
  112  continue
       iarr(3)=dble(ranmar())*npop+1
       if(iarr(3) == iarr(1) .or. iarr(3) == iarr(2)) goto 112
  113  continue
       iarr(4)=dble(ranmar())*npop+1
       if(iarr(4) == iarr(3) .or. iarr(4) == iarr(2) .or. iarr(4) == iarr(1)) goto 113
!      print*, i,'--'
!      print*, iarr(1),iarr(2),iarr(3),iarr(4)
       jbrr=dble(ranmar())*ndim+1 ; cr=ranmar() ; ff=2.d0*dble(ranmar())+0.d0 
       do j=1,ndim
       tmt=ranmar()
       if(jbrr == j .or. tmt <= cr)then
       xnew(j,i)=x(j,iarr(2))+ff*(x(j,iarr(3))-x(j,iarr(4)))
                                   endif
       if(jbrr /= j .and. tmt > cr)then
       xnew(j,i)=x(j,iarr(1))
                                   endif
       enddo
                    else
       tmt=-1.d0
       do j=1,ndim
       if(tmt < abs(x(j,i))) tmt=abs(x(j,i))
       enddo
       tmt=tmt/10.d0
       do j=1,ndim
       xnew(j,i)=x(j,i)+(ranmar()-0.5)*tmt
       enddo
                    endif
       enddo
       end
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine cross3(a,b,c)
       implicit none
       real*8 a(3),b(3),c(3)

       c(1)=a(2)*b(3)-a(3)*b(2)
       c(2)=a(3)*b(1)-a(1)*b(3)
       c(3)=a(1)*b(2)-a(2)*b(1)
       return
       end
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine tonormal(xxx)
       USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY : IEEE_IS_FINITE,IEEE_IS_NAN
       implicit none
       real*8 xxx
       real ranmar

!      if(isnan(xxx)) xxx=ranmar()
       if(IEEE_IS_NAN(xxx)) xxx=ranmar()
       if(.not. IEEE_IS_FINITE(xxx)) xxx=ranmar()
       end subroutine tonormal
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine csa_vasp_banner()
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
!      Written by In-Ho Lee, KRISS, September 11, 2013.
!      program csa_vasp5_2_12
       program csa_vasp5_4_1
       USE csa, ONLY : csa_initial,csa_final,csa_first_bank,csa_evolution
       implicit none
       character*8 fnnd ; character*10 fnnt

       call date_and_time(date=fnnd,time=fnnt)
       write(6,'(a10,2x,a8,2x,a10)') 'date,time ', fnnd,fnnt

       call csa_vasp_banner()
       call timestamp()
       call csa_initial()
       call csa_first_bank()
       call csa_evolution()
       call csa_final()
       call timestamp()
       end program csa_vasp5_4_1
!
