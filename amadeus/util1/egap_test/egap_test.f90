!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       module eigenvalues
       implicit none
       private
       save
       integer ne,nk,nbandi 
       real*8 ef,temper
       logical ldirect
       real*8, allocatable :: rkpt(:,:),eiv(:,:),wgt(:),cblowest(:),vbhighest(:)

       public :: ne,nk,nbandi,rkpt,eiv,wgt,ef,temper,cblowest,vbhighest,ldirect
       end module eigenvalues
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine egp_test(otname,einame,egp1,egp2,lfault)
       USE eigenvalues, ONLY : ldirect
       implicit none
       logical lfault
       character*280 otname,einame,cname1
       character*80 string0
       integer islc
       real*8 egp1,egp2
       logical lfault1,lfault2,lfault3,lfault7

       lfault=.false.
       lfault1=.false.
       lfault2=.false.
       lfault3=.false.
       egp1=-1.d0
       egp2=-1.d0
       call read_outcar1(otname,lfault1)
       if(lfault1) goto 911
       call read_eigenval(einame,lfault2)
       if(lfault2) goto 911
       call check_gap(egp1,egp2,lfault3)
       if(.not. lfault3)then
       if(ldirect)then
       islc=len_trim(otname) ; islc=islc-6 ; cname1=otname(1:islc) ; cname1=trim(cname1)//'CONTCAR'
       lfault7=.false.
       open(71,file=trim(cname1),form='formatted')
       do
       read(71,'(a80)',err=711,end=799) string0
       if(len_trim(string0) == 0) goto 799
       if(len_trim(string0) > 0) write(6,*) trim(string0)
       enddo
  711  continue
       lfault7=.true.
  799  continue
       close(71)
                  endif
                        endif
  911  continue
       if(lfault1 .or. lfault2 .or. lfault3) lfault=.true.
       if(lfault)then
       egp1=0.d0
       egp2=0.d0
                 endif
       if(egp1 <= -1.d0) egp1=0.d0
       if(egp2 <= -1.d0) egp2=0.d0
       end
!234567890
!      Written by In-Ho Lee, KRISS, July 3, 2014.
       subroutine eds_test(otname,einame,test,lfault)
       implicit none
       logical lfault
       character*280 otname,einame
       real*8 test
       logical lfault1,lfault2,lfault3
       real*8 test0,test1

       lfault=.false.
       lfault1=.false.
       lfault2=.false.
       lfault3=.false.
       test0=0.d0
       test1=0.d0
       call read_outcar1(otname,lfault1)
       if(lfault1) goto 911
       call read_eigenval(einame,lfault2)
       if(lfault2) goto 911
       call check_dos(test0,test1,lfault3)
  911  continue
       test=test0
       if(lfault1 .or. lfault2 .or. lfault3) lfault=.true.
       end
!234567890
!      Written by In-Ho Lee, KRISS, July 3, 2014.
       subroutine eds_test1(otname,einame,test,lfault)
       implicit none
       logical lfault
       character*280 otname,einame
       real*8 test
       logical lfault1,lfault2,lfault3
       real*8 test0,test1

       lfault=.false.
       lfault1=.false.
       lfault2=.false.
       lfault3=.false.
       test0=0.d0
       test1=0.d0
       call read_outcar1(otname,lfault1)
       if(lfault1) goto 911
       call read_eigenval(einame,lfault2)
       if(lfault2) goto 911
       call check_dos(test0,test1,lfault3)
  911  continue
       test=test1
       if(lfault1 .or. lfault2 .or. lfault3) lfault=.true.
       end
!234567890
!      Written by In-Ho Lee, KRISS, December 2, 2014.
       subroutine eds_test2(otname,einame,test,lfault)
       implicit none
       logical lfault
       character*280 otname,einame
       real*8 test
       logical lfault1,lfault2,lfault3
       real*8 test0,test1

       lfault=.false.
       lfault1=.false.
       lfault2=.false.
       lfault3=.false.
       test0=0.d0
       test1=0.d0
       call read_outcar1(otname,lfault1)
       if(lfault1) goto 911
       call read_eigenval(einame,lfault2)
       if(lfault2) goto 911
       call check_emass(test0,test1,lfault3)
  911  continue
       test=test1 ; if(test > test0) test=test0
       if(lfault1 .or. lfault2 .or. lfault3) lfault=.true.
       end
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine emass_test(otname,einame,tstm1,tstm2,gapsize,lfault)
       USE strings, ONLY : parse,value
       implicit none
       character*280 otname,einame
       real*8 tstm1,tstm2,gapsize
       logical lfault,lexist
       character*280 cemcout
       character*80 string0
       integer islc,kount
       logical lfault1,lfault2
       integer ios,nargs
       character*200 args(40)
       character*20 delims

       lfault=.false.
       lfault1=.false. ; lfault2=.false.
       gapsize=-1.d22
       call read_outcar1(otname,lfault1)
       if(lfault1) goto 911
       islc=len_trim(otname) ; islc=islc-6 ; cemcout=otname(1:islc) ; cemcout=trim(cemcout)//'EMCOUT'
       kount=0
       inquire(file=trim(cemcout),exist=lexist)
       if(.not. lexist)then
       lfault1=.true.
       goto 911
                       endif
       open(72,file=trim(cemcout),form='formatted')
       do
       kount=kount+1
       read(72,'(a80)',err=711,end=799) string0
       delims=' '
       call parse(string0,delims,args,nargs)
       if(nargs >= 2)then
       if(args(1) == 'serious')then
       if(args(2) == 'problem')then
       tstm1=1.d6 ; tstm2=1.d6
       goto 711
                               endif
                               endif
                     endif
       if(nargs > 0)then
       if(kount == 1)then
       call value(args(1),tstm1,ios)
       if(ios /= 0) tstm1=1.d6
                     endif
       if(kount == 2)then
       call value(args(1),tstm2,ios)
       if(ios /= 0) tstm2=1.d6
                     endif
                    endif
       if(len_trim(string0) == 0) goto 799
       if(len_trim(string0) > 0) write(6,*) trim(string0)
       if(kount == 2) goto 799
       enddo
  711  continue
       lfault2=.true.
  799  continue
       read(72,*) gapsize
       write(6,*) gapsize,' gapsize'
       close(72)
  911  continue
       if(lfault1) lfault=.true.
       if(lfault2) lfault=.true.
       end
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine read_outcar1(otname,lfault1)
       USE eigenvalues, ONLY : ef
       implicit none
       character*280 otname
       logical lfault1
       character*7 ctest7
       logical lfault
  
       lfault=.false.
       ef=1.d19
       open(81,file=trim(otname),form='formatted')
       do 
       read(81,*,err=911,end=999) ctest7
       if(ctest7 == 'E-fermi')then
       backspace(81)
       read(81,101,err=911,end=999) ef
                              endif
  101  format(10x,f9.4)
       enddo
  911  continue
       lfault=.true.
  999  continue
       close(81)
!      write(6,*) ef,' ef from OUTCAR'
!
       if(lfault)then
       ef=1.d19
       write(6,*) 'there is a falut with OUTCAR'
                 endif
       lfault1=lfault
       end
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine read_eigenval(einame,lfault1)
       USE eigenvalues, ONLY : ef,rkpt,wgt,eiv,ne,nk,nbandi,temper
       implicit none
       character*280 einame
       logical lfault1
       integer j,i,ik
       real*8 tmpx,tmpy,tmpz,sol1,sol2
       logical lfault
       real*8, external ::  fdft,zeroin

       lfault=.false.
       ne=0
       nk=0
       nbandi=0
       open(81,file=trim(einame),form='formatted')
       read(81,*,err=911,end=999)
       read(81,*,err=911,end=999)
       read(81,*,err=911,end=999)
       read(81,*,err=911,end=999)
       read(81,*,err=911,end=999)
!      read(81,'(i5,i5,i5)',err=911,end=999) ne,nk,nbandi
! version 5.4.1
       read(81,'(i7,i7,i7)',err=911,end=999) ne,nk,nbandi
       if(nk > 100000   .or. nk <= 0) goto 911
       if(nbandi > 2000 .or. nbandi <= 0) goto 911
       if(ne > 4000     .or. ne <= 0) goto 911
       allocate(eiv(nbandi,nk))
       allocate(wgt(nk))
       allocate(rkpt(3,nk))
       do ik=1,nk
!      write(6,*) ik
       read(81,*,err=911,end=999)
       read(81,*,err=911,end=999) tmpx,tmpy,tmpz,wgt(ik)
       rkpt(1,ik)=tmpx
       rkpt(2,ik)=tmpy
       rkpt(3,ik)=tmpz
       do i=1,nbandi
       read(81,*,err=911,end=999) j,eiv(i,ik)
       enddo
       enddo
       goto 999
  911  continue
       lfault=.true.
  999  continue
       close(81)
!      just skip new calculation, just use the ef value from OUTCAR
!      write(6,*) sum(wgt)
!      1.d-5 eV = 3.1577464d0 / (2.d0*13.6058d0)  K
!      au2kelvin= 3.1577464d5
       temper=1.d-4
       sol1=ef-1.d0
       sol2=ef+1.d0
!      write(6,*) sol1,fdft(sol1)
!      write(6,*) sol2,fdft(sol2)
!      tol=1.d-16 ; ef=zeroin(sol1,sol2,fdft,tol)
!      write(6,*) ef,' ef from calculaton, EIGENVAL'
       lfault1=lfault
!      write(6,*) 'in eigenval',nk
       end
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       real*8 function fdft(x)
       USE eigenvalues, ONLY : ne,wgt,eiv,temper,nbandi,nk
       implicit none
       real*8 x
       integer ik,i
       real*8 xxr

       fdft=float(ne)
       do ik=1,nk
       do i=1,nbandi
       xxr=(eiv(i,ik)-x)/temper ; if(xxr >  50.d0) xxr=50.d0 ; if(xxr < -50.d0) xxr=-50.d0
       fdft=fdft-(2.d0*wgt(ik))/(1.d0+exp(xxr))
       enddo
       enddo
       return
       end
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine check_gap(egp1,egp2,lfault1)
       USE eigenvalues, ONLY : temper,ef,rkpt,wgt,eiv,nbandi,nk,cblowest,vbhighest,ldirect
       implicit none
       logical lfault1
       real*8 egp1,egp2
       integer ik,i,i1,i2
       real*8 t1,t2,test,tst1,tst2,occtol,smallest,dgap,xxr
       logical lfault
       real*8, allocatable :: wrk7(:),wrk8(:)
       integer, allocatable :: iwrk7(:),iwrk8(:)

       lfault=.false.
       occtol=1.d-10
!      occtol=1.d-8
!      occtol=1.d-4
       write(6,*) ef,' ef'
       dgap=1.d19
       if(nk <= 0)then
       lfault=.true.
                  endif
       if(nbandi <= 0)then
       lfault=.true.
                      endif
       if(lfault) goto 119
       allocate(cblowest(nk),vbhighest(nk))
       do ik=1,nk
       i2=nbandi
       do i=nbandi,1,-1
       if(eiv(i,ik) >= ef)then
       xxr=(eiv(i,ik)-ef)/temper ; if(xxr >  50.d0) xxr=50.d0 ; if(xxr < -50.d0) xxr=-50.d0
       t2=(2.d0*wgt(ik))/(1.d0+exp(xxr))
       if(t2 < occtol) i2=i
                          endif
       enddo
       i1=1
       do i=1,nbandi
       if(eiv(i,ik) <= ef)then
       xxr=(eiv(i,ik)-ef)/temper ; if(xxr >  50.d0) xxr=50.d0 ; if(xxr < -50.d0) xxr=-50.d0
       t1=(2.d0*wgt(ik))/(1.d0+exp(xxr))
       if(t1 > occtol) i1=i
                          endif
       enddo
!
       if(i1 <1 .or. i1 > nbandi)then
       i1=1
       lfault=.true.
       write(6,*) 'something went wrong in egap.f90'
                                 endif
       if(i2 <1 .or. i2 > nbandi)then
       i2=1
       lfault=.true.
       write(6,*) 'something went wrong in egap.f90'
                                 endif
       if(lfault) goto 119
!
       cblowest(ik)=eiv(i2,ik)
       vbhighest(ik)=eiv(i1,ik)
!      t1=(2.d0        )/(1.d0+exp((cblowest(ik)-ef)/temper))
!      t2=(2.d0        )/(1.d0+exp((vbhighest(ik)-ef)/temper))
!      write(6,*) t1,t2
!      write(6,*) cblowest(ik)-vbhighest(ik)
       test=cblowest(ik)-vbhighest(ik)
       if(dgap > test)then
       dgap=test
                      endif
       enddo
       allocate(wrk7(nk),wrk8(nk)) ; allocate(iwrk7(nk),iwrk8(nk))
       wrk7=cblowest
       wrk8=vbhighest
       call sortnr(nk,wrk7,iwrk7)
       call sortnr(nk,wrk8,iwrk8)
       write(6,*) 'cb minimum'
       write(6,'(20f12.6)') (wrk7(iwrk7(i)),i=1,min(10,nk))
       write(6,*) 'vb maximum'
       write(6,'(20f12.6)') (wrk8(iwrk8(nk-i+1)),i=1,min(10,nk))
       if(nk >= 8)then
       write(6,'(20f12.6)') wrk7(iwrk7(1))-wrk8(iwrk8(nk)), wrk7(iwrk7(2))-wrk8(iwrk8(nk)), wrk7(iwrk7(3))-wrk8(iwrk8(nk)), &
       wrk7(iwrk7(4))-wrk8(iwrk8(nk)), wrk7(iwrk7(5))-wrk8(iwrk8(nk)),&
       wrk7(iwrk7(6))-wrk8(iwrk8(nk)), wrk7(iwrk7(7))-wrk8(iwrk8(nk)), wrk7(iwrk7(8))-wrk8(iwrk8(nk))
                  endif
       smallest=minval(cblowest)-maxval(vbhighest)

       tst1=wrk7(iwrk7(1))-wrk8(iwrk8(nk))
       tst2=wrk7(iwrk7(2))-wrk8(iwrk8(nk))
       if(abs(tst1-dgap) <1.d-8)then
       write(6,*) 'direct gap',tst1,'indirect gap',tst2
       egp1=tst2
       egp2=tst1
       smallest=tst2
                                endif
       if(smallest < dgap )then
       write(6,*) 'indirect gap, direct gap ', smallest,dgap
       egp1=smallest
       egp2=dgap
                           endif
       if(smallest <  dgap)then
       ldirect=.false.
       write(6,'(a17,2x,2f16.8,2x,i6,1x,a7)') 'indirect band gap',dgap,smallest,nk,'nk fine'
                           else
       write(6,'(a15,2x,2f16.8,2x,i6,1x,a11,1x,a7)') 'direct band gap',dgap,smallest,nk,'+++++++++++','nk fine'
       ldirect=.true.
                           endif
  119  continue
       if(allocated(wrk7)) deallocate(wrk7) 
       if(allocated(wrk8)) deallocate(wrk8) 
       if(allocated(iwrk7)) deallocate(iwrk7)
       if(allocated(iwrk8)) deallocate(iwrk8)
       if(allocated(cblowest)) deallocate(cblowest)
       if(allocated(vbhighest)) deallocate(vbhighest)
       if(allocated(eiv)) deallocate(eiv)
       if(allocated(wgt)) deallocate(wgt)
       if(allocated(rkpt)) deallocate(rkpt)
       if(lfault)then
       egp1=0.d0
       egp2=0.d0
                 endif
       call flush(6)
       lfault1=lfault
       end
!234567890
!      Written by In-Ho Lee, KRISS, July 3, 2014.
       subroutine check_dos(test0,test1,lfault1)
       USE eigenvalues, ONLY : ef,rkpt,wgt,eiv,nbandi,nk,cblowest,vbhighest
       implicit none
       logical lfault1
       real*8 test0,test1
       logical lfault
       logical linsulator
       real*8 tmq,tmr,pi,dee,epl,yyr,se,sf,ds,ddss(4000)
       real*8 cbm,vbm,rqpt(3,2)
       integer ic,ig,ng,ik,i

       lfault=.false.
       linsulator=.false.
       if(nk <= 0)then
       lfault=.true.
                  endif
       if(nbandi <= 0)then
       lfault=.true.
                      endif
       if(lfault) goto 119
       pi=4.d0*atan(1.d0)
       epl=0.005d0
       ng=4000
       ddss=0.d0 
       sf=maxval(eiv)+1.0d0
       se=minval(eiv)-1.0d0
       ds=(sf-se)/float(ng-1)
       dee=5.0d0*ds
       do ik=1,nk
       do i=1,nbandi
       do ig=1,ng
       yyr=se+ds*float(ig-1)
       yyr=(yyr-eiv(i,ik))/sqrt(2.d0*epl)
       if(yyr < -12.d0) yyr=-12.d0
       if(yyr >  12.d0) yyr= 12.d0
       ddss(ig)=ddss(ig)+2.d0*wgt(ik)/sqrt(2.d0*pi*epl)*exp(-yyr**2)
       enddo
       enddo
       enddo
!
       cbm=1.d19 ; vbm=-1.d19
       do ik=1,nk
       do i=1,nbandi
       if(eiv(i,ik) >= ef .and. eiv(i,ik) <= cbm)then 
       cbm=eiv(i,ik)
       rqpt(:,2)=rkpt(:,ik)
                                                 endif
       if(eiv(i,ik) <  ef .and. eiv(i,ik) >= vbm)then 
       vbm=eiv(i,ik)
       rqpt(:,1)=rkpt(:,ik)
                                                 endif
       do ig=1,ng
       yyr=se+ds*float(ig-1)
       yyr=(yyr-eiv(i,ik))/sqrt(2.d0*epl)
       if(yyr < -12.d0) yyr=-12.d0
       if(yyr >  12.d0) yyr= 12.d0
       ddss(ig)=ddss(ig)+2.d0*wgt(ik)/sqrt(2.d0*pi*epl)*exp(-yyr**2)
       enddo
       enddo
       enddo
       if(abs(cbm-vbm) > 5.0d0*ds)then
       linsulator=.true.
       goto 119
                                  endif
!
       tmr=0.d0
       tmq=0.d0
       test0=0.d0
       do ig=1,ng
       yyr=se+ds*float(ig-1)
       if(yyr >= ef-dee/2.0 .and. yyr  <= ef+dee/2.d0)then
       test0=test0+ddss(ig)
       tmq=tmq+1.d0
                                                      endif
       if(yyr <= ef) tmr=tmr+ddss(ig)
       enddo
       tmr=tmr*ds
       test0=test0/tmq
       test0=-abs(test0/tmr)
!
       ic=1
       tmq=1.d19
       do ig=1,ng
       yyr=se+ds*float(ig-1)
       if(tmq > abs(yyr-ef))then
       tmq=abs(yyr-ef)
       ic=ig
                             endif
       enddo
       if(ic+2 <= 4000 .and. ic-2 >= 1)then
       test1=(-ddss(ic+2)+8.d0*ddss(ic+1)-8.d0*ddss(ic-1)+ddss(ic-2))/(12.d0*ds)
       test1=test1/tmr
       test1=-abs(test1)
                                       else
       test1=0.d0
                                       endif
!
  119  continue
       if(allocated(eiv)) deallocate(eiv)
       if(allocated(wgt)) deallocate(wgt)
       if(allocated(rkpt)) deallocate(rkpt)
       if(allocated(cblowest)) deallocate(cblowest)
       if(allocated(vbhighest)) deallocate(vbhighest)
       if(lfault)then
       test0=0.d0
       test1=0.d0
                 endif
       if(linsulator)then
       test0=0.d0
       test1=0.d0
                     endif
!      call flush(6)
       lfault1=lfault
       end
!234567890
!      Written by In-Ho Lee, KRISS, December 2, 2014.
       subroutine check_emass(test0,test1,lfault1)
       USE eigenvalues, ONLY : ef,rkpt,wgt,eiv,nbandi,nk,cblowest,vbhighest
       implicit none
       logical lfault1
       real*8 test0,test1
       logical lfault
       logical lmetal
       real*8 tmp,tmq,tmr,xxr,pi,dee,epl,yyr,se,sf,ds,ddss(4000)
       real*8 cbm,vbm,zzz(2),rqpt(3,2)
       integer ic,iv,ig,ng,ik,i,ide

       lmetal=.false.
       lfault=.false.
       if(nk <= 0)then
       lfault=.true.
                  endif
       if(nbandi <= 0)then
       lfault=.true.
                      endif
       if(lfault) goto 119
       pi=4.d0*atan(1.d0)
       epl=0.005d0
       ng=4000
       ddss=0.d0 
       sf=maxval(eiv)+1.0d0
       se=minval(eiv)-1.0d0
       ds=(sf-se)/float(ng-1)
       cbm=1.d19 ; vbm=-1.d19
       do ik=1,nk
       do i=1,nbandi
       if(eiv(i,ik) >= ef .and. eiv(i,ik) <= cbm)then 
       cbm=eiv(i,ik)
       rqpt(:,2)=rkpt(:,ik)
                                                 endif
       if(eiv(i,ik) <  ef .and. eiv(i,ik) >= vbm)then 
       vbm=eiv(i,ik)
       rqpt(:,1)=rkpt(:,ik)
                                                 endif
       do ig=1,ng
       yyr=se+ds*float(ig-1)
       yyr=(yyr-eiv(i,ik))/sqrt(2.d0*epl)
       if(yyr < -12.d0) yyr=-12.d0
       if(yyr >  12.d0) yyr= 12.d0
       ddss(ig)=ddss(ig)+2.d0*wgt(ik)/sqrt(2.d0*pi*epl)*exp(-yyr**2)
       enddo
       enddo
       enddo
       if(abs(cbm-vbm) < 5.0d0*ds)then
       lmetal=.true.
       goto 119
                                  endif
       tmr=0.d0
       tmp=1.d19 ; ic=ng
       tmq=1.d19 ; iv=1
       do ig=1,ng
       yyr=se+ds*float(ig-1)
       if(yyr <= ef) tmr=tmr+ddss(ig)
       if(abs(yyr-cbm) <= tmp)then
       tmp=abs(yyr-cbm)
       ic=ig
                              endif
       if(abs(yyr-vbm) <= tmq)then
       tmq=abs(yyr-vbm)
       iv=ig
                              endif
       enddo
       tmr=tmr*ds
!
       dee=0.1d0
       ide=dee/ds ; if(ide <=1) ide=1       ; ig=ic+1+ide 
       if(ig <= 4000 .and. ig >=1)then
       yyr=se+ds*float(ig-1) ; xxr=((ddss(ig)/tmr)/sqrt(yyr-cbm))**(2.d0/3.d0)
                                  else
       xxr=0.d0
                                  endif
       zzz(1)=-xxr
       ide=2.0d0*dee/ds ; if(ide <=1) ide=1 ; ig=ic+1+ide 
       if(ig <= 4000 .and. ig >=1)then
       yyr=se+ds*float(ig-1) ; xxr=((ddss(ig)/tmr)/sqrt(yyr-cbm))**(2.d0/3.d0)
                                  else
       xxr=0.d0
                                  endif
       zzz(2)=-xxr
       test1=(zzz(1)+zzz(2))/2.d0
       if(abs(zzz(1)) < 1.d-8 .or. abs(zzz(2)) < 1.d-8 ) test1=0.d0
       ide=dee/ds ; if(ide <=1) ide=1       ; ig=iv-1-ide 
       if(ig <= 4000 .and. ig >=1)then
       yyr=se+ds*float(ig-1) ; xxr=((ddss(ig)/tmr)/sqrt(vbm-yyr))**(2.d0/3.d0)
                                  else
       xxr=0.d0
                                  endif
       zzz(1)=-xxr
       ide=2.0d0*dee/ds ; if(ide <=1) ide=1 ; ig=iv-1-ide 
       if(ig <= 4000 .and. ig >=1)then
       yyr=se+ds*float(ig-1) ; xxr=((ddss(ig)/tmr)/sqrt(vbm-yyr))**(2.d0/3.d0)
                                  else
       xxr=0.d0
                                  endif
       zzz(2)=-xxr
       test0=(zzz(1)+zzz(2))/2.d0
       if(abs(zzz(1)) < 1.d-8 .or. abs(zzz(2)) < 1.d-8 ) test0=0.d0
!      write(6,'(3f18.10)') vbm,ef,cbm
!      write(6,'(3f20.12)') rqpt(1,1),rqpt(2,1),rqpt(3,1)
!      write(6,'(3f20.12)') rqpt(1,2),rqpt(2,2),rqpt(3,2)
!
  119  continue
       if(allocated(eiv)) deallocate(eiv)
       if(allocated(wgt)) deallocate(wgt)
       if(allocated(rkpt)) deallocate(rkpt)
       if(allocated(cblowest)) deallocate(cblowest)
       if(allocated(vbhighest)) deallocate(vbhighest)
       if(lfault)then
       test0=0.d0
       test1=0.d0
                 endif
       if(lmetal)then
       test0=0.d0
       test1=0.d0
                 endif
!      call flush(6)
       lfault1=lfault
       end
!
       program test
       implicit none
       character*280 otname,einame
       real*8 egp1,egp2
       logical lfault

       otname='OUTCAR'
       einame='EIGENVAL'
       call egp_test(otname,einame,egp1,egp2,lfault)

       end program test
