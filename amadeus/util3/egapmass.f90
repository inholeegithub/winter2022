!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       module eigenvalues1
       implicit none
       private
       save
       integer ne,nk,nbandi 
       real*8 ef,temper
       logical ldirect
       real*8, allocatable :: rkpt(:,:),eiv(:,:),wgt(:),cblowest(:),vbhighest(:)

       public :: ne,nk,nbandi,rkpt,eiv,wgt,ef,temper,cblowest,vbhighest,ldirect
       end module eigenvalues1
!234567890
!      Written by In-Ho Lee, KRISS, April 20, 2016.
       subroutine get_emassprincipal(einame,objtfn,lfault)
       implicit none
       character*280 einame
       real*8 objtfn
       logical lfault
       character*280 fname
       real*8 effmassv(3),effmassc(3),tst1,tst2
       integer ne0,nk0,nbandi0,ibvbm0,kvbm0,ibcbm0,kcbm0
       real*8 ef0,a10(3),a20(3),a30(3),rkptv1,rkptv2,rkptv3,rkptc1,rkptc2,rkptc3,evbm0,ecbm0
       real*8 gapsize0
       logical lfault1,lfault2

       fname=trim(einame)
!      i=len_trim(fname) ; fname=trim(fname(1:i-8))//'inp_vbm' ; fname=trim(fname)
       fname='inp_vbm' ; fname=trim(fname)
       lfault1=.false.
       call get_emass(fname,effmassv,lfault1)
       tst1=minval(effmassv)
!      write(6,*) tst1,'tst1'
       write(6,'(3f22.8)') effmassv
       fname=trim(einame)
!      i=len_trim(fname) ; fname=trim(fname(1:i-8))//'inp_cbm' ; fname=trim(fname)
       fname='inp_cbm' ; fname=trim(fname)
       lfault2=.false.
       call get_emass(fname,effmassc,lfault2)
       tst2=minval(effmassc)
!      write(6,*) tst2,'tst2'
       write(6,'(3f22.8)') effmassc
!      i=len_trim(fname) ; fname=trim(fname(1:i-8))//'els_info' ; fname=trim(fname)
       fname='els_info' ; fname=trim(fname)
       open(77,file=trim(fname),form='formatted')
       read(77,*) ne0,nk0,nbandi0
       read(77,*) ef0
       read(77,*) ibvbm0,kvbm0,rkptv1,rkptv2,rkptv3,evbm0
       read(77,*) ibcbm0,kcbm0,rkptc1,rkptc2,rkptc3,ecbm0
       read(77,*) a10(1),a10(2),a10(3)
       read(77,*) a20(1),a20(2),a20(3)
       read(77,*) a30(1),a30(2),a30(3)
       close(77)
       gapsize0=ecbm0-evbm0
       write(6,*) gapsize0
       objtfn=min(tst1,tst2)-gapsize0
       if(lfault1) lfault=.true.
       if(lfault2) lfault=.true.
       if(.not. lfault)then
       write(6,*) tst1,tst2
       write(6,*) objtfn
                       else
       tst1=1.d6 ; tst2=1.d6 ; objtfn=1.d6
       write(6,*) tst1,tst2
       write(6,*) objtfn
                       endif
       end
!234567890
!      Written by In-Ho Lee, KRISS, April 19, 2016.
       subroutine get_emass(fname,effmass,lfault)
       USE strings, ONLY : parse,value
       implicit none
       character*280 fname
       real*8 effmass(3)
       logical lfault
       real*8 rkpt0(3),a1(3),a2(3),a3(3)
       integer iband
       real*8 dk,tmp
       integer n,i,ndr
       integer, allocatable :: iwrk(:)
       real*8, allocatable :: wrk(:)
       integer ios,nargs
       character*200 str1
       character*200 args(40)
       character*20 delims
       logical lexist
       character*280 gname

       effmass=1.d6
       lfault=.false.
       inquire(file='inp_vbm',exist=lexist)
       if(.not. lexist)then
       write(6,*) 'inp_vbm is not present.'
       lfault=.true.
       goto 111
                       endif
       inquire(file='inp_cbm',exist=lexist)
       if(.not. lexist)then
       write(6,*) 'inp_cbm is not present.'
       lfault=.true.
       goto 111
                       endif
       inquire(file='out_vbm',exist=lexist)
       if(.not. lexist)then
       write(6,*) 'out_vbm is not present.'
       lfault=.true.
       goto 111
                       endif
       inquire(file='out_cbm',exist=lexist)
       if(.not. lexist)then
       write(6,*) 'out_cbm is not present.'
       lfault=.true.
       goto 111
                       endif
  111  continue
       if(lfault)then
       write(6,*) 'serious problem'
       return
                 endif
       lfault=.false.
       open(1,file=trim(fname),form='formatted')
       do 
       read(1,*,err=711,end=799) rkpt0(1),rkpt0(2),rkpt0(3)
       read(1,*,err=711,end=799) dk
       read(1,*,err=711,end=799) iband
       read(1,*,err=711,end=799)
       read(1,*,err=711,end=799) a1(1),a1(2),a1(3)
       read(1,*,err=711,end=799) a2(1),a2(2),a2(3)
       read(1,*,err=711,end=799) a3(1),a3(2),a3(3)
       enddo
  711  continue
       lfault=.true.
  799  continue
       close(1)
       if(lfault)then
       write(6,*) 'inp file is not present'
       write(6,*) 'it is serious'
       return
                 endif
       effmass=1.d6
       lfault=.false.
       ndr=0
!      i=len_trim(fname) ; gname=trim(fname(1:i-7))//'emc_calc.log' ; gname=trim(gname)
       gname='emc_calc.log' ; gname=trim(gname)
       if(trim(fname) == 'inp_vbm') gname='out_vbm'
       if(trim(fname) == 'inp_cbm') gname='out_cbm'
       open(2,file=trim(gname),form='formatted')
       do
       read(2,'(a200)',err=911,end=999) str1
       delims=' '
       call parse(str1,delims,args,nargs)
       if(nargs >= 3)then
       if(args(1) == 'Effective' .and. args(2) == 'mass:')then
       call value(args(3),tmp,ios)
       if(ios /= 0) tmp=1.d6
       if(args(3) == '+Infinity') tmp=1.d6
       if(args(3) == '-Infinity') tmp=1.d6
       if(args(3) == 'Infinity') tmp=1.d6
       if(args(3) == 'infinity') tmp=1.d6
       if(args(3) == '-infinity') tmp=1.d6
       if(args(3) == 'NaN') tmp=1.d6
       if(args(3) == 'nan') tmp=1.d6
       if(args(3) == 'NAN') tmp=1.d6
       ndr=ndr+1
       effmass(ndr)=abs(tmp)
!      write(6,*) effmass(ndr)
       if(ndr == 3) goto 999
                                                          endif
                     endif
       enddo
  911  continue
       lfault=.true.
  999  continue
       close(2)
       if(ndr /= 3)then
       write(6,*) 'serious problem',ndr            
       lfault=.true.
                   endif
       if(.not. lfault)then
       n=3
       allocate(iwrk(n)) ; allocate(wrk(n))
       do i=1,n
       wrk(i)=effmass(i)
       enddo
       call sortnr(n,wrk,iwrk)
       do i=1,n
       effmass(i)=wrk(iwrk(i))
       enddo
       deallocate(iwrk) ; deallocate(wrk)
                       endif
       end
!234567890
!      Written by In-Ho Lee, KRISS, April 20, 2016.
       subroutine emassprincipal_inp(otname,einame,egp1,egp2,lfault)
       USE eigenvalues1, ONLY : rkpt,ne,nk,nbandi,ef
       implicit none
       logical lfault
       character*280 otname,einame
       real*8 egp1,egp2,evbm,ecbm,rkptv(3),rkptc(3),a1(3),a2(3),a3(3),scale0,cmat(3,3),vtest
       integer kvbm,kcbm,ibvbm,ibcbm
       logical lfault1,lfault2,lfault3
       character*280 fname
       integer i

       lfault=.false.
       lfault1=.false.
       lfault2=.false.
       lfault3=.false.
       egp1=-1.d0
       egp2=-1.d0
!      call read_outcar1(otname,lfault1)
       call read_efermi1(otname,lfault1)
       if(lfault1) goto 911
       call read_eigenval1(einame,lfault2)
       if(lfault2) goto 911
       call check_gapmassprep(egp1,egp2,kvbm,kcbm,ibvbm,ibcbm,evbm,ecbm,lfault3)
  911  continue
       if(lfault1 .or. lfault2 .or. lfault3) lfault=.true.
       if(lfault)then
       egp1=0.d0 ; egp2=0.d0
                 endif
       if(egp1 <= -1.d0) egp1=0.d0
       if(egp2 <= -1.d0) egp2=0.d0
!
       fname=trim(einame)
       i=len_trim(fname) ; fname=trim(fname(1:i-8))//'POSCAR' ; fname=trim(fname)
       fname='POSCAR' ; fname=trim(fname)
       open(14,file=trim(fname),form='formatted')
       read(14,*)
       read(14,*) scale0
       read(14,*) a1(1),a1(2),a1(3)
       read(14,*) a2(1),a2(2),a2(3)
       read(14,*) a3(1),a3(2),a3(3)
       if(scale0 > 0.d0)then
       a1=a1*scale0 ; a2=a2*scale0 ; a3=a3*scale0
                        endif
       if(scale0 < 0.d0)then
       cmat(1,:)=a1(:) ; cmat(2,:)=a2(:) ; cmat(3,:)=a3(:)
       vtest=(cmat(1,2)*cmat(2,3)-cmat(1,3)*cmat(2,2))*cmat(3,1) &
            +(cmat(1,3)*cmat(2,1)-cmat(1,1)*cmat(2,3))*cmat(3,2) &
            +(cmat(1,1)*cmat(2,2)-cmat(1,2)*cmat(2,1))*cmat(3,3)
       vtest=(abs(scale0)/abs(vtest))**(1.d0/3.d0)
       cmat=cmat*vtest
       a1(:)=cmat(1,:) ; a2(:)=cmat(2,:) ; a3(:)=cmat(3,:)
                        endif
       close(14)
       fname=trim(einame)
!      i=len_trim(fname) ; fname=trim(fname(1:i-8))//'els_info' ; fname=trim(fname)
       fname='els_info' ; fname=trim(fname)
       open(77,file=trim(fname),form='formatted')
       write(77,'(3i7)') ne,nk,nbandi
       write(77,'(f18.8)') ef
       write(77,'(2i7,1x,3f22.12,1x,f18.8)') ibvbm,kvbm,rkpt(1,kvbm),rkpt(2,kvbm),rkpt(3,kvbm),evbm
       write(77,'(2i7,1x,3f22.12,1x,f18.8)') ibcbm,kcbm,rkpt(1,kcbm),rkpt(2,kcbm),rkpt(3,kcbm),ecbm
       write(77,'(3f22.12)') a1(1),a1(2),a1(3)
       write(77,'(3f22.12)') a2(1),a2(2),a2(3)
       write(77,'(3f22.12)') a3(1),a3(2),a3(3)
       close(77)
       fname=trim(einame)
       i=len_trim(fname) ; fname=trim(fname(1:i-8))//'inp_vbm' ; fname=trim(fname)
       fname='inp_vbm' ; fname=trim(fname)
       rkptv(:)=rkpt(:,kvbm)
       call gen_inp(fname,ibvbm,rkptv,a1,a2,a3)
       fname=trim(einame)
       i=len_trim(fname) ; fname=trim(fname(1:i-8))//'inp_cbm' ; fname=trim(fname) 
       fname='inp_cbm' ; fname=trim(fname) 
       rkptc(:)=rkpt(:,kcbm)
       call gen_inp(fname,ibcbm,rkptc,a1,a2,a3)
       end
!234567890
!      Written by In-Ho Lee, KRISS, April 19, 2016.
       subroutine gen_inp(fname,iband,rkpt0,a1,a2,a3)
       implicit none
       character*280 fname
       integer iband
       real*8 rkpt0(3),a1(3),a2(3),a3(3)
       real*8 dk

       dk=0.01d0
       open(78,file=trim(fname),form='formatted')
       write(78,'(3f22.12)') rkpt0(1),rkpt0(2),rkpt0(3)
       write(78,'(f22.12)') dk
       write(78,'(i8)') iband
       write(78,'(a)') 'V'
       write(78,'(3f22.12)') a1(1),a1(2),a1(3)
       write(78,'(3f22.12)') a2(1),a2(2),a2(3)
       write(78,'(3f22.12)') a3(1),a3(2),a3(3)
       close(78)
       end
!234567890
!      Written by In-Ho Lee, KRISS, April 20, 2016.
       subroutine check_gapmassprep(egp1,egp2,kvbm,kcbm,ibvbm,ibcbm,evbm,ecbm,lfault3)
       USE eigenvalues1, ONLY : temper,ef,rkpt,wgt,eiv,nbandi,nk,cblowest,vbhighest,ldirect
       implicit none
       integer kvbm,kcbm,ibvbm,ibcbm
       logical lfault3
       real*8 egp1,egp2,evbm,ecbm
       integer ik,i,i1,i2
       real*8 t1,t2,test,tst1,tst2,occtol,smallest,dgap,xxr
       logical lfault
       logical lexist
       real*8, allocatable :: wrk7(:),wrk8(:)
       integer, allocatable :: iwrk7(:),iwrk8(:)

       lfault=.false.
       occtol=1.d-10
!      write(6,*) ef,' ef'
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
       write(6,*) 'something went wrong in egapmass.f90'
                                 endif
       if(i2 <1 .or. i2 > nbandi)then
       i2=1
       lfault=.true.
       write(6,*) 'something went wrong in egapmass.f90'
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
       kcbm=iwrk7(1)
       kvbm=iwrk8(nk)
       xxr=1.d20
       do i=1,nbandi
       tst1=abs(eiv(i,kcbm)-wrk7(kcbm))
       if(xxr > tst1)then
       xxr=tst1
       ibcbm=i
                     endif
       enddo
       xxr=1.d20
       do i=1,nbandi
       tst2=abs(eiv(i,kvbm)-wrk8(kvbm))
       if(xxr > tst2)then
       xxr=tst2
       ibvbm=i
                     endif
       enddo
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
       inquire(file='directbandgap',exist=lexist)
       if(lexist)then
       open(91,file='directbandgap')
       close(91,status='delete')
                 endif
                           else
       ldirect=.true.
       open(91,file='directbandgap',form='formatted')
       write(91,*) dgap,nk
       close(91)
       write(6,'(a15,2x,2f16.8,2x,i6,1x,a11,1x,a7)') 'direct band gap',dgap,smallest,nk,'+++++++++++','nk fine'
                           endif
       evbm=eiv(ibvbm,kvbm) ; ecbm=eiv(ibcbm,kcbm)
       write(6,*) 'iband,nk, rkpt0,eiv'
       write(6,'(2i7,1x,3f22.12,1x,f18.8)') ibvbm,kvbm,rkpt(1,kvbm),rkpt(2,kvbm),rkpt(3,kvbm),eiv(ibvbm,kvbm)
       write(6,'(2i7,1x,3f22.12,1x,f18.8)') ibcbm,kcbm,rkpt(1,kcbm),rkpt(2,kcbm),rkpt(3,kcbm),eiv(ibcbm,kcbm)
  119  continue
       if(allocated(wrk7)) deallocate(wrk7) 
       if(allocated(wrk8)) deallocate(wrk8) 
       if(allocated(iwrk7)) deallocate(iwrk7)
       if(allocated(iwrk8)) deallocate(iwrk8)
       if(allocated(cblowest)) deallocate(cblowest)
       if(allocated(vbhighest)) deallocate(vbhighest)
       if(lfault)then
       egp1=0.d0 ; egp2=0.d0
                 endif
       call flush(6)
       lfault3=lfault
       end
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine read_efermi1(otname,lfault4)
       USE strings, ONLY : parse,value
       USE eigenvalues1, ONLY : ef
       implicit none
       character*280 otname
       logical lfault4
       logical lfault
       integer ios,nargs
       character*200 str1
       character*200 args(40)
       character*20 delims

       lfault=.false.
       ef=1.d19
       open(81,file=trim(otname),form='formatted')
       do
       read(81,'(a200)',err=911,end=999) str1
       delims=' '
       call parse(str1,delims,args,nargs)
       if(nargs >= 3)then
       if(args(1) == 'E-fermi')then
       if(args(2) == ':')then
       call value(args(3),ef,ios)
       if(ios /= 0) ef=1.d18
                         endif
                               endif
                     endif
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
       lfault4=lfault
       end
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine read_eigenval1(einame,lfault1)
       USE eigenvalues1, ONLY : ef,rkpt,wgt,eiv,ne,nk,nbandi,temper
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
       read(81,'(i7,i7,i7)',err=911,end=999) ne,nk,nbandi
       if(nk > 100000   .or. nk <= 0) goto 911
       if(nbandi > 2000 .or. nbandi <= 0) goto 911
       if(ne > 4000     .or. ne <= 0) goto 911
       if(.not. allocated(eiv)) allocate(eiv(nbandi,nk))
       if(.not. allocated(wgt)) allocate(wgt(nk))
       if(.not. allocated(rkpt)) allocate(rkpt(3,nk))
       do ik=1,nk
!      write(6,*) ik
       read(81,*,err=911,end=999)
       read(81,*,err=911,end=999) tmpx,tmpy,tmpz,wgt(ik)
       rkpt(1,ik)=tmpx ; rkpt(2,ik)=tmpy ; rkpt(3,ik)=tmpz
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
       program egapmass_test
       USE eigenvalues1, ONLY : rkpt,eiv,wgt
       implicit none
       real*8 egp1,egp2,objtfn
       character*280 otname,einame
       character*80 cmd
       logical lfault1,lfault2,lfault
       logical lexist

!      otname='/home/ihlee/emc-1.50/fortran/test/OUTCAR'
!      einame='/home/ihlee/emc-1.50/fortran/test/EIGENVAL'
       otname='OUTCAR'
       einame='EIGENVAL'
       lfault=.false.
       inquire(file='EMCINP',exist=lexist)
       if(.not. lexist)then
       lfault2=.false.
       call get_emassprincipal(einame,objtfn,lfault2)
       if(lfault2) lfault=.true.
                       else
       lfault1=.false.
       call emassprincipal_inp(otname,einame,egp1,egp2,lfault1)
       if(lfault1) lfault=.true.
                       endif
       if(allocated(eiv)) deallocate(eiv)
       if(allocated(wgt)) deallocate(wgt)
       if(allocated(rkpt)) deallocate(rkpt)
       if(lfault)then
       cmd='echo "ERROR: while" >> stdout.log'
       cmd=trim(cmd) ; call system(cmd)
                 endif
       end program egapmass_test
