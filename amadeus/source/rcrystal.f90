!234567890
!      Written by In-Ho Lee, KRISS, October 12, 2015.
       subroutine gen_latt_site(ispg,ndeg,nspecies,nelements,symbl,sigmamatrix,voltol,refvol,qqq,lpbc,lvcs,lflag)
       implicit none
       integer ispg,ndeg,nspecies,nelements(nspecies)
       real*8 qqq(ndeg),sigmamatrix(nspecies,nspecies),voltol,refvol
       character*2 symbl(nspecies)
       logical lpbc,lvcs,lflag
       real*8 amat(3,3),a6(6),b6(6)
       real*8 voltol00
       real*8, allocatable :: d1(:,:),d2(:,:),d3(:,:)
       real*8, allocatable :: wrk11(:)
       integer, allocatable :: iwrk11(:)
       integer i,j,j1,k,ish,ij0,ifile
       real ranmar

       lflag=.true.
       ish=ndeg-6 
       do i=1,6
       b6(i)=qqq(ish+i)
       enddo
       ij0=nelements(1)
       do i=1,nspecies
       ij0=max(nelements(i),ij0)
       enddo
       allocate(d1(nspecies,ij0),d2(nspecies,ij0),d3(nspecies,ij0))
       ifile=-1
       voltol00=voltol
       if(.not. lvcs) voltol00=0.0d0
       call gensites(ispg,refvol,voltol00,sigmamatrix,symbl,nspecies,nelements,ij0,d1,d2,d3,a6,amat,ifile)
       if(lvcs)then
       do i=1,6
       qqq(ish+i)=a6(i)
       enddo
               endif
       if(.not. lvcs)then
       do i=1,6
       qqq(ish+i)=b6(i)
       enddo
                     endif
       if(lvcs)then
       k=0
       do i=1,nspecies
       do j=1,nelements(i)
       k=k+1
       qqq(3*(k-1)+1)=d1(i,j) ; qqq(3*(k-1)+2)=d2(i,j) ; qqq(3*(k-1)+3)=d3(i,j)
       enddo
       enddo
               endif
       if(.not. lvcs)then
       allocate(iwrk11(ij0)) ; allocate(wrk11(ij0))
       k=0
       do i=1,nspecies
       do j1=1,nelements(i)
       wrk11(j1)=ranmar()
       enddo
       j1=nelements(i) ; call sortnr(j1,wrk11,iwrk11)
       do j=1,nelements(i)
       k=k+1
       j1=iwrk11(j)
       qqq(3*(k-1)+1)=qqq(3*(k-1)+1)+(d1(i,j1)-0.5d0)*0.02d0
       qqq(3*(k-1)+2)=qqq(3*(k-1)+2)+(d2(i,j1)-0.5d0)*0.02d0
       qqq(3*(k-1)+3)=qqq(3*(k-1)+3)+(d3(i,j1)-0.5d0)*0.02d0
       qqq(3*(k-1)+1)=qqq(3*(k-1)+1)-anint(qqq(3*(k-1)+1))
       qqq(3*(k-1)+2)=qqq(3*(k-1)+2)-anint(qqq(3*(k-1)+2))
       qqq(3*(k-1)+3)=qqq(3*(k-1)+3)-anint(qqq(3*(k-1)+3))
       if(qqq(3*(k-1)+1) < 0.0d0) qqq(3*(k-1)+1)=qqq(3*(k-1)+1)+1.0d0
       if(qqq(3*(k-1)+2) < 0.0d0) qqq(3*(k-1)+2)=qqq(3*(k-1)+2)+1.0d0
       if(qqq(3*(k-1)+3) < 0.0d0) qqq(3*(k-1)+3)=qqq(3*(k-1)+3)+1.0d0
       enddo
       enddo
       deallocate(iwrk11) ; deallocate(wrk11)
                     endif
       deallocate(d1,d2,d3)
       write(6,'(1x,a26,1x,i4)') 'symm-introduction:mutation',ispg
       end
!234567890
!      Written by In-Ho Lee, KRISS, October 12, 2015.
       subroutine gensites(ispg0,refvol,voltol,sigmamatrix,symbl,nspecies0,nelements0,ij0,d1,d2,d3,alat,amat,ifile)
       use modmain, ONLY : atposl,avec
       implicit none
       integer ispg0,ij0,nspecies0,nelements0(nspecies0),ifile
       real*8 alat(6),amat(3,3),d1(nspecies0,ij0),d2(nspecies0,ij0),d3(nspecies0,ij0)
       real*8 sigmamatrix(nspecies0,nspecies0),refvol,voltol
       real*8 vtest,voltol00
       character*2 symbl(nspecies0)
       integer i,j,itry,mtry,ispg
       real ranmar
       logical lclash

       mtry=10000000
       mtry=10000
       mtry=10000*ij0
       j=sum(nelements0)
       if(j > 80 ) mtry=100
       voltol00=voltol
       ispg=ispg0
       if(ispg < 1 .or. ispg > 230) ispg=dble(ranmar())*230+1
       itry=0
  111  continue
       vtest=refvol*(1.0d0+voltol00*(ranmar()-0.5)*2.0d0)
       call gen_sg_lat(ispg,vtest,amat)
       call latmat(alat,amat,0)
       call preparation(ispg,alat,nspecies0,nelements0,symbl)
       call gencrystal
       call ddcheck(sigmamatrix,nspecies0,nelements0,lclash)
       itry=itry+1
       if(itry > mtry/2)then
       ispg=dble(ranmar())*230+1
       write(6,'(a14,1x,i4,1x,a2,1x,i4)') 'warning: SG id',ispg0,'->',ispg
                        endif
       if(itry > mtry)then
       write(6,'(a36,1x,i10)') 'warning: SG id, unusual case, failed',itry
       goto 222
                      endif
       if(lclash) goto 111
  222  continue
       write(6,'(i8)') itry
!------
       amat=transpose(avec) ; call latmat(alat,amat,0)
       do i=1,nspecies0
       do j=1,nelements0(i)
       d1(i,j)=atposl(1,j,i) ; d2(i,j)=atposl(2,j,i) ; d3(i,j)=atposl(3,j,i)
       if(itry > mtry)then
       d1(i,j)=ranmar() ; d2(i,j)=ranmar() ; d3(i,j)=ranmar()
                      endif
       enddo
       enddo
       call gen_poscar(ispg,nspecies0,nelements0,ifile)
       ispg0=ispg
       end 
!234567890
!      Written by In-Ho Lee, KRISS, October 12, 2015.
       subroutine ddcheck(sigmamatrix,nspecies0,nelements0,lclash)
       use modmain, ONLY : atposl,avec
       implicit none
       integer nspecies0,nelements0(nspecies0)
       real*8 sigmamatrix(nspecies0,nspecies0)
       logical lclash
       real*8 x,y,z,v1,v2,v3,dist,aa1(3),aa2(3),aa3(3),amat(3,3)
       integer n1,n2,n3,i,j,k,m,i1,j1,natot,natot_ext
       real*8, allocatable :: sites(:,:),sites_ext(:,:)
       integer, allocatable :: itype(:),itype_ext(:)
       logical lnan

!
       lnan=.false.
       do i=1,nspecies0
       do j=1,nelements0(i)
       if(isnan(atposl(1,j,i)) .or. isnan(atposl(2,j,i)) .or. isnan(atposl(2,j,i))) lnan=.true.
       enddo
       enddo
       if(isnan(avec(1,1)) .or. isnan(avec(1,2)) .or. isnan(avec(1,3))) lnan=.true.
       if(isnan(avec(2,1)) .or. isnan(avec(2,2)) .or. isnan(avec(2,3))) lnan=.true.
       if(isnan(avec(3,1)) .or. isnan(avec(3,2)) .or. isnan(avec(3,3))) lnan=.true.
       if(lnan)then
       lclash=.true.
               return
               endif
!
       lclash=.false.
       n1=2 ; n2=2 ; n3=2 ; amat=transpose(avec)
       aa1(:)=amat(1,:)*dble(n1) ; aa2(:)=amat(2,:)*dble(n2) ; aa3(:)=amat(3,:)*dble(n3)
       k=0
       do i=1,nspecies0
       do j=1,nelements0(i)
       k=k+1
       enddo
       enddo
       allocate(itype(k)) ; allocate(sites(3,k))
       natot=0
       do i=1,nspecies0
       do j=1,nelements0(i)
       natot=natot+1 ; sites(:,natot)=atposl(:,j,i) ; itype(natot)=i
       enddo
       enddo
       k=natot*(n1*n2*n3) ; allocate(itype_ext(k)) ; allocate(sites_ext(3,k))
       natot_ext=0
       do m=1,natot
       do i=0,n1-1
       do j=0,n2-1
       do k=0,n3-1
       natot_ext=natot_ext+1 ; itype_ext(natot_ext)=itype(m)
       sites_ext(1,natot_ext)=sites(1,m)+dble(i)
       sites_ext(2,natot_ext)=sites(2,m)+dble(j)
       sites_ext(3,natot_ext)=sites(3,m)+dble(k)
       enddo
       enddo
       enddo
       enddo
       do i=1,natot_ext
       sites_ext(1,i)=sites_ext(1,i)/dble(n1)
       sites_ext(2,i)=sites_ext(2,i)/dble(n2)
       sites_ext(3,i)=sites_ext(3,i)/dble(n3)
       enddo
       do i=1,natot_ext-1
       do j=i+1,natot_ext
       v1=sites_ext(1,i)-sites_ext(1,j)
       v2=sites_ext(2,i)-sites_ext(2,j)
       v3=sites_ext(3,i)-sites_ext(3,j)
       v1=v1-anint(v1) ; v2=v2-anint(v2) ; v3=v3-anint(v3)
       x=v1*aa1(1)+v2*aa2(1)+v3*aa3(1)
       y=v1*aa1(2)+v2*aa2(2)+v3*aa3(2)
       z=v1*aa1(3)+v2*aa2(3)+v3*aa3(3)
       dist=sqrt(x*x+y*y+z*z)
       i1=itype_ext(i) ; j1=itype_ext(j)
       if(dist < sigmamatrix(i1,j1))then
       lclash=.true. ; goto 999
                                    endif
       enddo
       enddo
  999  continue
       deallocate(itype,itype_ext) ; deallocate(sites,sites_ext)
       end
!234567890
!      Written by In-Ho Lee, KRISS, October 12, 2015.
       subroutine preparation(ispg,alat,nspecies0,nelements0,symbl)
       use modmain, ONLY : a,b,c,ab,ac,bc,ncell,primcell,nspecies,natoms
       use modmain, ONLY : maxatoms,maxspecies,maxwpos,spsymb,nwpos,wpos
       use modmain, ONLY : hrmg,num,schn,hall
       implicit none
       integer ispg,nspecies0,nelements0(nspecies0)
       real*8 alat(6)
       character*2 symbl(nspecies0)
       real*8 tmp,pi
       character*20 kndx(530)
       character*200 str1
       integer is,ip,kk,indexsg
       integer nwyc
       real*8, allocatable :: wrk(:,:)
       real ranmar

       nwyc=100
       if(maxwpos < nwyc)then
       write(6,*) 'increase maxwpos'
                         stop
                         endif
       allocate(wrk(3,nwyc))
       call hrmgencoding(kndx)
       do while(.true.)
       kk=dble(ranmar())*530+1 
! set the Hermann-Mauguin symbol
       hrmg=kndx(kk) ; hrmg=adjustl(hrmg)
! determine the Hall symbol from the Hermann-Mauguin symbol
       call sgsymb(hrmg,num,schn,hall)
       str1=trim(num) ; call getindex(str1,indexsg)
       if(indexsg < 1 .or. indexsg > 230)then
       write(6,*) 'something went wrong'
                                         stop
                                         endif
       if(indexsg == ispg) goto 990
       enddo
  990  continue
! set lattice vector lengths (Angstrom) and angles (rad) : this is for gencrystal
       pi=4.0d0*atan(1.0d0) ; tmp=180.0d0/pi
       a=alat(1) ; b=alat(2) ; c=alat(3) ; bc=alat(4)*tmp ; ac=alat(5)*tmp ; ab=alat(6)*tmp
       ncell(1)=1 ; ncell(2)=1 ; ncell(3)=1 ; primcell=.false.
       nspecies=nspecies0
       if(maxspecies < nspecies0)then
       write(6,*) 'increase maxspecies'
                                 stop
                                 endif
       do is=1,nspecies0
       natoms(is)=nelements0(is)
       enddo
       ip=natoms(1)
       do is=1,nspecies0
       ip=max(natoms(is),ip)
       enddo
       if(maxatoms < ip)then
       write(6,*) 'increase maxatoms'
                        stop
                        endif
       do is=1,nspecies
       spsymb(is)=symbl(is)
       call rwycpos(indexsg,wrk,nwyc)
       nwpos(is)=nwyc
       nwpos(is)=min(natoms(is),nwyc)
       do ip=1,nwpos(is)
       wpos(:,ip,is)=wrk(:,ip)
       enddo
       enddo
       deallocate(wrk)
       end
!234567890
!      Written by In-Ho Lee, KRISS, October 12, 2015.
       subroutine getindex(str1,iidd)
       use strings
       implicit none
       integer iidd
       character*200 str1
       integer ios,nargs
       character*200 args(40)
       character*20 delims

       delims=':'
       call parse(str1,delims,args,nargs)
       if(nargs >= 1)then
       call value(args(1),iidd,ios)
                     endif
       if(iidd < 1 .or. iidd > 230)then
       write(6,*) 'something went wrong'
                                   stop
                                   endif
       end
!234567890
!      Written by In-Ho Lee, KRISS, October 12, 2015.
       subroutine gen_poscar(ispg,nspecies0,nelements0,ifile)
       use modmain, ONLY : nspecies,natoms,atposl,avec,spsymb
       implicit none
       integer ifile,ispg,nspecies0,nelements0(nspecies0)
       real*8 alat(6),amat(3,3),vtest,pi
       integer i,j,isize
       character*280 fname

       if(nspecies0 /= nspecies)then
       write(6,*) 'something went wrong'
                                stop
                                endif
       write(6,'(10(1x,a2))') (spsymb(i),i=1,nspecies)
       write(6,'(10(1x,i6))') (natoms(i),i=1,nspecies)
       amat=transpose(avec) ; call latmat(alat,amat,0)
       vtest=(amat(1,2)*amat(2,3)-amat(1,3)*amat(2,2))*amat(3,1)  &
            +(amat(1,3)*amat(2,1)-amat(1,1)*amat(2,3))*amat(3,2)  &
            +(amat(1,1)*amat(2,2)-amat(1,2)*amat(2,1))*amat(3,3)
       vtest=abs(vtest)
       write(6,'(f18.8)') vtest
       pi=4.0d0*atan(1.0d0)
       fname='POSCAR_tmp'
       if(ifile > 0)then
       isize=7
       call xnumeral(ifile,fname,isize)
       fname='./deposit/POSCAR_'//trim(fname) ; fname=trim(fname)
                    endif
       open(11,file=trim(fname),form='formatted')
       write(11,'(i4,1x,6f10.4,1x,f18.4)') ispg,(alat(i),i=1,3),(alat(i)*180.0d0/pi,i=4,6),vtest
       write(11,*) '1.0'
       write(11,'(3f22.13)') amat(1,1),amat(1,2),amat(1,3)
       write(11,'(3f22.13)') amat(2,1),amat(2,2),amat(2,3)
       write(11,'(3f22.13)') amat(3,1),amat(3,2),amat(3,3)
       write(11,'(20(2x,a2,1x))') (spsymb(i),i=1,nspecies)
       write(11,'(20(i4,1x))') (nelements0(i),i=1,nspecies)
       write(11,'(a6)') "Direct"
       do i=1,nspecies0
       do j=1,nelements0(i)
       write(11,'(3f22.13)') atposl(1,j,i),atposl(2,j,i),atposl(3,j,i)
       enddo
       enddo
       close(11)
       end
!234567890
!      Written by In-Ho Lee, KRISS, October 12, 2015.
       subroutine hrmgencoding(kndx)
       character*20 kndx(530)

       kndx(1)='P1'
       kndx(2)='P-1'
       kndx(3)='P2:b'
       kndx(4)='P2:c'
       kndx(5)='P2:a'
       kndx(6)='P21:b'
       kndx(7)='P21:c'
       kndx(8)='P21:a'
       kndx(9)='C2:b1'
       kndx(10)='C2:b2'
       kndx(11)='C2:b3'
       kndx(12)='C2:c1'
       kndx(13)='C2:c2'
       kndx(14)='C2:c3'
       kndx(15)='C2:a1'
       kndx(16)='C2:a2'
       kndx(17)='C2:a3'
       kndx(18)='Pm:b'
       kndx(19)='Pm:c'
       kndx(20)='Pm:a'
       kndx(21)='Pc:b1'
       kndx(22)='Pc:b2'
       kndx(23)='Pc:b3'
       kndx(24)='Pc:c1'
       kndx(25)='Pc:c2'
       kndx(26)='Pc:c3'
       kndx(27)='Pc:a1'
       kndx(28)='Pc:a2'
       kndx(29)='Pc:a3'
       kndx(30)='Cm:b1'
       kndx(31)='Cm:b2'
       kndx(32)='Cm:b3'
       kndx(33)='Cm:c1'
       kndx(34)='Cm:c2'
       kndx(35)='Cm:c3'
       kndx(36)='Cm:a1'
       kndx(37)='Cm:a2'
       kndx(38)='Cm:a3'
       kndx(39)='Cc:b1'
       kndx(40)='Cc:b2'
       kndx(41)='Cc:b3'
       kndx(42)='Cc:-b1'
       kndx(43)='Cc:-b2'
       kndx(44)='Cc:-b3'
       kndx(45)='Cc:c1'
       kndx(46)='Cc:c2'
       kndx(47)='Cc:c3'
       kndx(48)='Cc:-c1'
       kndx(49)='Cc:-c2'
       kndx(50)='Cc:-c3'
       kndx(51)='Cc:a1'
       kndx(52)='Cc:a2'
       kndx(53)='Cc:a3'
       kndx(54)='Cc:-a1'
       kndx(55)='Cc:-a2'
       kndx(56)='Cc:-a3'
       kndx(57)='P2/m:b'
       kndx(58)='P2/m:c'
       kndx(59)='P2/m:a'
       kndx(60)='P21/m:b'
       kndx(61)='P21/m:c'
       kndx(62)='P21/m:a'
       kndx(63)='C2/m:b1'
       kndx(64)='C2/m:b2'
       kndx(65)='C2/m:b3'
       kndx(66)='C2/m:c1'
       kndx(67)='C2/m:c2'
       kndx(68)='C2/m:c3'
       kndx(69)='C2/m:a1'
       kndx(70)='C2/m:a2'
       kndx(71)='C2/m:a3'
       kndx(72)='P2/c:b1'
       kndx(73)='P2/c:b2'
       kndx(74)='P2/c:b3'
       kndx(75)='P2/c:c1'
       kndx(76)='P2/c:c2'
       kndx(77)='P2/c:c3'
       kndx(78)='P2/c:a1'
       kndx(79)='P2/c:a2'
       kndx(80)='P2/c:a3'
       kndx(81)='P21/c:b1'
       kndx(82)='P21/c:b2'
       kndx(83)='P21/c:b3'
       kndx(84)='P21/c:c1'
       kndx(85)='P21/c:c2'
       kndx(86)='P21/c:c3'
       kndx(87)='P21/c:a1'
       kndx(88)='P21/c:a2'
       kndx(89)='P21/c:a3'
       kndx(90)='C2/c:b1'
       kndx(91)='C2/c:b2'
       kndx(92)='C2/c:b3'
       kndx(93)='C2/c:-b1'
       kndx(94)='C2/c:-b2'
       kndx(95)='C2/c:-b3'
       kndx(96)='C2/c:c1'
       kndx(97)='C2/c:c2'
       kndx(98)='C2/c:c3'
       kndx(99)='C2/c:-c1'
       kndx(100)='C2/c:-c2'
       kndx(101)='C2/c:-c3'
       kndx(102)='C2/c:a1'
       kndx(103)='C2/c:a2'
       kndx(104)='C2/c:a3'
       kndx(105)='C2/c:-a1'
       kndx(106)='C2/c:-a2'
       kndx(107)='C2/c:-a3'
       kndx(108)='P222'
       kndx(109)='P2221'
       kndx(110)='P2122'
       kndx(111)='P2212'
       kndx(112)='P21212'
       kndx(113)='P22121'
       kndx(114)='P21221'
       kndx(115)='P212121'
       kndx(116)='C2221'
       kndx(117)='A2122'
       kndx(118)='B2212'
       kndx(119)='C222'
       kndx(120)='A222'
       kndx(121)='B222'
       kndx(122)='F222'
       kndx(123)='I222'
       kndx(124)='I212121'
       kndx(125)='Pmm2'
       kndx(126)='P2mm'
       kndx(127)='Pm2m'
       kndx(128)='Pmc21'
       kndx(129)='Pcm21'
       kndx(130)='P21ma'
       kndx(131)='P21am'
       kndx(132)='Pb21m'
       kndx(133)='Pm21b'
       kndx(134)='Pcc2'
       kndx(135)='P2aa'
       kndx(136)='Pb2b'
       kndx(137)='Pma2'
       kndx(138)='Pbm2'
       kndx(139)='P2mb'
       kndx(140)='P2cm'
       kndx(141)='Pc2m'
       kndx(142)='Pm2a'
       kndx(143)='Pca21'
       kndx(144)='Pbc21'
       kndx(145)='P21ab'
       kndx(146)='P21ca'
       kndx(147)='Pc21b'
       kndx(148)='Pb21a'
       kndx(149)='Pnc2'
       kndx(150)='Pcn2'
       kndx(151)='P2na'
       kndx(152)='P2an'
       kndx(153)='Pb2n'
       kndx(154)='Pn2b'
       kndx(155)='Pmn21'
       kndx(156)='Pnm21'
       kndx(157)='P21mn'
       kndx(158)='P21nm'
       kndx(159)='Pn21m'
       kndx(160)='Pm21n'
       kndx(161)='Pba2'
       kndx(162)='P2cb'
       kndx(163)='Pc2a'
       kndx(164)='Pna21'
       kndx(165)='Pbn21'
       kndx(166)='P21nb'
       kndx(167)='P21cn'
       kndx(168)='Pc21n'
       kndx(169)='Pn21a'
       kndx(170)='Pnn2'
       kndx(171)='P2nn'
       kndx(172)='Pn2n'
       kndx(173)='Cmm2'
       kndx(174)='A2mm'
       kndx(175)='Bm2m'
       kndx(176)='Cmc21'
       kndx(177)='Ccm21'
       kndx(178)='A21ma'
       kndx(179)='A21am'
       kndx(180)='Bb21m'
       kndx(181)='Bm21b'
       kndx(182)='Ccc2'
       kndx(183)='A2aa'
       kndx(184)='Bb2b'
       kndx(185)='Amm2'
       kndx(186)='Bmm2'
       kndx(187)='B2mm'
       kndx(188)='C2mm'
       kndx(189)='Cm2m'
       kndx(190)='Am2m'
       kndx(191)='Abm2'
       kndx(192)='Bma2'
       kndx(193)='B2cm'
       kndx(194)='C2mb'
       kndx(195)='Cm2a'
       kndx(196)='Ac2m'
       kndx(197)='Ama2'
       kndx(198)='Bbm2'
       kndx(199)='B2mb'
       kndx(200)='C2cm'
       kndx(201)='Cc2m'
       kndx(202)='Am2a'
       kndx(203)='Aba2'
       kndx(204)='Bba2'
       kndx(205)='B2cb'
       kndx(206)='C2cb'
       kndx(207)='Cc2a'
       kndx(208)='Ac2a'
       kndx(209)='Fmm2'
       kndx(210)='F2mm'
       kndx(211)='Fm2m'
       kndx(212)='Fdd2'
       kndx(213)='F2dd'
       kndx(214)='Fd2d'
       kndx(215)='Imm2'
       kndx(216)='I2mm'
       kndx(217)='Im2m'
       kndx(218)='Iba2'
       kndx(219)='I2cb'
       kndx(220)='Ic2a'
       kndx(221)='Ima2'
       kndx(222)='Ibm2'
       kndx(223)='I2mb'
       kndx(224)='I2cm'
       kndx(225)='Ic2m'
       kndx(226)='Im2a'
       kndx(227)='Pmmm'
       kndx(228)='Pnnn:1'
       kndx(229)='Pnnn:2'
       kndx(230)='Pccm'
       kndx(231)='Pmaa'
       kndx(232)='Pbmb'
       kndx(233)='Pban:1'
       kndx(234)='Pban:2'
       kndx(235)='Pncb:1'
       kndx(236)='Pncb:2'
       kndx(237)='Pcna:1'
       kndx(238)='Pcna:2'
       kndx(239)='Pmma'
       kndx(240)='Pmmb'
       kndx(241)='Pbmm'
       kndx(242)='Pcmm'
       kndx(243)='Pmcm'
       kndx(244)='Pmam'
       kndx(245)='Pnna'
       kndx(246)='Pnnb'
       kndx(247)='Pbnn'
       kndx(248)='Pcnn'
       kndx(249)='Pncn'
       kndx(250)='Pnan'
       kndx(251)='Pmna'
       kndx(252)='Pnmb'
       kndx(253)='Pbmn'
       kndx(254)='Pcnm'
       kndx(255)='Pncm'
       kndx(256)='Pman'
       kndx(257)='Pcca'
       kndx(258)='Pccb'
       kndx(259)='Pbaa'
       kndx(260)='Pcaa'
       kndx(261)='Pbcb'
       kndx(262)='Pbab'
       kndx(263)='Pbam'
       kndx(264)='Pmcb'
       kndx(265)='Pcma'
       kndx(266)='Pccn'
       kndx(267)='Pnaa'
       kndx(268)='Pbnb'
       kndx(269)='Pbcm'
       kndx(270)='Pcam'
       kndx(271)='Pmca'
       kndx(272)='Pmab'
       kndx(273)='Pbma'
       kndx(274)='Pcmb'
       kndx(275)='Pnnm'
       kndx(276)='Pmnn'
       kndx(277)='Pnmn'
       kndx(278)='Pmmn:1'
       kndx(279)='Pmmn:2'
       kndx(280)='Pnmm:1'
       kndx(281)='Pnmm:2'
       kndx(282)='Pmnm:1'
       kndx(283)='Pmnm:2'
       kndx(284)='Pbcn'
       kndx(285)='Pcan'
       kndx(286)='Pnca'
       kndx(287)='Pnab'
       kndx(288)='Pbna'
       kndx(289)='Pcnb'
       kndx(290)='Pbca'
       kndx(291)='Pcab'
       kndx(292)='Pnma'
       kndx(293)='Pmnb'
       kndx(294)='Pbnm'
       kndx(295)='Pcmn'
       kndx(296)='Pmcn'
       kndx(297)='Pnam'
       kndx(298)='Cmcm'
       kndx(299)='Ccmm'
       kndx(300)='Amma'
       kndx(301)='Amam'
       kndx(302)='Bbmm'
       kndx(303)='Bmmb'
       kndx(304)='Cmca'
       kndx(305)='Ccmb'
       kndx(306)='Abma'
       kndx(307)='Acam'
       kndx(308)='Bbcm'
       kndx(309)='Bmab'
       kndx(310)='Cmmm'
       kndx(311)='Ammm'
       kndx(312)='Bmmm'
       kndx(313)='Cccm'
       kndx(314)='Amaa'
       kndx(315)='Bbmb'
       kndx(316)='Cmma'
       kndx(317)='Cmmb'
       kndx(318)='Abmm'
       kndx(319)='Acmm'
       kndx(320)='Bmcm'
       kndx(321)='Bmam'
       kndx(322)='Ccca:1'
       kndx(323)='Ccca:2'
       kndx(324)='Cccb:1'
       kndx(325)='Cccb:2'
       kndx(326)='Abaa:1'
       kndx(327)='Abaa:2'
       kndx(328)='Acaa:1'
       kndx(329)='Acaa:2'
       kndx(330)='Bbcb:1'
       kndx(331)='Bbcb:2'
       kndx(332)='Bbab:1'
       kndx(333)='Bbab:2'
       kndx(334)='Fmmm'
       kndx(335)='Fddd:1'
       kndx(336)='Fddd:2'
       kndx(337)='Immm'
       kndx(338)='Ibam'
       kndx(339)='Imcb'
       kndx(340)='Icma'
       kndx(341)='Ibca'
       kndx(342)='Icab'
       kndx(343)='Imma'
       kndx(344)='Immb'
       kndx(345)='Ibmm'
       kndx(346)='Icmm'
       kndx(347)='Imcm'
       kndx(348)='Imam'
       kndx(349)='P4'
       kndx(350)='P41'
       kndx(351)='P42'
       kndx(352)='P43'
       kndx(353)='I4'
       kndx(354)='I41'
       kndx(355)='P-4'
       kndx(356)='I-4'
       kndx(357)='P4/m'
       kndx(358)='P42/m'
       kndx(359)='P4/n:1'
       kndx(360)='P4/n:2'
       kndx(361)='P42/n:1'
       kndx(362)='P42/n:2'
       kndx(363)='I4/m'
       kndx(364)='I41/a:1'
       kndx(365)='I41/a:2'
       kndx(366)='P422'
       kndx(367)='P4212'
       kndx(368)='P4122'
       kndx(369)='P41212'
       kndx(370)='P4222'
       kndx(371)='P42212'
       kndx(372)='P4322'
       kndx(373)='P43212'
       kndx(374)='I422'
       kndx(375)='I4122'
       kndx(376)='P4mm'
       kndx(377)='P4bm'
       kndx(378)='P42cm'
       kndx(379)='P42nm'
       kndx(380)='P4cc'
       kndx(381)='P4nc'
       kndx(382)='P42mc'
       kndx(383)='P42bc'
       kndx(384)='I4mm'
       kndx(385)='I4cm'
       kndx(386)='I41md'
       kndx(387)='I41cd'
       kndx(388)='P-42m'
       kndx(389)='P-42c'
       kndx(390)='P-421m'
       kndx(391)='P-421c'
       kndx(392)='P-4m2'
       kndx(393)='P-4c2'
       kndx(394)='P-4b2'
       kndx(395)='P-4n2'
       kndx(396)='I-4m2'
       kndx(397)='I-4c2'
       kndx(398)='I-42m'
       kndx(399)='I-42d'
       kndx(400)='P4/mmm'
       kndx(401)='P4/mcc'
       kndx(402)='P4/nbm:1'
       kndx(403)='P4/nbm:2'
       kndx(404)='P4/nnc:1'
       kndx(405)='P4/nnc:2'
       kndx(406)='P4/mbm'
       kndx(407)='P4/mnc'
       kndx(408)='P4/nmm:1'
       kndx(409)='P4/nmm:2'
       kndx(410)='P4/ncc:1'
       kndx(411)='P4/ncc:2'
       kndx(412)='P42/mmc'
       kndx(413)='P42/mcm'
       kndx(414)='P42/nbc:1'
       kndx(415)='P42/nbc:2'
       kndx(416)='P42/nnm:1'
       kndx(417)='P42/nnm:2'
       kndx(418)='P42/mbc'
       kndx(419)='P42/mnm'
       kndx(420)='P42/nmc:1'
       kndx(421)='P42/nmc:2'
       kndx(422)='P42/ncm:1'
       kndx(423)='P42/ncm:2'
       kndx(424)='I4/mmm'
       kndx(425)='I4/mcm'
       kndx(426)='I41/amd:1'
       kndx(427)='I41/amd:2'
       kndx(428)='I41/acd:1'
       kndx(429)='I41/acd:2'
       kndx(430)='P3'
       kndx(431)='P31'
       kndx(432)='P32'
       kndx(433)='R3:H'
       kndx(434)='R3:R'
       kndx(435)='P-3'
       kndx(436)='R-3:H'
       kndx(437)='R-3:R'
       kndx(438)='P312'
       kndx(439)='P321'
       kndx(440)='P3112'
       kndx(441)='P3121'
       kndx(442)='P3212'
       kndx(443)='P3221'
       kndx(444)='R32:H'
       kndx(445)='R32:R'
       kndx(446)='P3m1'
       kndx(447)='P31m'
       kndx(448)='P3c1'
       kndx(449)='P31c'
       kndx(450)='R3m:H'
       kndx(451)='R3m:R'
       kndx(452)='R3c:H'
       kndx(453)='R3c:R'
       kndx(454)='P-31m'
       kndx(455)='P-31c'
       kndx(456)='P-3m1'
       kndx(457)='P-3c1'
       kndx(458)='R-3m:H'
       kndx(459)='R-3m:R'
       kndx(460)='R-3c:H'
       kndx(461)='R-3c:R'
       kndx(462)='P6'
       kndx(463)='P61'
       kndx(464)='P65'
       kndx(465)='P62'
       kndx(466)='P64'
       kndx(467)='P63'
       kndx(468)='P-6'
       kndx(469)='P6/m'
       kndx(470)='P63/m'
       kndx(471)='P622'
       kndx(472)='P6122'
       kndx(473)='P6522'
       kndx(474)='P6222'
       kndx(475)='P6422'
       kndx(476)='P6322'
       kndx(477)='P6mm'
       kndx(478)='P6cc'
       kndx(479)='P63cm'
       kndx(480)='P63mc'
       kndx(481)='P-6m2'
       kndx(482)='P-6c2'
       kndx(483)='P-62m'
       kndx(484)='P-62c'
       kndx(485)='P6/mmm'
       kndx(486)='P6/mcc'
       kndx(487)='P63/mcm'
       kndx(488)='P63/mmc'
       kndx(489)='P23'
       kndx(490)='F23'
       kndx(491)='I23'
       kndx(492)='P213'
       kndx(493)='I213'
       kndx(494)='Pm-3'
       kndx(495)='Pn-3:1'
       kndx(496)='Pn-3:2'
       kndx(497)='Fm-3'
       kndx(498)='Fd-3:1'
       kndx(499)='Fd-3:2'
       kndx(500)='Im-3'
       kndx(501)='Pa-3'
       kndx(502)='Ia-3'
       kndx(503)='P432'
       kndx(504)='P4232'
       kndx(505)='F432'
       kndx(506)='F4132'
       kndx(507)='I432'
       kndx(508)='P4332'
       kndx(509)='P4132'
       kndx(510)='I4132'
       kndx(511)='P-43m'
       kndx(512)='F-43m'
       kndx(513)='I-43m'
       kndx(514)='P-43n'
       kndx(515)='F-43c'
       kndx(516)='I-43d'
       kndx(517)='Pm-3m'
       kndx(518)='Pn-3n:1'
       kndx(519)='Pn-3n:2'
       kndx(520)='Pm-3n'
       kndx(521)='Pn-3m:1'
       kndx(522)='Pn-3m:2'
       kndx(523)='Fm-3m'
       kndx(524)='Fm-3c'
       kndx(525)='Fd-3m:1'
       kndx(526)='Fd-3m:2'
       kndx(527)='Fd-3c:1'
       kndx(528)='Fd-3c:2'
       kndx(529)='Im-3m'
       kndx(530)='Ia-3d'
       end
!234567890
!      Written by In-Ho Lee, KRISS, November 17, 2017.
       real*8 function atomicmass(symbl)
       implicit none
       character*2 symbl
       real*8 arr(109)
       character*2 symbl0(109)
       integer i

       atomicmass=1.0079d0
       arr(  1)=  1.0079   ; symbl0(  1)='H'
       arr(  2)=  4.0026   ; symbl0(  2)='He'
       arr(  3)=  6.941    ; symbl0(  3)='Li'
       arr(  4)=  9.0122   ; symbl0(  4)='Be'
       arr(  5)=  10.811   ; symbl0(  5)='B'
       arr(  6)=  12.0107  ; symbl0(  6)='C'
       arr(  7)=  14.0067  ; symbl0(  7)='N'
       arr(  8)=  15.9994  ; symbl0(  8)='O'
       arr(  9)=  18.9984  ; symbl0(  9)='F'
       arr( 10)=  20.1797  ; symbl0( 10)='Ne'
       arr( 11)=  22.9897  ; symbl0( 11)='Na'
       arr( 12)=  24.305   ; symbl0( 12)='Mg'
       arr( 13)=  26.9815  ; symbl0( 13)='Al'
       arr( 14)=  28.0855  ; symbl0( 14)='Si'
       arr( 15)=  30.9738  ; symbl0( 15)='P'
       arr( 16)=  32.065   ; symbl0( 16)='S'
       arr( 17)=  35.453   ; symbl0( 17)='Cl'
       arr( 18)=  39.948   ; symbl0( 18)='Ar'
       arr( 19)=  39.0983  ; symbl0( 19)='K'
       arr( 20)=  40.078   ; symbl0( 20)='Ca'
       arr( 21)=  44.9559  ; symbl0( 21)='Sc'
       arr( 22)=  47.867   ; symbl0( 22)='Ti'
       arr( 23)=  50.9415  ; symbl0( 23)='V'
       arr( 24)=  51.9961  ; symbl0( 24)='Cr'
       arr( 25)=  54.938   ; symbl0( 25)='Mn'
       arr( 26)=  55.845   ; symbl0( 26)='Fe'
       arr( 27)=  58.9332  ; symbl0( 27)='Co'
       arr( 28)=  58.6934  ; symbl0( 28)='Ni'
       arr( 29)=  63.546   ; symbl0( 29)='Cu'
       arr( 30)=  65.39    ; symbl0( 30)='Zn'
       arr( 31)=  69.723   ; symbl0( 31)='Ga'
       arr( 32)=  72.64    ; symbl0( 32)='Ge'
       arr( 33)=  74.9216  ; symbl0( 33)='As'
       arr( 34)=  78.96    ; symbl0( 34)='Se'
       arr( 35)=  79.904   ; symbl0( 35)='Br'
       arr( 36)=  83.8     ; symbl0( 36)='Kr'
       arr( 37)=  85.4678  ; symbl0( 37)='Rb'
       arr( 38)=  87.62    ; symbl0( 38)='Sr'
       arr( 39)=  88.9059  ; symbl0( 39)='Y'
       arr( 40)=  91.224   ; symbl0( 40)='Zr'
       arr( 41)=  92.9064  ; symbl0( 41)='Nb'
       arr( 42)=  95.94    ; symbl0( 42)='Mo'
       arr( 43)=  98.00    ; symbl0( 43)='Tc'
       arr( 44)=  101.07   ; symbl0( 44)='Ru'
       arr( 45)=  102.9055 ; symbl0( 45)='Rh'
       arr( 46)=  106.42   ; symbl0( 46)='Pd'
       arr( 47)=  107.8682 ; symbl0( 47)='Ag'
       arr( 48)=  112.411  ; symbl0( 48)='Cd'
       arr( 49)=  114.818  ; symbl0( 49)='In'
       arr( 50)=  118.71   ; symbl0( 50)='Sn'
       arr( 51)=  121.76   ; symbl0( 51)='Sb'
       arr( 52)=  127.6    ; symbl0( 52)='Te'
       arr( 53)=  126.9045 ; symbl0( 53)='I'
       arr( 54)=  131.293  ; symbl0( 54)='Xe'
       arr( 55)=  132.9055 ; symbl0( 55)='Cs'
       arr( 56)=  137.327  ; symbl0( 56)='Ba'
       arr( 57)=  138.9055 ; symbl0( 57)='La'
       arr( 58)=  140.116  ; symbl0( 58)='Ce'
       arr( 59)=  140.9077 ; symbl0( 59)='Pr'
       arr( 60)=  144.24   ; symbl0( 60)='Nd'
       arr( 61)=  145.00   ; symbl0( 61)='Pm'
       arr( 62)=  150.36   ; symbl0( 62)='Sm'
       arr( 63)=  151.964  ; symbl0( 63)='Eu'
       arr( 64)=  157.25   ; symbl0( 64)='Gd'
       arr( 65)=  158.9253 ; symbl0( 65)='Tb'
       arr( 66)=  162.5    ; symbl0( 66)='Dy'
       arr( 67)=  164.9303 ; symbl0( 67)='Ho'
       arr( 68)=  167.259  ; symbl0( 68)='Er'
       arr( 69)=  168.9342 ; symbl0( 69)='Tm'
       arr( 70)=  173.04   ; symbl0( 70)='Yb'
       arr( 71)=  174.967  ; symbl0( 71)='Lu'
       arr( 72)=  178.49   ; symbl0( 72)='Hf'
       arr( 73)=  180.9479 ; symbl0( 73)='Ta'
       arr( 74)=  183.84   ; symbl0( 74)='W'
       arr( 75)=  186.207  ; symbl0( 75)='Re'
       arr( 76)=  190.23   ; symbl0( 76)='Os'
       arr( 77)=  192.217  ; symbl0( 77)='Ir'
       arr( 78)=  195.078  ; symbl0( 78)='Pt'
       arr( 79)=  196.9665 ; symbl0( 79)='Au'
       arr( 80)=  200.59   ; symbl0( 80)='Hg'
       arr( 81)=  204.3833 ; symbl0( 81)='Tl'
       arr( 82)=  207.2    ; symbl0( 82)='Pb'
       arr( 83)=  208.9804 ; symbl0( 83)='Bi'
       arr( 84)=  209.     ; symbl0( 84)='Po'
       arr( 85)=  210.     ; symbl0( 85)='At'
       arr( 86)=  222.     ; symbl0( 86)='Rn'
       arr( 87)=  223.     ; symbl0( 87)='Fr'
       arr( 88)=  226.     ; symbl0( 88)='Ra'
       arr( 89)=  227.     ; symbl0( 89)='Ac'
       arr( 90)=  232.0381 ; symbl0( 90)='Th'
       arr( 91)=  231.0359 ; symbl0( 91)='Pa'
       arr( 92)=  238.0289 ; symbl0( 92)='U'
       arr( 93)=  237.     ; symbl0( 93)='Np'
       arr( 94)=  244.     ; symbl0( 94)='Pu'
       arr( 95)=  243.     ; symbl0( 95)='Am'
       arr( 96)=  247.     ; symbl0( 96)='Cm'
       arr( 97)=  247.     ; symbl0( 97)='Bk'
       arr( 98)=  251.     ; symbl0( 98)='Cf'
       arr( 99)=  252.     ; symbl0( 99)='Es'
       arr(100)=  257.     ; symbl0(100)='Fm'
       arr(101)=  258.     ; symbl0(101)='Md'
       arr(102)=  259.     ; symbl0(102)='No'
       arr(103)=  262.     ; symbl0(103)='Lr'
       arr(104)=  261.     ; symbl0(104)='Rf'
       arr(105)=  262.     ; symbl0(105)='Db'
       arr(106)=  266.     ; symbl0(106)='Sg'
       arr(107)=  264.     ; symbl0(107)='Bh'
       arr(108)=  277.     ; symbl0(108)='Hs'
       arr(109)=  268.     ; symbl0(109)='Mt'
       do i=1,109
       if(trim(adjustl(symbl)) == trim(adjustl(symbl0(i))))then
       atomicmass=arr(i)
                                                           exit
                                                           endif
       enddo
       end
!234567890
!      https://en.wikipedia.org/wiki/Covalent_radius
!      Written by In-Ho Lee, KRISS, November 16, 2015.
       real*8 function covlaentrr(symbol)
       implicit none
       character*2 symbol
       real*8 rr
    
       rr=1.0d0
       if(trim(symbol) == 'H')  rr=0.31d0
       if(trim(symbol) == 'He') rr=0.28d0
       if(trim(symbol) == 'Li') rr=1.28d0
       if(trim(symbol) == 'Be') rr=0.96d0
       if(trim(symbol) == 'B')  rr=0.84d0
       if(trim(symbol) == 'C')  rr=0.69d0
       if(trim(symbol) == 'N')  rr=0.71d0
       if(trim(symbol) == 'O')  rr=0.66d0
       if(trim(symbol) == 'F')  rr=0.57d0
       if(trim(symbol) == 'Ne') rr=0.58d0
       if(trim(symbol) == 'Na') rr=1.66d0
       if(trim(symbol) == 'Mg') rr=1.41d0
       if(trim(symbol) == 'Al') rr=1.21d0
       if(trim(symbol) == 'Si') rr=1.11d0
       if(trim(symbol) == 'P')  rr=1.07d0
       if(trim(symbol) == 'S')  rr=1.05d0
       if(trim(symbol) == 'Cl') rr=1.02d0
       if(trim(symbol) == 'Ar') rr=1.06d0
       if(trim(symbol) == 'K')  rr=2.03d0
       if(trim(symbol) == 'Ca') rr=1.76d0
       if(trim(symbol) == 'Sc') rr=1.70d0
       if(trim(symbol) == 'Ti') rr=1.60d0
       if(trim(symbol) == 'V')  rr=1.53d0
       if(trim(symbol) == 'Cr') rr=1.39d0
       if(trim(symbol) == 'Mn') rr=1.39d0
       if(trim(symbol) == 'Fe') rr=1.32d0
       if(trim(symbol) == 'Co') rr=1.26d0
       if(trim(symbol) == 'Ni') rr=1.24d0
       if(trim(symbol) == 'Cu') rr=1.32d0
       if(trim(symbol) == 'Zn') rr=1.22d0
       if(trim(symbol) == 'Ga') rr=1.22d0
       if(trim(symbol) == 'Ge') rr=1.20d0
       if(trim(symbol) == 'As') rr=1.19d0
       if(trim(symbol) == 'Se') rr=1.20d0
       if(trim(symbol) == 'Br') rr=1.20d0
       if(trim(symbol) == 'Kr') rr=1.16d0
       if(trim(symbol) == 'Rb') rr=2.20d0
       if(trim(symbol) == 'Sr') rr=1.95d0
       if(trim(symbol) == 'Y')  rr=1.90d0
       if(trim(symbol) == 'Zr') rr=1.75d0
       if(trim(symbol) == 'Nb') rr=1.64d0
       if(trim(symbol) == 'Mo') rr=1.54d0
       if(trim(symbol) == 'Tc') rr=1.47d0
       if(trim(symbol) == 'Ru') rr=1.46d0
       if(trim(symbol) == 'Rh') rr=1.42d0
       if(trim(symbol) == 'Pd') rr=1.39d0
       if(trim(symbol) == 'Ag') rr=1.45d0
       if(trim(symbol) == 'Cd') rr=1.44d0
       if(trim(symbol) == 'In') rr=1.42d0
       if(trim(symbol) == 'Sn') rr=1.39d0
       if(trim(symbol) == 'Sb') rr=1.39d0
       if(trim(symbol) == 'Te') rr=1.38d0
       if(trim(symbol) == 'I' ) rr=1.39d0
       if(trim(symbol) == 'Xe') rr=1.40d0
       if(trim(symbol) == 'Cs') rr=2.44d0
       if(trim(symbol) == 'Ba') rr=2.15d0
       if(trim(symbol) == 'La') rr=2.07d0
       if(trim(symbol) == 'Lu') rr=1.87d0
       if(trim(symbol) == 'Hf') rr=1.75d0
       if(trim(symbol) == 'Ta') rr=1.70d0
       if(trim(symbol) == 'W')  rr=1.62d0
       if(trim(symbol) == 'Re') rr=1.51d0
       if(trim(symbol) == 'Os') rr=1.44d0
       if(trim(symbol) == 'Ir') rr=1.41d0
       if(trim(symbol) == 'Pt') rr=1.36d0
       if(trim(symbol) == 'Au') rr=1.36d0
       if(trim(symbol) == 'Hg') rr=1.32d0
       if(trim(symbol) == 'Tl') rr=1.45d0
       if(trim(symbol) == 'Pb') rr=1.46d0
       if(trim(symbol) == 'Bi') rr=1.48d0
       if(trim(symbol) == 'Po') rr=1.40d0
       if(trim(symbol) == 'At') rr=1.50d0
       if(trim(symbol) == 'Rn') rr=1.50d0
       if(trim(symbol) == 'Fr') rr=2.60d0
       if(trim(symbol) == 'Ra') rr=2.21d0
       if(trim(symbol) == 'Ce') rr=2.04d0
       if(trim(symbol) == 'Pr') rr=2.03d0
       if(trim(symbol) == 'Nd') rr=2.01d0
       if(trim(symbol) == 'Pm') rr=1.99d0
       if(trim(symbol) == 'Sm') rr=1.98d0
       if(trim(symbol) == 'Eu') rr=1.98d0
       if(trim(symbol) == 'Gd') rr=1.96d0
       if(trim(symbol) == 'Tb') rr=1.94d0
       if(trim(symbol) == 'Dy') rr=1.92d0
       if(trim(symbol) == 'Ho') rr=1.92d0
       if(trim(symbol) == 'Er') rr=1.89d0
       if(trim(symbol) == 'Tm') rr=1.90d0
       if(trim(symbol) == 'Yb') rr=1.87d0
       if(trim(symbol) == 'Ac') rr=2.15d0
       if(trim(symbol) == 'Th') rr=2.06d0
       if(trim(symbol) == 'Pa') rr=2.00d0
       if(trim(symbol) == 'U')  rr=1.96d0
       if(trim(symbol) == 'Np') rr=1.90d0
       if(trim(symbol) == 'Pu') rr=1.87d0
       if(trim(symbol) == 'Am') rr=1.80d0
       if(trim(symbol) == 'Cm') rr=1.69d0
       covlaentrr=rr
       end
!234567890
