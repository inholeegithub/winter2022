!      ifort toposcar.f90 numeral.f cmds.f 
!      This is useful especially for enthalpy optimization case.
!      ./a.out < csa.in 
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       program read_fort1
       implicit none
       integer ndeg,npop,npop1,nspecies
       real*8 davg,energy_best
       real*8, allocatable :: posi(:,:),posi1(:,:),energy_sorted(:),energy_sorted1(:),posi_best(:)
       real*8, allocatable :: dir(:,:),exdir(:,:)
       integer i,j,na,ish,ipop
       integer n1,n2,n3
       real*8 rcut
       real*8 r6(6),cmat(3,3),a1(3),a2(3),a3(3),celvol,tmp
       character*2, allocatable :: symbl(:)
       integer, allocatable :: nelements(:)
       integer isize
       character*280 fname,gname
       real*8, external :: atomicmass

       read(5,*) nspecies
       allocate(symbl(nspecies),nelements(nspecies))
       read(5,*)(symbl(i),i=1,nspecies)
       do i=1,nspecies
       symbl(i)=adjustl(symbl(i))
       enddo
       read(5,*)(nelements(i),i=1,nspecies)
       open(1,file='fort.1',form='formatted')
       read(1,*) ndeg,npop,npop1
       allocate(posi(ndeg,npop)) ; allocate(posi1(ndeg,npop1))
       allocate(energy_sorted(npop)) ; allocate(energy_sorted1(npop1))
       allocate(posi_best(ndeg))
       read(1,*) posi1,energy_sorted1
       read(1,*) posi,energy_sorted
       read(1,*) posi_best,energy_best,davg
       close(1)
       rcut=6.0d0
       n1=3
       n2=3
       n3=3
       na=0
       do i=1,nspecies
       do j=1,nelements(i)
       na=na+1
       enddo
       enddo
       allocate(dir(na,3))
       allocate(exdir(n1*n2*n3*na,3))

       ish=ndeg-6
       do ipop=1,npop
       do i=1,6
       r6(i)=posi(ish+i,ipop)
       enddo
       call latmat(r6,cmat,1)
       a1(:)=cmat(1,:) ; a2(:)=cmat(2,:) ; a3(:)=cmat(3,:)
       celvol=       (cmat(1,2)*cmat(2,3)-cmat(1,3)*cmat(2,2))*cmat(3,1)
       celvol=celvol+(cmat(1,3)*cmat(2,1)-cmat(1,1)*cmat(2,3))*cmat(3,2)
       celvol=celvol+(cmat(1,1)*cmat(2,2)-cmat(1,2)*cmat(2,1))*cmat(3,3)
       celvol=abs(celvol)
!
       tmp=0.d0
       do i=1,nspecies
       tmp=tmp+atomicmass(symbl(i))*nelements(i)
       enddo
       tmp=tmp/celvol*(1.660539d-24)/(1.d-24)
!
       isize=4 ; call xnumeral(ipop,fname,isize)
       fname='POSCAR_'//trim(fname)
       isize=4 ; call xnumeral(ipop,gname,isize)
       gname='QOSCAR_'//trim(gname)
       write(6,*) trim(fname)
       open(71,file=trim(fname),form='formatted')
       write(71,'(2e23.15,1x,f20.4)') energy_sorted(ipop),celvol,tmp
       write(71,'(a3)') '1.0'
       write(71,'(3f23.16)') a1(1),a1(2),a1(3)
       write(71,'(3f23.16)') a2(1),a2(2),a2(3)
       write(71,'(3f23.16)') a3(1),a3(2),a3(3)
       write(71,'(20(2x,a2,1x))') (symbl(i),i=1,nspecies)
       write(71,'(20(i4,1x))') (nelements(i),i=1,nspecies)
       write(71,'(a6)') "Direct"
       na=0
       do i=1,nspecies
       do j=1,nelements(i)
       na=na+1
       write(71,'(3f20.16)') posi(3*(na-1)+1,ipop),posi(3*(na-1)+2,ipop),posi(3*(na-1)+3,ipop)
       dir(na,1)=posi(3*(na-1)+1,ipop)
       dir(na,2)=posi(3*(na-1)+2,ipop)
       dir(na,3)=posi(3*(na-1)+3,ipop)
       enddo
       enddo
       close(71)
       call expposcar(n1,n2,n3,nspecies,symbl,nelements,na,a1,a2,a3,dir,exdir,gname)
       enddo
       write(6,*)
       do ipop=1,npop
       do i=1,6
       r6(i)=posi(ish+i,ipop)
       enddo
       call latmat(r6,cmat,1)
       a1(:)=cmat(1,:) ; a2(:)=cmat(2,:) ; a3(:)=cmat(3,:)
       celvol=       (cmat(1,2)*cmat(2,3)-cmat(1,3)*cmat(2,2))*cmat(3,1)
       celvol=celvol+(cmat(1,3)*cmat(2,1)-cmat(1,1)*cmat(2,3))*cmat(3,2)
       celvol=celvol+(cmat(1,1)*cmat(2,2)-cmat(1,2)*cmat(2,1))*cmat(3,3)
       celvol=abs(celvol)
       write(6,'(2e23.15)') celvol,energy_sorted(ipop)
       enddo
       write(6,*)
       do ipop=1,npop
       do i=1,6
       r6(i)=posi(ish+i,ipop)
       enddo
       call latmat(r6,cmat,1)
       a1(:)=cmat(1,:) ; a2(:)=cmat(2,:) ; a3(:)=cmat(3,:)
       celvol=       (cmat(1,2)*cmat(2,3)-cmat(1,3)*cmat(2,2))*cmat(3,1)
       celvol=celvol+(cmat(1,3)*cmat(2,1)-cmat(1,1)*cmat(2,3))*cmat(3,2)
       celvol=celvol+(cmat(1,1)*cmat(2,2)-cmat(1,2)*cmat(2,1))*cmat(3,3)
       celvol=abs(celvol)
       write(6,'(2e23.15)') celvol/float(na),energy_sorted(ipop)/float(na)
       enddo
       write(6,*)
       tmp=minval(energy_sorted(:))
       write(6,*) tmp
       tmp=energy_sorted(1)
       write(6,*) tmp
       do ipop=1,npop
       do i=1,6
       r6(i)=posi(ish+i,ipop)
       enddo
       call latmat(r6,cmat,1)
       a1(:)=cmat(1,:) ; a2(:)=cmat(2,:) ; a3(:)=cmat(3,:)
       celvol=       (cmat(1,2)*cmat(2,3)-cmat(1,3)*cmat(2,2))*cmat(3,1)
       celvol=celvol+(cmat(1,3)*cmat(2,1)-cmat(1,1)*cmat(2,3))*cmat(3,2)
       celvol=celvol+(cmat(1,1)*cmat(2,2)-cmat(1,2)*cmat(2,1))*cmat(3,3)
       celvol=abs(celvol)
       write(6,'(2e23.15)') celvol/float(na),(energy_sorted(ipop)-tmp)/float(na)
       enddo
       deallocate(dir)
       deallocate(exdir)
       deallocate(symbl,nelements)
       deallocate(posi,posi1,energy_sorted,energy_sorted1)
       deallocate(posi_best)
       write(6,*) n1,n2,n3
       write(6,*) nspecies
       write(6,*) npop
       write(6,*) rcut
       call mdsdriver(n1,n2,n3,nspecies,npop,rcut)
       stop
       end program read_fort1
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine latmat(rlat,wmat,ksign)
       implicit none
       integer ksign
       real*8 rlat(6),wmat(3,3)
       real*8 ra,rb,rc,cosinea,cosineb,cosinec
       real*8 epslat,tmr
       integer i,j

       epslat=1.0d-6
       if(ksign == 1)then
       wmat=0.0d0
       wmat(1,1)=rlat(1)
       wmat(2,1)=rlat(2)*cos(rlat(6))
       wmat(2,2)=rlat(2)*sin(rlat(6))
       wmat(3,1)=rlat(3)*cos(rlat(5))
       wmat(3,2)=rlat(3)*cos(rlat(4))*sin(rlat(6)) &
        -((rlat(3)*cos(rlat(5))-rlat(3)*cos(rlat(4))*cos(rlat(6)))/tan(rlat(6)))
       tmr=rlat(3)**2-wmat(3,1)**2-wmat(3,2)**2  ; if(tmr <= 1.0d-12) tmr=0.0d0
       wmat(3,3)=sqrt(tmr)
       do i=1,3
       do j=1,3
       if(abs(wmat(i,j)) < epslat) wmat(i,j)=0.0d0
       enddo
       enddo
                     else
       rlat=0.0d0
       ra=sqrt(wmat(1,1)**2+wmat(1,2)**2+wmat(1,3)**2)
       rb=sqrt(wmat(2,1)**2+wmat(2,2)**2+wmat(2,3)**2)
       rc=sqrt(wmat(3,1)**2+wmat(3,2)**2+wmat(3,3)**2)
       cosinea=(wmat(2,1)*wmat(3,1)+wmat(2,2)*wmat(3,2)+wmat(2,3)*wmat(3,3))/rb/rc
       cosineb=(wmat(1,1)*wmat(3,1)+wmat(1,2)*wmat(3,2)+wmat(1,3)*wmat(3,3))/rc/ra
       cosinec=(wmat(1,1)*wmat(2,1)+wmat(1,2)*wmat(2,2)+wmat(1,3)*wmat(2,3))/ra/rb
       rlat(1)=ra ; rlat(2)=rb ; rlat(3)=rc
       rlat(4)=acos(cosinea) ; rlat(5)=acos(cosineb) ; rlat(6)=acos(cosinec)
                     endif
       return
       end
!234567890
       subroutine expposcar(n1,n2,n3,nspecies,symbl,nelements,na,a1,a2,a3,dir,exdir,gname)
       implicit none
       character*2 symbl(nspecies)
       character*280 gname
       integer n1,n2,n3,nspecies,na,nelements(nspecies)
       real*8 a1(3),a2(3),a3(3),dir(na,3),exdir(na*n1*n2*n3,3)
       integer i,j,k,ia,is
       real*8 cmatrix(3,3)

       is=0
       do ia=1,na
       do i=1,n1
       do j=1,n2
       do k=1,n3
       is=is+1
       exdir(is,1)=dir(ia,1)+dble(i-1)
       exdir(is,2)=dir(ia,2)+dble(j-1)
       exdir(is,3)=dir(ia,3)+dble(k-1)
       enddo
       enddo
       enddo
       enddo
       do ia=1,is
       exdir(ia,1)=exdir(ia,1)/dble(n1)
       exdir(ia,2)=exdir(ia,2)/dble(n2)
       exdir(ia,3)=exdir(ia,3)/dble(n3)
       enddo
       cmatrix(1,:)=n1*a1(:) ; cmatrix(2,:)=n2*a2(:) ; cmatrix(3,:)=n3*a3(:)
       open(22,file=trim(gname),form='formatted')
       write(22,*) n1,n2,n3
       write(22,*) '1.0'
       write(22,'(3f24.12)') cmatrix(1,1),cmatrix(1,2),cmatrix(1,3)
       write(22,'(3f24.12)') cmatrix(2,1),cmatrix(2,2),cmatrix(2,3)
       write(22,'(3f24.12)') cmatrix(3,1),cmatrix(3,2),cmatrix(3,3)
       write(22,'(20(2x,a2,1x))') (symbl(i),i=1,nspecies)
       write(22,'(20(i4,1x))') (n1*n2*n3*nelements(i),i=1,nspecies)
       write(22,*) 'direct'
       do ia=1,is
       write(22,'(3f22.12)') exdir(ia,1),exdir(ia,2),exdir(ia,3)
       enddo
       close(22)
       end
!234567890
       subroutine mdsdriver(n1,n2,n3,nspecies,npop,rcut)
       implicit none
       integer n1,n2,n3,npop,nspecies
       real*8 rcut
       integer nsamples,mdim
       integer i,i2,j,k,ierr,iprint
       integer na,ipop,isize
       real*8 cmat(3,3),d1,d2,d3,x,y,z,r
       character*2, allocatable :: symbl1(:)
       integer, allocatable :: nelements1(:)
       real*8, allocatable :: bdldist(:,:),exdir(:,:)
       character*280 gname
       real*4 valm
       real*4, allocatable :: dissim(:,:),aa(:,:),work1(:),work2(:),a2work(:,:)

       mdim=int(rcut/0.1d0)+1
       allocate(bdldist(mdim,npop)) ; bdldist=0.d0
       allocate(symbl1(nspecies)) ; allocate(nelements1(nspecies))
       do ipop=1,npop
       isize=4 ; call xnumeral(ipop,gname,isize)
       gname='QOSCAR_'//trim(gname)
       open(11,file=trim(gname),form='formatted')
       read(11,*)
       read(11,*)
       read(11,*) cmat(1,1),cmat(1,2),cmat(1,3)
       read(11,*) cmat(2,1),cmat(2,2),cmat(2,3)
       read(11,*) cmat(3,1),cmat(3,2),cmat(3,3)
       read(11,*) (symbl1(i),i=1,nspecies)
       read(11,*) (nelements1(i),i=1,nspecies)
       read(11,*)
       if(ipop == 1)then
       na=0
       do i=1,nspecies
       do j=1,nelements1(i)
       na=na+1
       enddo
       enddo
       allocate(exdir(na,3))
                    endif
       do i=1,na
       read(11,*) exdir(i,1),exdir(i,2),exdir(i,3)
       enddo
       close(11)
       do i=1,na-1
       do j=i+1,na
       d1=exdir(i,1)-exdir(j,1)
       d2=exdir(i,2)-exdir(j,2)
       d3=exdir(i,3)-exdir(j,3)
       d1=d1-anint(d1)
       d2=d2-anint(d2)
       d3=d3-anint(d3)
       x=d1*cmat(1,1)+d2*cmat(2,1)+d3*cmat(3,1)
       y=d1*cmat(1,2)+d2*cmat(2,2)+d3*cmat(3,2)
       z=d1*cmat(1,3)+d2*cmat(2,3)+d3*cmat(3,3)
       r=sqrt(x*x+y*y+z*z)
       if(r < rcut)then
       k=int(r/0.1d0)+1
       bdldist(k,ipop)=bdldist(k,ipop)+1.d0
                   endif
       enddo
       enddo
       enddo
       deallocate(exdir) ; deallocate(symbl1,nelements1)
!
       nsamples=npop
        write(6,*) 'original data',nsamples,mdim
        allocate(aa(nsamples,mdim))
        allocate(dissim(nsamples,nsamples))
        allocate(a2work(nsamples,nsamples))
        allocate(work1(nsamples),work2(nsamples))
        do i=1,nsamples
        do j=1,mdim
        aa(i,j)=bdldist(j,i)
        enddo
        enddo
       deallocate(bdldist)
!
!-----  Construct distances  ----------------------------------------
!
        valm = 0.0
        do 3 i = 1,nsamples
           do 2 i2 = 1,nsamples
              dissim(i,i2) = 0.0
                 do 1 j = 1,mdim
                    dissim(i,i2)=dissim(i,i2)+(aa(i,j)-aa(i2,j))**2
   1             continue
              dissim(i,i2)=sqrt(dissim(i,i2))
           if (dissim(i,i2).gt.valm) valm = dissim(i,i2)
   2       continue
   3    continue
!
!-----  Normalize to maximum distance = unity  ----------------------
!
        do 5 i = 1,nsamples
           do 4 j = 1,nsamples
              dissim(i,j) = dissim(i,j)/valm
   4       continue
   5    continue
!
        ierr = 0
        iprint = 2
        call cmds(nsamples,dissim,iprint,work1,work2,a2work,ierr)
        if (ierr.ne.0) goto 9000
!
        goto 9900
 9000   write (6,*) ' ABNORMAL END: IERR =', IERR
 9900   continue
       deallocate(aa,work1,work2,a2work,dissim)
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
