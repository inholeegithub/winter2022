!234567890 
!      ifort -o xsf2poscar.x  strings.f90 xsf2poscar.f90
!      Written by In-Ho Lee, KRISS, March 27, 2016
       program xsf2poscar 
       USE strings, ONLY : parse,value,is_digit
       implicit none
       integer nspecies
       real*8 a1(3),a2(3),a3(3)
       integer, allocatable :: natoms(:)
       real*8, allocatable :: dir(:,:),car(:,:)
       character*2, allocatable :: symbl(:)
       character*200 char1
       character*2 jsym(1000000)
       integer ival(1000000)
       real*8 xyz(1000000,3)
       real*8 rnum
       integer i,j,k,ktotal,i0
       integer ilen,ipos
       integer ios,nargs
       character*200 str1
       character*200 args(40)
       character*20 delims
       logical lnew
       integer ispecies
       logical lfault
       character*2, external :: sy

       lfault=.false.
       ival=0
       do 
       read(5,'(a200)',err=911,end=999) str1
       delims=' '
       call parse(str1,delims,args,nargs)
       if(nargs > 0)then
       char1=trim(adjustl(args(1))) ; if(char1(1:2) == '##') args(1)='#'
       if(trim(args(1)) == '#') cycle
       if(args(1) == 'CRYSTAL') cycle
       if(args(1) == 'ATOMS')then
       write(6,*) 'ATOMS are not supported here.'
                             stop
                             endif
       if(args(1) == 'SLAB')then
       write(6,*) 'SLAB is not supported here.'
                            stop
                            endif
       if(args(1) == 'ANIMSTEPS')then
       write(6,*) 'ANIMSTEPS are not supported here.'
                                 stop
                                 endif
       if(args(1) == 'PRIMVEC')then
       read(5,*) a1(1),a1(2),a1(3)
       read(5,*) a2(1),a2(2),a2(3)
       read(5,*) a3(1),a3(2),a3(3)
                               endif
       if(args(1) == 'CONVVEC')then
       read(5,*) 
       read(5,*) 
       read(5,*) 
                               endif
       if(args(1) == 'PRIMCOORD')then
       read(5,*) nspecies
       k=0
       do 
!      read(5,*,err=911,end=999) jsym(k),xyz(k,1),xyz(k,2),xyz(k,3)
       read(5,'(a200)',err=911,end=999) str1
       delims=' '
       call parse(str1,delims,args,nargs)
       if(nargs > 1)then
       char1=trim(adjustl(args(1))) ; if(char1(1:2) == '##') args(1)='#'
       if(trim(args(1)) == '#') cycle
       k=k+1
       ios=0
       call value(args(1),ival(k),ios)
       if(ios == 0) jsym(k)=sy(ival(k))
       if(ios /= 0) jsym(k)=trim(adjustl(trim(args(1))))
       call value(args(2),xyz(k,1),ios)
       call value(args(3),xyz(k,2),ios)
       call value(args(4),xyz(k,3),ios)
!      write(6,*) jsym(k),xyz(k,1),xyz(k,2),xyz(k,3)
                    endif
       enddo
                                 endif
                    endif
       enddo
  911  continue
       lfault=.true.
  999  continue
       ktotal=k
       allocate(dir(k,3)) ; allocate(car(k,3))
       allocate(natoms(nspecies)) ; allocate(symbl(nspecies))
       ispecies=1 ; symbl(ispecies)=trim(adjustl(jsym(1)))
       if(ispecies == nspecies) goto 111
       do 
       do j=1,ktotal
       if(ispecies < nspecies)then
       lnew=.true.
       do i0=1,ispecies
       if(trim(adjustl(jsym(j))) == trim(adjustl(symbl(i0)))) lnew=.false.
       enddo
       if(lnew)then
       ispecies=ispecies+1
       symbl(ispecies)=trim(adjustl(jsym(j)))
       if(ispecies == nspecies) goto 111
               endif
                              endif
       enddo
       enddo
  111  continue
       k=0
       do ispecies=1,nspecies
       natoms(ispecies)=0
       do j=1,ktotal
       if(trim(adjustl(symbl(ispecies))) == trim(adjustl(jsym(j))))then
       natoms(ispecies)=natoms(ispecies)+1
       k=k+1
       car(k,:)=xyz(j,:)
                                                                   endif
       enddo
       enddo
       dir=car
       call tolaty(dir,a1,a2,a3,ktotal)
       write(6,*) 'head'
       write(6,*) '1.0'
       write(6,'(3f22.12)') a1(1),a1(2),a1(3)
       write(6,'(3f22.12)') a2(1),a2(2),a2(3)
       write(6,'(3f22.12)') a3(1),a3(2),a3(3)
       write(6,'(10(1x,a2,1x))') (symbl(i),i=1,nspecies)
       write(6,'(10i4)') (natoms(i),i=1,nspecies)
       k=1
       if(k > 0)then
       write(6,*) 'direct'
       do k=1,ktotal
       write(6,'(3f22.12)') dir(k,1),dir(k,2),dir(k,3)
       enddo
                else
       write(6,*) 'cartesian'
       do k=1,ktotal
       write(6,'(3f22.12)') car(k,1),car(k,2),car(k,3)
       enddo
                endif
       deallocate(symbl,natoms)
       deallocate(dir,car)
       stop
       end program xsf2poscar

!234567890
!      Written by In-Ho Lee, KRISS, March 27, 2016
       subroutine tolaty(xyz,a1,a2,a3,ktot)
       implicit none
       integer ktot
       real*8 xyz(ktot,3),a1(3),a2(3),a3(3)
       real*8 b(3,3),devid
       integer j,i
       real*8 d1,d2,d3

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
!      Written by In-Ho Lee, KRISS, March 27, 2016
       character*2 function sy(inti)
       implicit none
       integer inti

       if(inti ==  1) sy='H'
       if(inti ==  2) sy='He'
       if(inti ==  3) sy='Li'
       if(inti ==  4) sy='Be'
       if(inti ==  5) sy='B'
       if(inti ==  6) sy='C'
       if(inti ==  7) sy='N'
       if(inti ==  8) sy='O'
       if(inti ==  9) sy='F'
       if(inti == 10) sy='Ne'
       if(inti == 11) sy='Na'
       if(inti == 12) sy='Mg'
       if(inti == 13) sy='Al'
       if(inti == 14) sy='Si'
       if(inti == 15) sy='P'
       if(inti == 16) sy='S'
       if(inti == 17) sy='Cl'
       if(inti == 18) sy='Ar'
       if(inti == 19) sy='K'
       if(inti == 20) sy='Ca'
       if(inti == 21) sy='Sc'
       if(inti == 22) sy='Ti'
       if(inti == 23) sy='V'
       if(inti == 24) sy='Cr'
       if(inti == 25) sy='Mn'
       if(inti == 26) sy='Fe'
       if(inti == 27) sy='Co'
       if(inti == 28) sy='Ni'
       if(inti == 29) sy='Cu'
       if(inti == 30) sy='Zn'
       if(inti == 31) sy='Ga'
       if(inti == 32) sy='Ge'
       if(inti == 33) sy='As'
       if(inti == 34) sy='Se'
       if(inti == 35) sy='Br'
       if(inti == 36) sy='Kr'
       if(inti == 37) sy='Rb'
       if(inti == 38) sy='Sr'
       if(inti == 39) sy='Y'
       if(inti == 40) sy='Zr'
       if(inti == 41) sy='Nb'
       if(inti == 42) sy='Mo'
       if(inti == 43) sy='Tc'
       if(inti == 44) sy='Ru'
       if(inti == 45) sy='Rh'
       if(inti == 46) sy='Pd'
       if(inti == 47) sy='Ag'
       if(inti == 48) sy='Cd'
       if(inti == 49) sy='In'
       if(inti == 50) sy='Sn'
       if(inti == 51) sy='Sb'
       if(inti == 52) sy='Te'
       if(inti == 53) sy='I'
       if(inti == 54) sy='Xe'
       if(inti == 55) sy='Cs'
       if(inti == 56) sy='Ba'
       if(inti == 57) sy='La'
       if(inti == 58) sy='Ce'
       if(inti == 59) sy='Pr'
       if(inti == 60) sy='Nd'
       if(inti == 61) sy='Pm'
       if(inti == 62) sy='Sm'
       if(inti == 63) sy='Eu'
       if(inti == 64) sy='Gd'
       if(inti == 65) sy='Tb'
       if(inti == 66) sy='Dy'
       if(inti == 67) sy='Ho'
       if(inti == 68) sy='Er'
       if(inti == 69) sy='Tm'
       if(inti == 70) sy='Yb'
       if(inti == 71) sy='Lu'
       if(inti == 72) sy='Hf'
       if(inti == 73) sy='Ta'
       if(inti == 74) sy='W'
       if(inti == 75) sy='Re'
       if(inti == 76) sy='Os'
       if(inti == 77) sy='Ir'
       if(inti == 78) sy='Pt'
       if(inti == 79) sy='Au'
       if(inti == 80) sy='Hg'
       if(inti == 81) sy='Tl'
       if(inti == 82) sy='Pb'
       if(inti == 83) sy='Bi'
       if(inti == 84) sy='Po'
       if(inti == 85) sy='At'
       if(inti == 86) sy='Rn'
       if(inti == 87) sy='Fr'
       if(inti == 88) sy='Ra'
       if(inti == 89) sy='Ac'
       if(inti == 90) sy='Th'
       if(inti == 91) sy='Pa'
       if(inti == 92) sy='U'
       if(inti == 93) sy='Np'
       if(inti == 94) sy='Pu'
       if(inti == 95) sy='Am'
       if(inti == 96) sy='Cm'
       if(inti == 97) sy='Bk'
       if(inti == 98) sy='Cf'
       if(inti == 99) sy='Es'
       if(inti ==100) sy='Fm'
       if(inti ==101) sy='Md'
       if(inti ==102) sy='No'
       if(inti ==103) sy='Lr'
       if(inti ==104) sy='Rf'
       if(inti ==105) sy='Db'
       if(inti ==106) sy='Sg'
       if(inti ==107) sy='Bh'
       if(inti ==108) sy='Hs'
       if(inti ==109) sy='Mt'
       if(inti ==110) sy='Ds'
       if(inti ==111) sy='Rg'
       if(inti ==112) sy='Cn'
       return
       end 
