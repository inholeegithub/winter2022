!      Usage ifort -o calfinger.x calfinger.f90
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
!      Written by In-Ho Lee, KRISS, April 25, 2004.
       subroutine arb_grid(cdivid,npoint,x_start,x_finish,array_g)
       implicit none
       character(3) cdivid
       integer npoint
       real*8 x_start,x_finish,array_g(npoint)
       real*8 delta
       integer i

       if(npoint == 1)then
       array_g(1)=x_start
       return
       endif
       if(cdivid == 'lin' .or. cdivid =='LIN')then
       delta=(x_finish-x_start)/float(npoint-1)
       do i=1,npoint
       array_g(i)=x_start+delta*float(i-1)
       enddo
       elseif(cdivid == 'inv' .or. cdivid == 'INV')then
       delta=(1.d0/x_finish-1.d0/x_start)/float(npoint-1)
       do i=1,npoint
       array_g(i)=1.d0/(1.d0/x_start+delta*float(i-1))
       enddo
       elseif(cdivid == 'log' .or. cdivid == 'LOG')then
       delta=dlog(x_finish/x_start)/float(npoint-1)
       do i=1,npoint
       array_g(i)=exp(dlog(x_start)+delta*float(i-1))
       enddo
       else
       write(6,*) 'input error, check cdivid'
       stop
       endif
       end
!234567890
!      Written by In-Ho Lee, KRISS, September 21, 2020.
       real*8 function m33det(a)
       implicit none
       real*8, dimension(3,3), intent(in)  :: a

       m33det =  a(1,1)*a(2,2)*a(3,3)  &
               - a(1,1)*a(2,3)*a(3,2)  &
               - a(1,2)*a(2,1)*a(3,3)  &
               + a(1,2)*a(2,3)*a(3,1)  &
               + a(1,3)*a(2,1)*a(3,2)  &
               - a(1,3)*a(2,2)*a(3,1)
       end
!234567890
!      Written by In-Ho Lee, KRISS, September 21, 2020.
       subroutine tolaty0(xyz,a1,a2,a3,ktot)
       implicit none
       integer ktot
       real*8 xyz(ktot,3),a1(3),a2(3),a3(3)
       integer j
       real*8 b(3,3),devid,d1,d2,d3

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
       xyz(j,1)=d1 ; xyz(j,2)=d2 ; xyz(j,3)=d3
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
!      Written by In-Ho Lee, KRISS, September 21, 2020.
       subroutine getinv(natoms,nsigma,ninv,a1,a2,a3,dir,sktinv,sref,nrref,rref)
       implicit none
       integer natoms,nsigma,ninv,nrref,n1,n2,n3
       real*8 sktinv(ninv,nsigma,nrref,natoms),sref(nsigma),rref(nrref)
       real*8 a1(3),a2(3),a3(3),dir(natoms,3)
       real*8 skinv,vkinv,tk(3,3),amat(3,3),sigma
!      real*8 tkinv1,tkinv2,tkinv3
       real*8 ainv1,ainv2,ainv3,ainv4,ainv5,ainv6,ainv7
       real*8 xi,yi,zi,xj,yj,zj,trace1,trace2
       real*8 ck,pi,arg,brg,rij,x,y,z,tmp,deno
       integer iatom,jatom,ii,jj,kk,ksigma,iref
       real*8, external :: m33det

       tmp=2.d0*rref(nrref)
       n1=int(tmp/sqrt(dot_product(a1,a1)))+1
       n2=int(tmp/sqrt(dot_product(a2,a2)))+1
       n3=int(tmp/sqrt(dot_product(a3,a3)))+1
       sktinv(:,:,:,:)=0.d0
       pi=4.d0*atan(1.d0)
       do iref=1,nrref
       do ksigma=1,nsigma
       sigma=sref(ksigma)
       ck=1.d0/(2.d0*pi)**(0.5d0)/sigma
       do iatom=1,natoms
       xi=dir(iatom,1)*a1(1)+dir(iatom,2)*a2(1)+dir(iatom,3)*a3(1)
       yi=dir(iatom,1)*a1(2)+dir(iatom,2)*a2(2)+dir(iatom,3)*a3(2)
       zi=dir(iatom,1)*a1(3)+dir(iatom,2)*a2(3)+dir(iatom,3)*a3(3)
       deno=0.d0
       skinv=0.d0 ; vkinv=0.d0 ;  tk(:,:)=0.d0
       do jatom=1,natoms
       do ii=-n1,n1
       do jj=-n2,n2
       do kk=-n3,n3
       xj=dir(jatom,1)*a1(1)+dir(jatom,2)*a2(1)+dir(jatom,3)*a3(1)
       yj=dir(jatom,1)*a1(2)+dir(jatom,2)*a2(2)+dir(jatom,3)*a3(2)
       zj=dir(jatom,1)*a1(3)+dir(jatom,2)*a2(3)+dir(jatom,3)*a3(3)
       xj=xj+ii*a1(1)+jj*a2(1)+kk*a3(1)
       yj=yj+ii*a1(2)+jj*a2(2)+kk*a3(2)
       zj=zj+ii*a1(3)+jj*a2(3)+kk*a3(3)
       x=xi-xj ; y=yi-yj ; z=zi-zj 
       rij=sqrt(x*x+y*y+z*z)
       if(rij < 1.d-3) cycle
       arg=-((rij-rref(iref))/sigma)**2/2.d0 ; if(arg < -30.d0) arg=-30.d0
       brg=ck*exp(arg) 
       x=x/rij
       y=y/rij
       z=z/rij
       skinv=skinv+brg
       vkinv=vkinv+sqrt(x*x+y*y+z*z)*brg
       tk(1,1)=tk(1,1)+brg*x*x ; tk(1,2)=tk(1,2)+brg*x*y
       tk(1,3)=tk(1,3)+brg*x*z ; tk(2,1)=tk(2,1)+brg*y*x
       tk(2,2)=tk(2,2)+brg*y*y ; tk(2,3)=tk(2,3)+brg*y*z
       tk(3,1)=tk(3,1)+brg*z*x ; tk(3,2)=tk(3,2)+brg*z*y
       tk(3,3)=tk(3,3)+brg*z*z
       deno=deno+1.d0
       enddo
       enddo
       enddo
       enddo
       skinv=skinv/deno
       vkinv=vkinv/deno
       tk=tk/deno
       amat=matmul(tk,tk)
       trace1=tk(1,1)+tk(2,2)+tk(3,3)
       trace2=amat(1,1)+amat(2,2)+amat(3,3)
       ainv1=trace1 ; ainv2=(trace1**2-trace2)/2.d0 ; ainv2=max(ainv2,0.d0)
       ainv3=m33det(tk) ; ainv4=ainv1**2-2.d0*ainv2 ; ainv4=max(ainv4,0.d0)
       ainv5=ainv1**3-3.d0*ainv1*ainv2+3.d0*ainv3
       ainv6=2.d0/9.d0*(ainv1**2-3.d0*ainv2) ; ainv6=max(ainv6,0.d0)
       ainv7=2.d0/54.d0*(-9.d0*ainv1*ainv2+27.d0*ainv3+2.d0*ainv1**3)
       if(ninv == 9)then
       sktinv(1,ksigma,iref,iatom)=sktinv(1,ksigma,iref,iatom)+skinv
       sktinv(2,ksigma,iref,iatom)=sktinv(2,ksigma,iref,iatom)+vkinv
       sktinv(3,ksigma,iref,iatom)=sktinv(3,ksigma,iref,iatom)+ainv1
       sktinv(4,ksigma,iref,iatom)=sktinv(4,ksigma,iref,iatom)+sqrt(ainv2)
       sktinv(5,ksigma,iref,iatom)=sktinv(5,ksigma,iref,iatom)+sqrt(ainv4)
       sktinv(6,ksigma,iref,iatom)=sktinv(6,ksigma,iref,iatom)+sqrt(ainv6)
       sktinv(7,ksigma,iref,iatom)=sktinv(7,ksigma,iref,iatom)+(ainv3)**(1.d0/3.d0)
       sktinv(8,ksigma,iref,iatom)=sktinv(8,ksigma,iref,iatom)+(ainv5)**(1.d0/3.d0)
       sktinv(9,ksigma,iref,iatom)=sktinv(9,ksigma,iref,iatom)+(ainv7)**(1.d0/3.d0)
                    endif
       enddo
       enddo
       enddo
       end
!234567890
!      Written by In-Ho Lee, KRISS, September 21, 2020.
       subroutine invfeature(fname,gname,ninv0,nsigma0,nrref0)
       implicit none
       character*280 fname,gname
       integer ninv0,nsigma0,nrref0
       integer ninv,nsigma,nrref
       integer nspecies,natoms
       real*8 scale0,a1(3),a2(3),a3(3),vtest,tmp
       integer, allocatable :: nelements(:)
       character*2, allocatable :: symbl(:)
       real*8, allocatable :: dir(:,:)
       integer i,j,k,i0,j0,iref
       character*1 ch1
       real*8, allocatable :: sktinv(:,:,:,:)
       integer, parameter   :: nlen=1000
       integer              :: nwords,ipos
       character (len=nlen) :: text
       real*8 x_start,x_finish
       CHARACTER(3) cdivid
       real*8, allocatable :: sref(:),rref(:),work(:,:,:),vork(:)
       real*8, external :: atomicmass

       open(1,file=trim(fname),form='formatted')
       read(1,*)
       read(1,*)
       read(1,*)
       read(1,*)
       read(1,*)
       read(1,'(a1000)') text
       close(1)
       ipos = 1 ; nwords = 0
       loop: do
       i = verify(text(ipos:), ' ') !-- Find next non-blank.
       if (i == 0) exit loop        !-- No word found.
       nwords = nwords + 1          !-- Found something.
       ipos = ipos + i - 1          !-- Move to start of the word.
       i = scan(text(ipos:), ' ')   !-- Find next blank.
       if (i == 0) exit loop        !-- No blank found.
       ipos = ipos + i - 1          !-- Move to the blank.
       end do loop
       nspecies=nwords
       allocate(nelements(nspecies)) ; allocate(symbl(nspecies))
       open(1,file=trim(fname),form='formatted')
       read(1,*)
       read(1,*) scale0
       read(1,*) a1(1),a1(2),a1(3)
       read(1,*) a2(1),a2(2),a2(3)
       read(1,*) a3(1),a3(2),a3(3)
       if(scale0 > 0.d0)then
       a1=a1*scale0 ; a2=a2*scale0 ; a3=a3*scale0
                        endif
       if(scale0 < 0.d0)then
       vtest=(a1(2)*a2(3)-a1(3)*a2(2))*a3(1) &
            +(a1(3)*a2(1)-a1(1)*a2(3))*a3(2) &
            +(a1(1)*a2(2)-a1(2)*a2(1))*a3(3)
       vtest=abs(vtest)
       vtest=abs(scale0)/vtest ; vtest=vtest**(1.d0/3.d0)
       a1=a1*vtest ; a2=a2*vtest ; a3=a3*vtest
                        endif
       read(1,*) (symbl(i),i=1,nspecies)
       do i=1,nspecies
       symbl(i)=adjustl(symbl(i))
       enddo
       read(1,*) (nelements(i),i=1,nspecies)
       natoms=sum(nelements)
       allocate(dir(natoms,3))
  100  continue
       read(1,*) ch1
       if(ch1 == 'D') ch1='d'
       if(ch1 == 'C') ch1='c'
       if(ch1 == 'K') ch1='c'
       if(ch1 == 'k') ch1='c'
       if(ch1 == 'S') ch1='s'
       if(ch1 == 's') goto 100
       do i=1,natoms
       read(1,*) dir(i,1),dir(i,2),dir(i,3)
       enddo
       close(1)
!
       vtest=(a1(2)*a2(3)-a1(3)*a2(2))*a3(1) &
            +(a1(3)*a2(1)-a1(1)*a2(3))*a3(2) &
            +(a1(1)*a2(2)-a1(2)*a2(1))*a3(3)
       vtest=abs(vtest)
       tmp=0.d0
       do i=1,nspecies
       tmp=tmp+atomicmass(symbl(i))*nelements(i)
       enddo
       tmp=tmp/vtest*(1.660539d-24)/(1.d-24)
       write(6,*) tmp, ' g/cm^3'
!
       if(ch1 == 'c')then
       dir=dir*scale0
       call tolaty0(dir,a1,a2,a3,natoms)
                     endif
!      cdivid='inv'
!      cdivid='log'
       cdivid='lin'
       nrref=nrref0
       allocate(rref(nrref))
       x_start=0.00d0 ; x_finish=8.00d0
       x_start=0.50d0 ; x_finish=8.00d0
       call arb_grid(cdivid,nrref,x_start,x_finish,rref)
!      do i=1,nrref
!      write(6,*) rref(i)
!      enddo
!      cdivid='log'
!      cdivid='inv'
       cdivid='lin'
       nsigma=nsigma0
       allocate(sref(nsigma))
!      x_start=0.200d0  ; x_finish=1.20d0
!      x_start=0.050d0  ; x_finish=0.50d0
       x_start=0.100d0  ; x_finish=0.50d0
       call arb_grid(cdivid,nsigma,x_start,x_finish,sref)
!      do i=1,nsigma
!      write(6,*) sref(i)
!      enddo
       ninv=9
       ninv=ninv0
       allocate(sktinv(ninv,nsigma,nrref,natoms))
       call getinv(natoms,nsigma,ninv,a1,a2,a3,dir,sktinv,sref,nrref,rref)
       if(.true.)then
       allocate(work(ninv,nsigma,nrref))
       open(18,file='./fingers/'//trim(gname),form='formatted')
!      open(17,file='./fingers/'//trim(gname),form='formatted')
!      write(17,'(a1,2x,3i5)') '#',ninv,nsigma,nrref
!      open(7,file='fort.7',form='formatted')
       work(:,:,:)=0.d0
       do iref=1,nrref
       do j=1,nsigma
       do k=1,ninv
       i0=0
       do i=1,nspecies
       do j0=1,nelements(i)
       i0=i0+1
       work(k,j,iref)=work(k,j,iref)+sktinv(k,j,iref,i0)
       enddo
       enddo
       enddo
       enddo
       enddo
       work=work/dble(i0)
       allocate(vork(ninv*nsigma*nrref))
       j0=0
       do k=1,ninv
       do j=1,nsigma
!      write(7,'(a1,2x,2i5)') '#', k,j
       do iref=1,nrref
!!     write(7,'(i8,e20.8)') iref,work(k,j,iref)
!      write(7,'(e18.8,1x,e20.8)') rref(iref),work(k,j,iref)
       write(18,'(e18.8,1x,e20.8)') rref(iref),work(k,j,iref)
       j0=j0+1
       vork(j0)=work(k,j,iref)
!      write(17,*) vork(j0)
       enddo
!      write(7,*) '&'
       write(18,*) '&'
       enddo
       enddo
       close(18)
!      close(7)
       close(17)
       deallocate(work)
       deallocate(vork)
                  endif
       deallocate(symbl) ; deallocate(nelements) ; deallocate(dir)
       deallocate(sktinv)
       deallocate(sref)
       deallocate(rref)
       end
!234567890
!      Written by In-Ho Lee, KRISS, October 2, 2020.
       program calfinger
       implicit none
       integer i,j
       real*8 tmp
       character*280 fname,gname
       integer ninv0,nsigma0,nrref0
       integer ntargets
       logical lexist

       ninv0=9 ; nsigma0=1 ; nrref0=200
       fname='POSCAR_tmq' ; gname='POSCAR_tmq_fin'
       if(.true.)then
       call system('mkdir fingers')
       call system('sleep 0.1')
       open(19,file='targets',form='formatted')
       read(19,*) ntargets
       ntargets=iabs(ntargets)
       if( ntargets > 0)then
       write(6,*) ntargets
       do i=1,ntargets
       read(19,'(a280)') fname
       fname=adjustl(fname) ; fname=trim(fname)
       gname=trim(fname)//'_finplt' ; gname=trim(gname)
       write(6,*) trim(fname),' ', trim(gname)
       call invfeature(fname,gname,ninv0,nsigma0,nrref0)
       enddo
                        endif
       close(19)
                 endif
       inquire(file='fort.7',exist=lexist)
!      if(lexist)then
!      open(7,file='fort.7',form='formatted')
!      close(7,status='delete')
!                endif
       end program calfinger
!234567890
