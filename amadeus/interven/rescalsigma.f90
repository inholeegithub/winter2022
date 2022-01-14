!    purpose:
!    to generate POSCAR file with a better set of atomic positions
!    this is only useful for the crystal structures with the serial number greater than npop.
!    Input : ../csa.in, POSCAR
!    Output: POSCAR
!    usage :
!    mkdir 0000
!    cp ~/csa_vasp/interven/rescalsigma.f90 .
!    ifort rescalsigma.f90
!    ./a.out
!    OR
!    cd 0001
!    cp ~/csa_vasp/interven/rescalsigma.f90 .
!    ifort rescalsigma.f90
!    ./a.out   
!    echo ./0???  | xargs -n 1 cp a.out
!    vi ./0???/CSA_SOLDIER.pbs
!    INCAR_rlx : use a small value for POTIM = 0.1 (or 0.05) in case of steric hindrance
!
!234567890
!      Written by In-Ho Lee, KRISS, October 6, 2017.
       implicit none
       integer nspecies,npop
       real*8 a1(3),a2(3),a3(3),scale0
       real*8 factor,amp
       real*8, allocatable :: dir(:,:),sigmamatrix(:,:)
       real*8, allocatable :: sigmamatrixref(:,:)
       integer, allocatable :: nelements(:)
       character*2, allocatable :: symbl(:)
       integer na,istructure
       integer i,j
       character*9 ch9_1,ch9_2
       character*1 ch1
       real*8 cmat(3,3),vtest,x,y,z,d1,d2,d3
       logical lfault2,lexist
       real*8, external :: covlaentrr

       lfault2=.false.
       inquire(file='../csa.in',exist=lexist)
       if(.not. lexist)then
       write(6,*) 'somehow ../csa.in is not present'
       lfault2=.true.
                       stop
                       endif
       open(80,file='../csa.in',form='formatted')
       read(80,*) nspecies
       allocate(sigmamatrix(nspecies,nspecies))
       allocate(sigmamatrixref(nspecies,nspecies))
       allocate(symbl(nspecies)) ; allocate(nelements(nspecies))
       read(80,*) (symbl(i),i=1,nspecies)
       write(6,'(10(1x,a2,1x))') (symbl(i),i=1,nspecies)
       write(6,'(10f6.3)') (covlaentrr(symbl(i)),i=1,nspecies)
       do i=1,nspecies
       do j=1,nspecies
       sigmamatrixref(i,j)=(covlaentrr(symbl(i))+covlaentrr(symbl(j)))/2.d0
       enddo
       enddo
       write(6,*) 'reference sigmamatrix from covalent radii table (Angstrom)'
       do i=1,nspecies
       write(6,*) (sigmamatrixref(i,j),j=1,nspecies)
       enddo
       do i=1,5
       read(80,*)
       enddo
       do i=1,nspecies
       read(80,*) (sigmamatrix(i,j),j=1,nspecies)
       enddo
!      do i=1,nspecies
!      write(6,*) (sigmamatrix(i,j),j=1,nspecies)
!      enddo
       do i=1,4
       read(80,*)
       enddo
       read(80,*) npop
       close(80)
       inquire(file='POSCAR',exist=lexist)
       if(.not. lexist)then
       write(6,*) 'somehow POSCAR is not present'
       write(6,*) 'otherwise, it is a test run on a test directory, e.g. 0000'
       write(6,*) 'you may consider the pressure value for the sigmamatrix setting'
       write(6,*) 'cut and paste'
       lfault2=.true.
       goto 909
                       stop
                       endif
       open(81,file='POSCAR',form='formatted')
       read(81,*,err=911,end=999) istructure
       read(81,*,err=911,end=999) scale0
       read(81,*,err=911,end=999) a1(1),a1(2),a1(3)
       read(81,*,err=911,end=999) a2(1),a2(2),a2(3)
       read(81,*,err=911,end=999) a3(1),a3(2),a3(3)
       cmat(1,:)=a1(:) ; cmat(2,:)=a2(:) ; cmat(3,:)=a3(:)
       if(scale0 > 0.d0)then
       cmat=cmat*scale0
                        endif
       if(scale0 < 0.d0)then
       vtest=(cmat(1,2)*cmat(2,3)-cmat(1,3)*cmat(2,2))*cmat(3,1) &
            +(cmat(1,3)*cmat(2,1)-cmat(1,1)*cmat(2,3))*cmat(3,2) &
            +(cmat(1,1)*cmat(2,2)-cmat(1,2)*cmat(2,1))*cmat(3,3)
       vtest=(abs(scale0)/abs(vtest))**(1.d0/3.d0)
       cmat=cmat*vtest
                        endif
       a1(:)=cmat(1,:) ; a2(:)=cmat(2,:) ; a3(:)=cmat(3,:)
       read(81,*,err=911,end=999) (symbl(j),j=1,nspecies)
       do i=1,nspecies
       symbl(i)=adjustl(symbl(i))
       enddo
       read(81,*,err=911,end=999) (nelements(j),j=1,nspecies) 
       na=sum(nelements)
       allocate(dir(na,3)) 
       read(81,*,err=911,end=999) ch9_1
       ch1=ch9_1(1:1)
       if(ch1 == 'S') ch9_1='Selective' ; if(ch1 == 's') ch9_1='Selective'
       if(ch1 == 'C') ch9_1='Cartesian' ; if(ch1 == 'c') ch9_1='Cartesian'
       if(ch1 == 'D') ch9_1='Direct'    ; if(ch1 == 'd') ch9_1='Direct'
       if(ch9_1 == 'Selective' .or. ch9_1 == 'selective')then
       read(81,*,err=911,end=999) ch9_2
       ch1=ch9_2(1:1)
       if(ch1 == 'C') ch9_2='Cartesian' ; if(ch1 == 'c') ch9_2='Cartesian'
       if(ch1 == 'D') ch9_2='Direct'    ; if(ch1 == 'd') ch9_2='Direct'
       if(ch9_2 == 'Cartesian' .or. ch9_2 == 'cartesian')then
       do i=1,na
       read(81,*,err=911,end=999) x,y,z
       if(scale0 > 0.d0)then
       x=x*scale0 ; y=y*scale0 ; z=z*scale0
                        endif
       dir(i,1)=x ; dir(i,2)=y ; dir(i,3)=z
       enddo
       call tolaty(na,a1,a2,a3,dir)
       goto 900
                                                         else
       do i=1,na
       read(81,*,err=911,end=999) d1,d2,d3
       dir(i,1)=d1 ; dir(i,2)=d2 ; dir(i,3)=d3
       enddo
       goto 900
                                                         endif
                                                         else
       if(ch9_1 == 'Cartesian' .or. ch9_1 == 'cartesian')then
       do i=1,na
       read(81,*,err=911,end=999) x,y,z
       if(scale0 > 0.d0)then
       x=x*scale0 ; y=y*scale0 ; z=z*scale0
                        endif
       dir(i,1)=x ; dir(i,2)=y ; dir(i,3)=z
       enddo
       call tolaty(na,a1,a2,a3,dir)
       goto 900
                                                         else
       do i=1,na
       read(81,*,err=911,end=999) d1,d2,d3
       dir(i,1)=d1 ; dir(i,2)=d2 ; dir(i,3)=d3
       enddo
       goto 900
                                                         endif
                                                         endif
  911  continue
  999  continue
       lfault2=.true.
  900  continue
       close(81)
       factor=1.5d0
       call check_relax(istructure,npop,nspecies,symbl,nelements,na,amp,factor,sigmamatrix,dir,a1,a2,a3)
       factor=2.0d0
       factor=0.0d0
       call check_relax(istructure,npop,nspecies,symbl,nelements,na,amp,factor,sigmamatrix,dir,a1,a2,a3)
  909  continue
       if(allocated(dir)) deallocate(dir)
       if(allocated(symbl)) deallocate(symbl)
       if(allocated(nelements)) deallocate(nelements)
       if(allocated(sigmamatrix)) deallocate(sigmamatrix)
       stop
       end
!234567890
!      Written by In-Ho Lee, KRISS, October 6, 2017.
       subroutine check_relax(istructure,npop,nspecies,symbl,nelements,na,amp,factor,sigmamatrix,dir,a1,a2,a3)
       implicit none
       integer istructure,npop,nspecies,na,nelements(nspecies)
       character*2 symbl(nspecies)
       real*8 amp,factor,sigmamatrix(nspecies,nspecies),dir(na,3)
       real*8 a1(3),a2(3),a3(3)
       real*8 d1,d2,d3,x,y,z,aainv,bbinv,ccinv,sclz,test,tmp,vec(3),wec(3)
       integer i,k,k0,j,iter,niter,kprint
       integer, allocatable :: klist(:),itype(:)
       real*8, allocatable :: sigmamatrix0(:,:),dir0(:,:)
       real*8, external :: covlaentrr

       kprint=0
       kprint=1
       allocate(dir0(na,3)) ; allocate(itype(na))
       k=0
       do i=1,nspecies
       do j=1,nelements(i)
       k=k+1 ; itype(k)=i
       enddo
       enddo
       k=0
       do i=1,na-1
       do j=i+1,na
       d1=dir(i,1)-dir(j,1)
       d2=dir(i,2)-dir(j,2)
       d3=dir(i,3)-dir(j,3)
       d1=d1-anint(d1) ; d2=d2-anint(d2) ; d3=d3-anint(d3)
       x=d1*a1(1)+d2*a2(1)+d3*a3(1)
       y=d1*a1(2)+d2*a2(2)+d3*a3(2)
       z=d1*a1(3)+d2*a2(3)+d3*a3(3)
       test=sqrt(x**2+y**2+z**2)
       if(test < sigmamatrix(itype(i),itype(j)))then
       k=k+1
       if(kprint > 0) write(6,'(f13.4,1x,a1)') test, 'A'
                                                endif
       enddo
       enddo
       if(kprint > 0)then
       write(6,*) k
       if(k > 0) write(6,*) 'something went wrong in AMADEUS'
       if(k == 0) write(6,*) 'at least distances are properly treated in AMADEUS with sigmamatrix'
       do i=1,nspecies
       write(6,*) (sigmamatrix(i,j),j=1,nspecies)
       enddo
                     endif
       if(istructure <= npop)then
       if(kprint > 0) write(6,*) 'we stop here ',istructure
       goto 991
                             endif
!      factor=1.5d0
!      factor=2.0d0
       amp=1.d-1
       if(kprint > 0)then
       write(6,*) 'however, we want to have a modified sigmamatrix with factor',factor
                     endif
       allocate(sigmamatrix0(nspecies,nspecies))
       if(factor <= 0.d0)then
       do i=1,nspecies
       do j=1,nspecies
       sigmamatrix0(i,j)=(covlaentrr(symbl(i))+covlaentrr(symbl(j)))/2.d0
       enddo
       enddo
       write(6,*) 'a new generation of sigmamatrix (Angstrom)'
!      do i=1,nspecies
!      write(6,*) (sigmamatrix0(i,j),j=1,nspecies)
!      enddo
                         else
       sigmamatrix0=sigmamatrix*factor
                         endif
       if(kprint > 0)then
       do i=1,nspecies
       write(6,*) (sigmamatrix0(i,j),j=1,nspecies)
       enddo
                     endif
       k=0
       do i=1,na-1
       do j=i+1,na
       d1=dir(i,1)-dir(j,1) ; d2=dir(i,2)-dir(j,2) ; d3=dir(i,3)-dir(j,3)
       d1=d1-anint(d1) ; d2=d2-anint(d2) ; d3=d3-anint(d3)
       x=d1*a1(1)+d2*a2(1)+d3*a3(1)
       y=d1*a1(2)+d2*a2(2)+d3*a3(2)
       z=d1*a1(3)+d2*a2(3)+d3*a3(3)
       test=sqrt(x**2+y**2+z**2)
       if(test < sigmamatrix0(itype(i),itype(j)))then
       k=k+1
       if(kprint > 0) write(6,'(f13.4,1x,a1)') test, 'A'
                                                 endif
       enddo
       enddo
       if(kprint > 0) write(6,*) k,' steric hindrance at starting point'
       k0=k ; dir0=dir
!
       sclz=10.d0
       aainv=sclz/sqrt(dot_product(a1,a1))
       bbinv=sclz/sqrt(dot_product(a2,a2))
       ccinv=sclz/sqrt(dot_product(a3,a3))
       allocate(klist(na))
       call random_seed()
       niter=10000000
       do iter=1,niter
       k=0 ; klist=0
       do i=1,na-1
       do j=i+1,na
       vec(:)=dir(i,:)-dir(j,:) ; vec=vec-anint(vec)
       wec(1)=vec(1)*a1(1)+vec(2)*a2(1)+vec(3)*a3(1)
       wec(2)=vec(1)*a1(2)+vec(2)*a2(2)+vec(3)*a3(2)
       wec(3)=vec(1)*a1(3)+vec(2)*a2(3)+vec(3)*a3(3)
       test=sqrt(dot_product(wec,wec))
       if(test < sigmamatrix0(itype(i),itype(j)))then
       k=k+1 ; klist(i)=klist(i)+1 ; klist(j)=klist(j)+1
                                                 endif
       enddo
       enddo
       if(k < k0)then
       k0=k ; dir0=dir
       if(kprint > 0) write(6,*) k,' steric hindrance in a regular cell, factor',factor
                 endif
       if(k > 0)then
       do i=1,na
       if(klist(i) > 0)then
       call random_number(tmp) 
       dir(i,1)=dir0(i,1)+(tmp-0.5d0)*amp*aainv
       call random_number(tmp) 
       dir(i,2)=dir0(i,2)+(tmp-0.5d0)*amp*bbinv
       call random_number(tmp) 
       dir(i,3)=dir0(i,3)+(tmp-0.5d0)*amp*ccinv
                       endif
       enddo
!-----{
       call random_number(tmp) 
       if(tmp < 0.03)then
       do i=1,na
       call random_number(tmp) 
       dir(i,1)=dir(i,1)+(tmp-0.5d0)*amp*aainv
       call random_number(tmp) 
       dir(i,2)=dir(i,2)+(tmp-0.5d0)*amp*bbinv
       call random_number(tmp) 
       dir(i,3)=dir(i,3)+(tmp-0.5d0)*amp*ccinv
       enddo
                     endif
!-----}
                else
       goto 200
                endif
       enddo
  200  continue
       if(allocated(klist)) deallocate(klist)
       do i=1,na
       dir(i,1)=dir(i,1)-anint(dir(i,1))
       if(dir(i,1) < 0.d0) dir(i,1)=dir(i,1)+1.d0
       dir(i,2)=dir(i,2)-anint(dir(i,2))
       if(dir(i,2) < 0.d0) dir(i,2)=dir(i,2)+1.d0
       dir(i,3)=dir(i,3)-anint(dir(i,3))
       if(dir(i,3) < 0.d0) dir(i,3)=dir(i,3)+1.d0
       enddo
!
       open(7,file='POSCAR',form='formatted')
       write(7,'(i9)') istructure
       write(7,'(a3)') '1.0'
       write(7,'(3f23.16)') a1(1),a1(2),a1(3)
       write(7,'(3f23.16)') a2(1),a2(2),a2(3)
       write(7,'(3f23.16)') a3(1),a3(2),a3(3)
       write(7,'(20(2x,a2,1x))') (symbl(i),i=1,nspecies)
       write(7,'(20(i4,1x))') (nelements(i),i=1,nspecies)
       write(7,'(a6)') "Direct"
       do i=1,na
       write(7,'(3f20.16)') dir(i,1),dir(i,2),dir(i,3)
       enddo
       close(7)
!
  991  continue
       if(allocated(itype)) deallocate(itype)
       if(allocated(dir0)) deallocate(dir0)
       if(allocated(sigmamatrix0)) deallocate(sigmamatrix0)
       end
!234567890
!      Written by In-Ho Lee, KRISS, January 28, 2013.
       subroutine tolaty(na,a1,a2,a3,dir)
       implicit none
       integer na
       real*8 a1(3),a2(3),a3(3),dir(na,3)
       real*8 b(3,3),devid,d1,d2,d3
       integer j

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
       do j=1,na
       d1=b(1,1)*dir(j,1)+b(1,2)*dir(j,2)+b(1,3)*dir(j,3)
       d2=b(2,1)*dir(j,1)+b(2,2)*dir(j,2)+b(2,3)*dir(j,3)
       d3=b(3,1)*dir(j,1)+b(3,2)*dir(j,2)+b(3,3)*dir(j,3)
       dir(j,1)=d1 ; dir(j,2)=d2 ; dir(j,3)=d3
       enddo
       do j=1,na
       dir(j,1)=dir(j,1)-anint(dir(j,1))
       dir(j,2)=dir(j,2)-anint(dir(j,2))
       dir(j,3)=dir(j,3)-anint(dir(j,3))
       enddo
       do j=1,na
       if(dir(j,1) < 0.d0) dir(j,1)=dir(j,1)+1.d0
       if(dir(j,2) < 0.d0) dir(j,2)=dir(j,2)+1.d0
       if(dir(j,3) < 0.d0) dir(j,3)=dir(j,3)+1.d0
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
