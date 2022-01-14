!234567890
!      Written by In-Ho Lee, KRISS, April 10, 2016.
!      ifort strings.f90 numeral.f poscar_xsf_set.f90 
!      This is for VASP 5.2.2
       program poscar_xsf_set
       USE strings, ONLY : parse, value
       implicit none
       integer natot,nspecies
       integer, allocatable :: nelements(:)
       character*2, allocatable :: symbl(:)
       real*8 etot,deltat,time,tmp,a1(3),a2(3),a3(3)
       real*8, allocatable :: xyz(:,:),fxyz(:,:),dir(:,:)
       character*200 str1
       character*200 args(40)
       character*20 delims
       integer ios,nargs
       logical lsw
       logical lfault1,lfault2
       integer i,ifile
       
       lfault1=.false. 
!      open(2,file='POSCAR',form='formatted')
       open(2,file='CONTCAR',form='formatted')
       read(2,*)
       read(2,*)
       read(2,*)
       read(2,*)
       read(2,*)
       read(2,'(a200)',err=931,end=939) str1
       delims=' '
       call parse(str1,delims,args,nargs)
       nspecies=nargs
       allocate(nelements(nspecies))
       allocate(symbl(nspecies))
       do i=1,nspecies
       symbl(i)=trim(adjustl(args(i)))
       enddo
       read(2,*) (nelements(i),i=1,nspecies)
  931  continue
       lfault1=.true. 
  939  continue
       close(2)
!
       natot=sum(nelements)
       allocate(dir(natot,3))
       allocate(xyz(natot,3))
       allocate(fxyz(natot,3))
       deltat=2.d0  *(1.d-3)   ! in ps unit
       lsw=.false.
       time=0.d0
!
       ifile=0
       lfault2=.false.
       open(1,file='OUTCAR',form='formatted')
       do 
       read(1,'(a200)',err=921,end=929) str1
       delims=' '
       call parse(str1,delims,args,nargs)
       if(nargs == 6)then
       if(args(1) == 'direct')    then
       if(args(2) == 'lattice')   then
       if(args(3) == 'vectors')   then
       if(args(4) == 'reciprocal')then
       if(args(5) == 'lattice')   then
       if(args(6) == 'vectors')   then
       read(1,*) a1(1),a1(2),a1(3)
       read(1,*) a2(1),a2(2),a2(3)
       read(1,*) a3(1),a3(2),a3(3)
!      write(6,*) a1(1),a1(2),a1(3)
!      write(6,*) a2(1),a2(2),a2(3)
!      write(6,*) a3(1),a3(2),a3(3)
                                  endif
                                  endif
                                  endif
                                  endif
                                  endif
                                  endif
                     endif
       if(nargs == 3)then
       if(args(1) == 'POSITION')   then
       if(args(2) == 'TOTAL-FORCE')then
       if(args(3) == '(eV/Angst)') then
       read(1,*)
       do i=1,natot
       read(1,*) xyz(i,1),xyz(i,2),xyz(i,3), fxyz(i,1),fxyz(i,2),fxyz(i,3)
       enddo
                                   endif
                                   endif
                                   endif
                     endif
       if(nargs == 5) then
       if(args(1) == 'total'      )then
       if(args(2) == 'drift:'     )then
       lsw=.true.
                                   endif
                                   endif
                      endif
       if(lsw)then
       if(nargs ==     7        )then
       if(args(1) == 'energy'     )then
       if(args(2) == 'without'   )then
       if(args(3) == 'entropy='    )then
       if(args(5) == 'energy(sigma->0)'        )then
       if(args(6) == '='        )then
       time=time+deltat
       call value(args(7),etot,ios)
       ifile=ifile+1
!      tmp=etot/float(natot)
       tmp=etot
!      write(6,*) time,tmp
       dir=xyz
       call tolaty(dir,a1,a2,a3,natot)
!      call gen_poscar(ifile,etot,nspecies,nelements,symbl,natot,xyz,fxyz,dir,a1,a2,a3)
       call gen_xsf(ifile,etot,nspecies,nelements,symbl,natot,xyz,fxyz,a1,a2,a3)
       lsw=.false.
                                 endif
                                 endif
                                 endif
                                 endif
                                 endif
                                 endif
              endif
       enddo
  921  continue
       lfault2=.true. 
  929  continue
       close(1)
       deallocate(xyz,fxyz,dir)
       deallocate(nelements) ; deallocate(symbl)
       stop
       end program poscar_xsf_set
!234567890
!      Written by In-Ho Lee, KRISS, August 29, 2014.
       subroutine gen_poscar(ifile,etot,nspecies,nelements,symbl,natot,xyz,fxyz,dir,a1,a2,a3)
       implicit none
       integer ifile,nspecies,natot,nelements(nspecies)
       character*2 symbl(nspecies)
       real*8 etot,xyz(natot,3),fxyz(natot,3),a1(3),a2(3),a3(3),dir(natot,3)
       integer i,isize
       integer imode
       character*280 fname

       isize=5
       call xnumeral(ifile,fname,isize)
       fname='POSCAR_'//trim(fname)
       imode=1
       imode=0
       open(11,file=trim(fname),form='formatted')
       write(11,'(f22.8)') etot
       write(11,*) '1.0'
       write(11,'(3f22.12)') a1(1),a1(2),a1(3)
       write(11,'(3f22.12)') a2(1),a2(2),a2(3)
       write(11,'(3f22.12)') a3(1),a3(2),a3(3)
       write(11,'(10(1x,a2))') (symbl(i),i=1,nspecies)
       write(11,'(10(i4))') (nelements(i),i=1,nspecies)
       if(imode == 0)then
       write(11,*) 'cartesian'
       do i=1,natot
       write(11,'(3f22.12)') xyz(i,1),xyz(i,2),xyz(i,3)
       enddo
                     endif
       if(imode == 1)then
       write(11,*) 'direct'
       do i=1,natot
       write(11,'(3f22.12)') dir(i,1),dir(i,2),dir(i,3)
       enddo
                     endif
       close(11)
       end
!234567890
!      Written by In-Ho Lee, KRISS, April 10, 2016.
       subroutine gen_xsf(ifile,etot,nspecies,nelements,symbl,natot,xyz,fxyz,a1,a2,a3)
       implicit none
       integer ifile,nspecies,natot,nelements(nspecies)
       character*2 symbl(nspecies)
       real*8 etot,xyz(natot,3),fxyz(natot,3),a1(3),a2(3),a3(3)
       integer i,j,k,isize
       character*280 fname

       isize=5
       call xnumeral(ifile,fname,isize)
       fname=trim(fname)//'.xsf'
       fname=trim(fname)
       open(11,file=trim(fname),form='formatted')
!      aenet-xsf format
       write(11,'(a17,1x,f22.8,1x,a2)') '# total energy = ',etot,'eV'
       write(11,'(a7)') 'CRYSTAL'
       write(11,'(a7)') 'PRIMVEC'
       write(11,'(3f22.12)') a1(1),a1(2),a1(3)
       write(11,'(3f22.12)') a2(1),a2(2),a2(3)
       write(11,'(3f22.12)') a3(1),a3(2),a3(3)
       write(11,'(a9)') 'PRIMCOORD'
       write(11,'(i6,1x,i1)') natot, 1
       k=0
       do i=1,nspecies
       do j=1,nelements(i)
       k=k+1
!      write(11,'(1x,a2,1x,3f22.12)') symbl(i),xyz(k,1),xyz(k,2),xyz(k,3)
       write(11,'(1x,a2,1x,6f22.12)') symbl(i),xyz(k,1),xyz(k,2),xyz(k,3),fxyz(k,1),fxyz(k,2),fxyz(k,3)
       enddo
       enddo
       close(11)
       end
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
