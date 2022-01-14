!234567890
!      Written by In-Ho Lee, KRISS, August 29, 2014.
!      ifort strings.f90 etot_plot.f90
       program etot_plot
       USE strings, ONLY : parse, value
       implicit none
       integer natom
       real*8 etot,deltat,time,tmp
       character*200 str1
       character*200 args(40)
       character*20 delims
       integer ios,nargs
       logical lsw
       
!
       natom=128
       deltat=2.d0  *(1.d-3)   ! in ps unit
!
       open(1,file='OUTCAR',form='formatted')
       lsw=.false.
       time=0.d0
       do 
       read(1,'(a200)',err=911,end=999) str1
       delims=' '
       call parse(str1,delims,args,nargs)
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
       tmp=etot/float(natom)
       write(6,*) time,tmp
       lsw=.false.
                                 endif
                                 endif
                                 endif
                                 endif
                                 endif
                                 endif
              endif
       enddo
  911  continue
  999  continue
       close(1)
       stop
       end program etot_plot
