!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine object(egp1,egp2,obj,mcmd,lfault_stdout)
       implicit none
       logical lfault_stdout
       integer mcmd
       real*8 obj,egp2,egp1,penalty
       integer kcase

! default case, mcmd == 0   : enthalpy minimization 
! default case, mcmd == 4   : electronic DOS derived effective mass maximization 
!               mcmd == 1   : direct band gap optimization
!               mcmd == 2   : electronic DOS at Fermi level maximization
!               mcmd == 3   : electronic DOS slope at Fermi level maximization
!               mcmd == 5   : effective mass minimization
!               mcmd == 6   : special bond length preference
!               mcmd == 7   : xrd simulation vs experimental data
!
       if(mcmd ==1)then
        if(egp1 >= egp2)then
       obj=-egp1
                        else
       obj=abs(egp2-egp1)-egp1
                        endif
       kcase=0
       kcase=1
       if(kcase == 1)then
       penalty=0.d0
       if(egp2 < 0.5d0) penalty=0.5d0-egp2
       if(egp2 > 0.8d0) penalty=-0.8d0+egp2
       obj=obj+penalty
                     endif
                   endif
       if(mcmd ==2)then
       obj=-abs(obj)
                   endif
       if(mcmd ==3)then
       obj=-abs(obj)
                   endif
!      if(mcmd ==5)then
!      obj=min(egp1,egp2)
!                  endif

       if(lfault_stdout) obj=1.d19
       write(6,'(e22.10,2x,a7)') obj, 'obj ftn'
       return
       end
!
