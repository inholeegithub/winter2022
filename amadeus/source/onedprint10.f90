!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine onedprint10(xxx,npt)
       implicit none
       integer npt
       real*8 xxx(npt)
       integer i,j,k,irem

       do i=1,npt/10
       j=(i-1)*10
       write(6,112) (xxx(j+k),k=1,10)
       enddo
       irem=npt-(npt/10)*10
       j=(npt/10)*10
       write(6,112) (xxx(j+k),k=1,irem)
!112   format(10e11.4)
 112   format(10e12.4)

       return
       end
!
