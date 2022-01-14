!23456789
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine init_seed()
       implicit none
       integer n,ival(8),jv(3),i
       integer, allocatable :: kseed(:)

       call date_and_time(values=ival)
       jv(1) = ival(8) + 2048*ival(7)
       jv(2) = ival(6) + 64*ival(5) ! value(4) isn't really random
       jv(3) = ival(3) + 32*ival(2) + 32*8*ival(1)
       call random_seed(size=n)
       allocate(kseed(n))
       call random_seed()          ! Give the seed an implementation-dependent kick
       call random_seed(get=kseed)
       do i=1, n
       kseed(i) = kseed(i) + jv(mod(i-1, 3) + 1)
       enddo
       call random_seed(put=kseed)
       deallocate(kseed)
       end
!23456789
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine f90sleep(dtinsec)
       implicit none
       real*8 :: dtinsec         ! desired sleep interval [s]
       integer,dimension(8) :: t ! arguments for date_and_time
       integer :: ms1,ms2  ! start and end times [ms]
       real*8 :: dt              ! desired sleep interval [ms]
 
       dt=dtinsec*1.d3
       call date_and_time(values=t)
       ms1=(t(5)*3600+t(6)*60+t(7))*1000+t(8)

       do
         call date_and_time(values=t)
         ms2=(t(5)*3600+t(6)*60+t(7))*1000+t(8)
         if(ms2-ms1>=dt)exit
       enddo
       end
