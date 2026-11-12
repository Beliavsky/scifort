program trial
 use scifort_kinds,only: dp
 use scifort_signal_design,only: butter_result
 use scifort_signal_iir_extra,only: cheby2,iirnotch,iirpeak,iirfilter
 implicit none
 type(butter_result):: f
 integer :: n,i,style,st,k
 real(dp):: rs,fs,c(2)
 character(20):: type
 do
 read(*,*,iostat=st) n,rs,fs,style,c
 if(st/=0) exit
 select case (style)
 case (1)
  f=cheby2(n,rs,c(1),'lowpass',fs)
 case (2)
  f=cheby2(n,rs,c(1),'highpass',fs)
 case (3)
  f=cheby2(n,rs,c(:2),'bandpass',fs)
 case (4)
  f=cheby2(n,rs,c(:2),'bandstop',fs)
 case (5)
  f=iirnotch(c(1),rs,fs)
 case (6)
  f=iirpeak(c(1),rs,fs)
 end select
 print '(3(i0,1x))',f%status, n, style
 if (f%status == 0) then
  print '(i0)',size(f%b)
  print '(*(es26.17,1x))', f%b
  print '(*(es26.17,1x))', f%a
  print '(i0)',size(f%sos,1)
  do k=1,size(f%sos,1)
   print '(*(es26.17,1x))',f%sos(k,:)
  end do
 end if
 end do
end program
