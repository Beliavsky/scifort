program trial
 use scifort_kinds,only: dp
 use scifort_signal_design,only: butter_result
 use scifort_signal_iir_extra,only: buttord,cheb1ord,cheb2ord,iir_order_result,iirdesign
 implicit none
 type(butter_result):: f
 type(iir_order_result):: o
 integer :: st,style,flag
 real(dp):: gp,gs,fs,wp,ws
 do
 read(*,*,iostat=st) style,wp,ws,gp,gs,fs,flag
 if(st/=0) exit
 select case (style)
 case (1)
  o=buttord(wp,ws,gp,gs,fs)
  if (flag==1) f=iirdesign(wp,ws,gp,gs,'butter',fs)
 case (2)
  o=cheb1ord(wp,ws,gp,gs,fs)
  if (flag==1) f=iirdesign(wp,ws,gp,gs,'cheby1',fs)
 case (3)
  o=cheb2ord(wp,ws,gp,gs,fs)
  if (flag==1) f=iirdesign(wp,ws,gp,gs,'cheby2',fs)
 end select
 print '(i0,1x,i0,1x,es26.17)',o%status,o%order,o%cutoff
 if (flag==1 .and. o%status==0) then
   print '(i0,1x,i0)',f%status,size(f%sos,1)
   print '(*(es26.17,1x))', f%b
   print '(*(es26.17,1x))', f%a
 end if
 end do
end program
