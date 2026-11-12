! SPDX-License-Identifier: MIT
! SciPy comparison cases for tight/unitary STFT windows.
program tight_compare
    use scifort_kinds, only: dp
    use scifort_signal_short_time_fft, only: short_time_fft, short_fft_success
    implicit none
    type(short_time_fft) :: s
    real(dp), allocatable :: win(:)
    character(12):: scale
    integer :: n,h,k,j,c,nf
    do n=5,12
        allocate(win(n))
        do j=1,n
            win(j)=exp(-((real(j-1,dp)-(n-1)*0.5_dp)/2.0_dp)**2)
        end do
        do h=2,min(4,n-1)
            do c=1,4
                select case(c)
                case (1)
                    scale='none'
                case (2)
                    scale='unitary'
                case (3)
                    scale='magnitude'
                case (4)
                    scale='psd'
                end select
                nf=n+3
                call s%initialize_equal_dual(win,h,8.0_dp,mfft=nf, &
                    fft_mode='twosided', scale_to=trim(scale))
                if (s%status/=short_fft_success) error stop 'failure'
                print '(5(i0,1x),a)', n,h,nf,c,s%status,trim(scale)
                do k=1,n
                    print '(4(es24.16,1x))',real(s%win(k),dp),aimag(s%win(k)), &
                        real(s%dual_win(k),dp),aimag(s%dual_win(k))
                end do
            end do
        end do
        deallocate(win)
    end do
end program tight_compare
