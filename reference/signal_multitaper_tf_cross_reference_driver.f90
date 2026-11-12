! SPDX-License-Identifier: MIT
! Numerical reference driver: test/compare_signal_multitaper_tf_cross.py
program signal_multitaper_tf_cross_reference_driver
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_cross_spectrogram, multitaper_tf_coherence, &
        multitaper_cross_spectrogram_result, multitaper_tf_coherence_result
    implicit none
    integer :: cases, icase, n, seg, ov, k, nf, kind, center, i, j, f, t
    real(dp) :: nw, fs, a, b, c, d
    character(20) :: mode
    complex(dp), allocatable :: x(:), y(:)
    type(multitaper_cross_spectrogram_result) :: cross
    type(multitaper_tf_coherence_result) :: coherence

    read(*,*) cases
    do icase = 1, cases
        read(*,*) n,seg,ov,k,nf,kind,center,nw,fs,mode
        allocate(x(n),y(n))
        do i = 1, n
            read(*,*) a,b,c,d
            x(i) = cmplx(a,b,dp)
            y(i) = cmplx(c,d,dp)
        end do
        select case(kind)
        case(0)
            cross = multitaper_cross_spectrogram(real(x),real(y),nw,seg,ov,k,fs,nf,center/=0,trim(mode))
            coherence = multitaper_tf_coherence(real(x),real(y),nw,seg,ov,k,fs,nf,center/=0,trim(mode))
        case(1)
            cross = multitaper_cross_spectrogram(x,y,nw,seg,ov,k,fs,nf,center/=0,trim(mode))
            coherence = multitaper_tf_coherence(x,y,nw,seg,ov,k,fs,nf,center/=0,trim(mode))
        case(2)
            cross = multitaper_cross_spectrogram(real(x),y,nw,seg,ov,k,fs,nf,center/=0,trim(mode))
            coherence = multitaper_tf_coherence(real(x),y,nw,seg,ov,k,fs,nf,center/=0,trim(mode))
        case default
            cross = multitaper_cross_spectrogram(x,real(y),nw,seg,ov,k,fs,nf,center/=0,trim(mode))
            coherence = multitaper_tf_coherence(x,real(y),nw,seg,ov,k,fs,nf,center/=0,trim(mode))
        end select
        f = size(cross%frequencies)
        t = size(cross%times)
        print '(3(i0,1x))', cross%status, f, t
        if (cross%status == 0 .and. coherence%status == 0) then
            do i = 1, f
                print '(es25.16)', cross%frequencies(i)
            end do
            do j = 1, t
                print '(es25.16)', cross%times(j)
            end do
            do j = 1, t
                do i = 1, f
                    print '(5(es25.16,1x))', real(cross%power(i,j)),aimag(cross%power(i,j)), &
                        cross%auto_x(i,j),cross%auto_y(i,j),coherence%coherence(i,j)
                end do
            end do
        end if
        deallocate(x,y)
    end do
end program signal_multitaper_tf_cross_reference_driver
