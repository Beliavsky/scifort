! SPDX-License-Identifier: MIT
! Numerical reference driver: deterministic real/complex multitaper data.
program reference_signal_partial_shrinkage_cv
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_spectral_matrix, &
        multitaper_spectral_matrix_spectrogram, multitaper_matrix_result, &
        multitaper_spectral_matrix_spectrogram_result
    implicit none
    type(multitaper_matrix_result) :: static
    type(multitaper_spectral_matrix_spectrogram_result) :: moving
    real(dp) :: x(96,4), time
    complex(dp) :: z(96,4)
    character(8) :: weighting
    integer :: m, i, c, f, frame, variant

    do c = 1, 4
        do i = 1, 96
            time = real(i,dp)
            x(i,c) = sin(0.08_dp*time+0.5_dp*c) + 0.3_dp*cos(time*(0.11_dp+0.015_dp*c)) + &
                0.18_dp*sin(time*(0.33_dp+0.004_dp*c*c))
            z(i,c) = cmplx(x(i,c),0.26_dp*cos(time*(0.09_dp+0.014_dp*c)) + &
                0.17_dp*sin(time*0.23_dp+real(c,dp)),dp)
        end do
    end do
    do m = 1, 3
        select case (m)
        case (1)
            weighting = 'equal'
        case (2)
            weighting = 'eigen'
        case (3)
            weighting = 'adaptive'
        end select
        do variant = 1, 2
            if (variant == 1) then
                static = multitaper_spectral_matrix(x,3.5_dp,kmax=6,nfft=128, &
                    weighting=trim(weighting),cv_shrinkage=.true.)
                moving = multitaper_spectral_matrix_spectrogram(x,2.5_dp,32,noverlap=16,kmax=6, &
                    nfft=48,weighting=trim(weighting),cv_shrinkage=.true.)
                print '(a,1x,a,1x,a,1x,i0,1x,i0)', 'STATIC',trim(weighting),'real', &
                    size(static%power,1),1
            else
                static = multitaper_spectral_matrix(z,3.5_dp,kmax=6,nfft=129, &
                    weighting=trim(weighting),cv_shrinkage=.true.)
                moving = multitaper_spectral_matrix_spectrogram(z,2.5_dp,32,noverlap=16,kmax=6, &
                    nfft=48,weighting=trim(weighting),cv_shrinkage=.true.)
                print '(a,1x,a,1x,a,1x,i0,1x,i0)', 'STATIC',trim(weighting),'complex', &
                    size(static%power,1),1
            end if
            do f = 1, size(static%power,1)
                print '(5(es26.17,1x))', static%cv_shrinkage(f), real(static%power(f,1,1),dp), &
                    real(static%power(f,2,2),dp),real(static%power(f,1,2),dp), aimag(static%power(f,1,2))
            end do
            print '(a,1x,a,1x,a,1x,i0,1x,i0)', 'TF',trim(weighting), &
                merge('real   ','complex',variant==1),size(moving%power,1),size(moving%power,2)
            do frame = 1, size(moving%power,2)
                do f = 1, size(moving%power,1)
                    print '(5(es26.17,1x))', moving%cv_shrinkage(f,frame), &
                        real(moving%power(f,frame,1,1),dp),real(moving%power(f,frame,2,2),dp), &
                        real(moving%power(f,frame,1,2),dp),aimag(moving%power(f,frame,1,2))
                end do
            end do
        end do
    end do
end program reference_signal_partial_shrinkage_cv
