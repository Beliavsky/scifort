! SPDX-License-Identifier: MIT
! Independent-output driver for NumPy reference comparison.
program reference_signal_multiple_coherence
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_spectral_matrix_spectrogram_result
    use scifort_signal_multiple_coherence, only : multiple_coherence_tf_result, &
        multiple_coherence_tf_matrix
    use scifort_signal_partial_coherence, only : partial_coherence_success
    implicit none
    integer, parameter :: nf=7, nt=4, nc=5
    type(multitaper_spectral_matrix_spectrogram_result) :: mat
    type(multiple_coherence_tf_result) :: result
    complex(dp) :: z(nc,nc), v
    integer :: f, t, i, j, k, mode
    real(dp) :: angle, alpha, ridge

    allocate(mat%frequencies(nf),mat%times(nt),mat%power(nf,nt,nc,nc), &
        mat%cv_shrinkage(nf,nt))
    mat%status = 0
    mat%effective_tapers = 5.4_dp
    mat%frequencies = [(real(f-1,dp)/real(nf,dp),f=1,nf)]
    mat%times = [(real(t,dp),t=1,nt)]
    do t=1,nt
        do f=1,nf
            do j=1,nc
                do i=1,nc
                    angle = 0.13_dp*real(i,dp)+0.17_dp*real(j,dp)+ &
                        0.11_dp*real(f,dp)+0.07_dp*real(t,dp)
                    z(i,j) = cmplx(sin(angle)+0.17_dp*real(i,dp)*merge(1.0_dp,0.0_dp,i == j), &
                        cos(2.0_dp*angle)+0.11_dp*real(j,dp),dp)
                end do
            end do
            do j=1,nc
                do i=1,nc
                    v = (0.0_dp,0.0_dp)
                    do k=1,nc
                        v = v+conjg(z(k,i))*z(k,j)
                    end do
                    if (i == j) v = v+cmplx(0.15_dp+0.035_dp*real(i,dp),0.0_dp,dp)
                    mat%power(f,t,i,j) = v
                end do
            end do
            mat%cv_shrinkage(f,t) = 0.04_dp*real(f,dp)+0.03_dp*real(t,dp)
        end do
    end do
    do mode=1,5
        select case(mode)
        case (1)
            ridge = 0.0_dp
            alpha = 0.0_dp
            result = multiple_coherence_tf_matrix(mat,ridge=ridge,shrinkage=alpha)
        case (2)
            ridge = 0.04_dp
            alpha = 0.20_dp
            result = multiple_coherence_tf_matrix(mat,ridge=ridge,shrinkage=alpha)
        case (3)
            ridge = 0.10_dp
            alpha = 0.80_dp
            result = multiple_coherence_tf_matrix(mat,ridge=ridge,shrinkage=alpha)
        case (4)
            ridge = 0.05_dp
            result = multiple_coherence_tf_matrix(mat,ridge=ridge,cv_shrinkage=.true.)
        case (5)
            ridge = 0.02_dp
            result = multiple_coherence_tf_matrix(mat,ridge=ridge,auto_shrinkage=.true., &
                effective_tapers=mat%effective_tapers)
        end select
        if (result%status /= partial_coherence_success) &
            error stop 'multiple coherence reference failed'
        do t=1,nt
            do f=1,nf
                do i=1,nc
                    write(*,'(4(i0,1x),4(es25.16e3,1x))') mode,f,t,i, &
                        result%coherence(f,t,i),result%conditional_power(f,t,i), &
                        result%explained_power(f,t,i),result%shrinkage_used(f,t)
                end do
            end do
        end do
    end do
end program reference_signal_multiple_coherence
