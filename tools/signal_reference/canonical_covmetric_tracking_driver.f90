! SPDX-License-Identifier: MIT
! Text-protocol driver for independent NumPy/SciPy subspace-tracking checks.
program canonical_covmetric_tracking_driver
    use scifort_kinds, only : dp
    use scifort_signal_canonical_coherence, only : canonical_coherence_result
    use scifort_signal_canonical_subspace_tracking, only : track_canonical_subspaces, &
        canonical_subspace_tracks_result, canonical_subspace_success
    implicit none
    type(canonical_coherence_result) :: source
    type(canonical_subspace_tracks_result) :: tracked
    integer :: cases, case_id, nf, nm, na, nb, f, k, i
    integer, allocatable :: mask(:)
    real(dp) :: gap, minimum

    read(*,*) cases
    do case_id = 1, cases
        read(*,*) nf,nm,na,nb,gap,minimum
        allocate(source%frequencies(nf),source%valid(nf),source%squared_modes(nf,nm), &
            source%coefficients_a(nf,na,nm),source%coefficients_b(nf,nb,nm), &
            source%covariance_a(nf,na,na),source%covariance_b(nf,nb,nb),mask(nf))
        read(*,*) source%frequencies
        read(*,*) mask
        source%valid = (mask /= 0)
        source%status = 0
        do f = 1, nf
            read(*,*) source%squared_modes(f,:)
            do k = 1, nm
                do i = 1, na
                    read(*,*) source%coefficients_a(f,i,k)
                end do
                do i = 1, nb
                    read(*,*) source%coefficients_b(f,i,k)
                end do
            end do
            do k = 1, na
                do i = 1, na
                    read(*,*) source%covariance_a(f,i,k)
                end do
            end do
            do k = 1, nb
                do i = 1, nb
                    read(*,*) source%covariance_b(f,i,k)
                end do
            end do
        end do
        tracked = track_canonical_subspaces(source,gap_tolerance=gap, &
            min_similarity=minimum,covariance_metric=.true.)
        if (tracked%status /= canonical_subspace_success) error stop 1
        write(*,'(2(i0,1x))') nf,nm
        do f = 1, nf
            do k = 1, nm
                write(*,'(7(i0,1x),2(es26.17e3,1x))',advance='no') &
                    f,k,tracked%track_id(f,k),tracked%block_first(f,k), &
                    tracked%block_size(f,k),merge(1,0,tracked%matched(f,k)), &
                    merge(1,0,tracked%segment_start(f,k)), &
                    tracked%principal_cosines(f,k),tracked%similarity(f,k)
                do i = 1, na
                    write(*,'(2(es26.17e3,1x))',advance='no') &
                        real(tracked%basis_a(f,i,k),dp),aimag(tracked%basis_a(f,i,k))
                end do
                do i = 1, nb
                    write(*,'(2(es26.17e3,1x))',advance='no') &
                        real(tracked%basis_b(f,i,k),dp),aimag(tracked%basis_b(f,i,k))
                end do
                write(*,*)
            end do
        end do
        deallocate(source%frequencies,source%valid,source%squared_modes, &
            source%coefficients_a,source%coefficients_b, &
            source%covariance_a,source%covariance_b,mask)
    end do
end program canonical_covmetric_tracking_driver
