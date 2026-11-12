! SPDX-License-Identifier: MIT
! Test selected tridiagonal DPSS eigenpairs and legacy QL parity.
program test_signal_dpss_partial
    use scifort_kinds, only : dp
    use scifort_signal_dpss, only : dpss, dpss_result, &
        dpss_status_success, dpss_status_invalid_input
    implicit none
    type(dpss_result) :: selected, ql, automatic
    integer :: m, k, j, h
    real(dp) :: overlap

    do m = 32, 128, 32
        do k = 1, 6
            selected = dpss(m, 3.5_dp, kmax=k, solver='partial')
            ql = dpss(m, 3.5_dp, kmax=k, solver='ql')
            if (selected%status /= dpss_status_success) error stop 'partial solver failed'
            if (ql%status /= dpss_status_success) error stop 'ql solver failed'
            if (maxval(abs(selected%tapers-ql%tapers)) > 2.0e-9_dp) then
                error stop 'partial versus QL mismatch'
            end if
            if (maxval(abs(selected%ratios-ql%ratios)) > 2.0e-10_dp) then
                error stop 'concentration ratio mismatch'
            end if
            do j = 1, k
                do h = 1, k
                    overlap = dot_product(selected%tapers(:,j),selected%tapers(:,h))
                    if (j == h) overlap = overlap - 1.0_dp
                    if (abs(overlap) > 5.0e-11_dp) error stop 'nonorthogonal partial tapers'
                end do
            end do
        end do
    end do

    selected = dpss(513, 4.0_dp, kmax=6, sym=.false., &
        norm='subsample', solver='partial')
    ql = dpss(513, 4.0_dp, kmax=6, sym=.false., &
        norm='subsample', solver='ql')
    if (selected%status /= 0 .or. ql%status /= 0) error stop 'periodic dpss status'
    if (maxval(abs(selected%tapers-ql%tapers)) > 1.0e-7_dp) then
        error stop 'periodic subsample parity'
    end if
    if (maxval(abs(selected%ratios-ql%ratios)) > 1.0e-9_dp) then
        error stop 'periodic concentration parity'
    end if
    automatic = dpss(513, 4.0_dp, kmax=6, sym=.false., norm='subsample')
    if (automatic%status /= 0) error stop 'auto dpss status'
    if (maxval(abs(automatic%tapers-selected%tapers)) > 1.0e-13_dp) then
        error stop 'auto did not use selected solver'
    end if

    selected = dpss(64, 3.5_dp, kmax=6, solver='invalid')
    if (selected%status /= dpss_status_invalid_input) error stop 'invalid solver accepted'
    selected = dpss(2048, 3.5_dp, kmax=5, solver='partial')
    if (selected%status /= 0) error stop 'long partial DPSS failed'
    if (maxval(abs(matmul(transpose(selected%tapers),selected%tapers) - &
        unit_matrix(5))) > 5.0e-10_dp) error stop 'long tapers not orthonormal'
    if (any(selected%ratios < 0.0_dp) .or. any(selected%ratios > 1.0_dp)) then
        error stop 'invalid spectral concentration ratios'
    end if
    print '(a)', 'Selected tridiagonal DPSS eigenpair tests passed'
contains
    pure function unit_matrix(n) result(a)
        integer, intent(in) :: n !! order of identity matrix
        real(dp) :: a(n,n)
        integer :: i
        a = 0.0_dp
        do i = 1,n
            a(i,i) = 1.0_dp
        end do
    end function unit_matrix
end program test_signal_dpss_partial
