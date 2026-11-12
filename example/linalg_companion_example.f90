! SPDX-License-Identifier: MIT
program linalg_companion_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_companion, only : companion, fiedler_companion
    use scifort_linalg_fiedler, only : fiedler
    implicit none
    real(dp) :: a(5), c(4,4), fc(4,4), distances(3,3)
    integer :: status, j

    ! P(x) = (x-1)(x-3)(x-5)(x-7).
    a = [1.0_dp,-16.0_dp,86.0_dp,-176.0_dp,105.0_dp]
    call companion(a,c,status)
    if (status /= linalg_status_success) error stop 'companion failed'
    call fiedler_companion(a,fc,status)
    if (status /= linalg_status_success) error stop 'fiedler_companion failed'
    print '(a)', 'Ordinary companion matrix:'
    do j=1,size(c,1)
        print '(4f10.2)', c(j,:)
    end do
    print '(a)', 'Fiedler companion matrix (same characteristic polynomial):'
    do j=1,size(fc,1)
        print '(4f10.2)', fc(j,:)
    end do

    call fiedler([1.0_dp,4.0_dp,12.0_dp],distances,status)
    if (status /= linalg_status_success) error stop 'fiedler failed'
    print '(a)', 'Symmetric Fiedler distance matrix:'
    do j=1,size(distances,1)
        print '(3f10.2)', distances(j,:)
    end do
end program linalg_companion_example
