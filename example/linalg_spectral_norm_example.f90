! SPDX-License-Identifier: MIT
program linalg_spectral_norm_example
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_spectral_norm, only : estimate_spectral_norm, &
        estimate_spectral_norm_diff, estimate_spectral_norm_operator
    implicit none
    real(dp) :: a(3,2), b(3,2), norm_a, norm_diff, norm_operator
    integer :: status

    ! The columns of A are orthogonal, with norms 5 and 2.
    a(:,1) = [3.0_dp,4.0_dp,0.0_dp]
    a(:,2) = [0.0_dp,0.0_dp,2.0_dp]
    b = 0.0_dp
    b(:,1) = a(:,1)

    call estimate_spectral_norm(a,norm_a,status,seed=42_int64)
    if (status /= linalg_status_success) error stop 'spectral norm failed'
    call estimate_spectral_norm_diff(a,b,norm_diff,status,seed=42_int64)
    if (status /= linalg_status_success) error stop 'spectral difference failed'
    call estimate_spectral_norm_operator(3,2,matvec,rmatvec,norm_operator, &
        status,seed=42_int64)
    if (status /= linalg_status_success) error stop 'operator norm failed'

    print '(a,es16.8)', '||A||_2 (expected 5): ', norm_a
    print '(a,es16.8)', '||A-B||_2 (expected 2): ', norm_diff
    print '(a,es16.8)', 'matrix-free ||A||_2: ', norm_operator
contains
    subroutine matvec(x,y,status)
        real(dp), intent(in) :: x(:) !! input column vector
        real(dp), intent(out) :: y(:) !! output row vector
        integer, intent(out) :: status !! zero for success
        y = matmul(a,x)
        status = 0
    end subroutine matvec
    subroutine rmatvec(x,y,status)
        real(dp), intent(in) :: x(:) !! input row vector
        real(dp), intent(out) :: y(:) !! output column vector
        integer, intent(out) :: status !! zero for success
        y = matmul(transpose(a),x)
        status = 0
    end subroutine rmatvec
end program linalg_spectral_norm_example
