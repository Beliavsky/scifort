! SPDX-License-Identifier: MIT
program linalg_polar_frechet_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_polar, only : polar
    use scifort_linalg_expm_frechet, only : expm_frechet
    implicit none
    real(dp) :: a(3,2), u(3,2), p(2,2), polar_error, orth_error
    complex(dp) :: matrix(2,2), direction(2,2), exponential(2,2), derivative(2,2)
    integer :: status

    a(:,1) = [1.0_dp, 2.0_dp, 0.0_dp]
    a(:,2) = [2.0_dp, 0.0_dp, 1.0_dp]
    call polar(a, u, p, status, residual=polar_error, orthogonality=orth_error)
    if (status /= linalg_status_success) error stop 'polar decomposition failed'
    print '(a,es12.4)', 'polar reconstruction error: ', polar_error
    print '(a,es12.4)', 'polar orthogonality error: ', orth_error

    matrix(1,:) = [cmplx(0.2_dp, 0.3_dp, dp), cmplx(1.0_dp,0.0_dp,dp)]
    matrix(2,:) = [cmplx(0.0_dp,0.0_dp,dp), cmplx(-0.4_dp,0.2_dp,dp)]
    direction = cmplx(0.0_dp,0.0_dp,dp)
    direction(2,1) = cmplx(1.0_dp,-0.5_dp,dp)
    call expm_frechet(matrix,direction,exponential,derivative,status)
    if (status /= linalg_status_success) error stop 'expm_frechet failed'
    print '(a,2es16.7)', 'Frechet derivative element (1,1): ', &
        real(derivative(1,1),dp), aimag(derivative(1,1))
end program linalg_polar_frechet_example
