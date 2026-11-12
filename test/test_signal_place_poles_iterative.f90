program test_signal_place_poles_iterative
    use scifort_kinds, only : dp
    use scifort_signal, only : place_poles, pole_placement_result, &
        pole_status_success, pole_status_invalid_input
    implicit none
    real(dp) :: a(4,4), b(4,2), old_det
    complex(dp) :: complex_poles(4)
    type(pole_placement_result) :: fit

    a(1,:) = [1.380_dp, -0.2077_dp, 6.715_dp, -5.676_dp]
    a(2,:) = [-0.5814_dp, -4.290_dp, 0.0_dp, 0.675_dp]
    a(3,:) = [1.067_dp, 4.273_dp, -6.654_dp, 5.893_dp]
    a(4,:) = [0.048_dp, 4.273_dp, 1.343_dp, -2.104_dp]
    b(:,1) = [0.0_dp, 1.136_dp, 0.0_dp, -3.146_dp]
    b(:,2) = [5.679_dp, 1.136_dp, 0.0_dp, 0.0_dp]

    fit = place_poles(a,b,[-0.2_dp,-0.5_dp,-5.0566_dp,-8.6659_dp], &
        method='KNV0',rtol=1.0e-3_dp,maxiter=30)
    call check(fit,cmplx([-0.2_dp,-0.5_dp,-5.0566_dp,-8.6659_dp],0.0_dp,dp))
    if (fit%nb_iter < 1 .or. fit%nb_iter > 30) error stop 'knv0 iteration count'
    if (fit%det_transfer <= 0.0_dp) error stop 'knv0 determinant'
    fit = place_poles(a,b,[-0.2_dp,-0.5_dp,-5.0566_dp,-8.6659_dp], &
        method='YT',rtol=-1.0_dp,maxiter=5)
    call check(fit,cmplx([-0.2_dp,-0.5_dp,-5.0566_dp,-8.6659_dp],0.0_dp,dp))
    if (fit%nb_iter /= 5) error stop 'negative rtol forces full iteration budget'
    old_det = fit%det_transfer
    fit = place_poles(a,b,[-0.2_dp,-0.5_dp,-5.0566_dp,-8.6659_dp], &
        method='yt',rtol=-1.0_dp,maxiter=1)
    call check(fit,cmplx([-0.2_dp,-0.5_dp,-5.0566_dp,-8.6659_dp],0.0_dp,dp))
    if (fit%det_transfer > old_det*(1.0_dp+1.0e-8_dp)) &
        error stop 'yt determinant must not decrease across iterations'

    complex_poles = [cmplx(-1.0_dp,0.0_dp,dp), cmplx(-3.0_dp,0.0_dp,dp), &
        cmplx(-2.0_dp,1.0_dp,dp), cmplx(-2.0_dp,-1.0_dp,dp)]
    fit = place_poles(a,b,complex_poles,method='YT',maxiter=15)
    call check(fit,complex_poles)
    if (fit%nb_iter < 1 .or. fit%nb_iter > 15) error stop 'complex YT iterations'
    fit = place_poles(a,b,complex_poles,method='KNV0')
    if (fit%status /= pole_status_invalid_input) error stop 'knv0 rejects complex'
    fit = place_poles(a,b,real(complex_poles,dp),method='bogus')
    if (fit%status /= pole_status_invalid_input) error stop 'invalid method'
    fit = place_poles(a,b,real(complex_poles,dp),maxiter=0)
    if (fit%status /= pole_status_invalid_input) error stop 'invalid maxiter'
    fit = place_poles(a,b,real(complex_poles,dp),rtol=2.0_dp)
    if (fit%status /= pole_status_invalid_input) error stop 'invalid rtol'
    print *, 'iterative pole placement passed'
contains
    subroutine check(res, poles)
        type(pole_placement_result), intent(in) :: res !! fitted gains and diagnostics
        complex(dp), intent(in) :: poles(:) !! requested poles
        real(dp) :: closed(4,4), pchar(5), id(4,4), t(4,4), c
        complex(dp) :: target(5)
        integer :: i, j
        if (res%status /= pole_status_success) error stop 'pole placement failed'
        if (.not. allocated(res%gain_matrix)) error stop 'gain allocation'
        if (res%det_transfer <= 0.0_dp) error stop 'nonpositive transfer determinant'
        if (res%rtol < 0.0_dp .or. res%rtol > 1.0_dp) error stop 'rtol range'
        closed = a - matmul(b,res%gain_matrix)
        id = 0.0_dp
        do i=1,4
            id(i,i) = 1.0_dp
        end do
        t = id
        pchar = 0.0_dp
        pchar(1) = 1.0_dp
        do i=1,4
            t = matmul(closed,t)
            c = 0.0_dp
            do j=1,4
                c = c-t(j,j)
            end do
            c = c/real(i,dp)
            pchar(i+1) = c
            t = t+c*id
        end do
        target = (0.0_dp,0.0_dp)
        target(1) = (1.0_dp,0.0_dp)
        do i=1,4
            do j=i+1,2,-1
                target(j) = target(j)-poles(i)*target(j-1)
            end do
        end do
        if (maxval(abs(cmplx(pchar,0.0_dp,dp)-target)) > 1.0e-6_dp) &
            error stop 'closed loop eigenvalue polynomial mismatch'
    end subroutine check
end program test_signal_place_poles_iterative
