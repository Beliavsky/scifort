program test_signal_place_poles
    use scifort_kinds, only : dp
    use scifort_signal, only : place_poles, pole_placement_result, &
        pole_status_success, pole_status_invalid_input, &
        pole_status_uncontrollable
    implicit none
    type(pole_placement_result) :: fit
    real(dp) :: a2(2,2), b2(2,1), b_square(2,2)
    real(dp) :: a3(3,3), b3(3,1), badb(3,1), bsmall(3,2)

    a2(1,:) = [0.0_dp,1.0_dp]
    a2(2,:) = [-2.0_dp,-3.0_dp]
    b2(:,1) = [0.0_dp,1.0_dp]
    fit = place_poles(a2,b2,[-4.0_dp,-5.0_dp])
    if (fit%status /= pole_status_success) error stop 'siso 2 status'
    if (maxval(abs(fit%gain_matrix(1,:)-[18.0_dp,6.0_dp])) > 1.e-11_dp) error stop 'siso 2 gain'

    a3 = 0.0_dp
    a3(1,2) = 1.0_dp
    a3(2,3) = 1.0_dp
    b3 = 0.0_dp
    b3(3,1) = 1.0_dp
    fit = place_poles(a3,b3,[-1.0_dp,-2.0_dp,-3.0_dp])
    if (fit%status /= pole_status_success) error stop 'siso 3 status'
    if (maxval(abs(fit%gain_matrix(1,:)-[6.0_dp,11.0_dp,6.0_dp])) > 1.e-11_dp) &
        error stop 'siso 3 gain'
    fit = place_poles(a3,b3,[cmplx(-1.0_dp,0.0_dp,dp), &
        cmplx(-2.0_dp,3.0_dp,dp),cmplx(-2.0_dp,-3.0_dp,dp)])
    if (fit%status /= pole_status_success) error stop 'conjugate status'
    if (maxval(abs(fit%gain_matrix(1,:)-[13.0_dp,17.0_dp,5.0_dp])) > 1.e-11_dp) &
        error stop 'conjugate gain'
    b_square = 0.0_dp
    b_square(1,1) = 1.0_dp
    b_square(2,2) = 1.0_dp
    fit = place_poles(a2,b_square,[-1.0_dp,-2.0_dp])
    if (fit%status /= pole_status_success) error stop 'full actuation status'
    if (maxval(abs(matmul(b_square,fit%gain_matrix)-(a2- &
        reshape([-1.0_dp,0.0_dp,0.0_dp,-2.0_dp],[2,2])))) > 1.e-11_dp) &
        error stop 'full actuation gain'
    badb=0.0_dp
    fit = place_poles(a3,badb,[-1.0_dp,-2.0_dp,-3.0_dp])
    if (fit%status /= pole_status_uncontrollable) error stop 'uncontrollable'
    bsmall = 0.0_dp
    bsmall(2,1)=1.0_dp
    bsmall(3,2)=1.0_dp
    fit = place_poles(a3,bsmall,[-1.0_dp,-2.0_dp,-3.0_dp])
    if (fit%status /= pole_status_success) error stop 'mimo status'
    if (size(fit%gain_matrix,1) /= 2 .or. size(fit%gain_matrix,2) /= 3) error stop 'mimo shape'
    fit = place_poles(a2,b2,[-1.0_dp])
    if (fit%status /= pole_status_invalid_input) error stop 'bad shape'
    print *, 'signal place poles passed'
end program test_signal_place_poles
