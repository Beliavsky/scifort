program test_signal_place_poles_mimo
    use scifort_kinds, only : dp
    use scifort_signal, only : place_poles, pole_placement_result, &
        pole_status_success, pole_status_invalid_input, pole_status_uncontrollable
    implicit none
    real(dp) :: a(3,3), b(3,2), b3(3,3), bad(3,2)
    real(dp) :: a4(4,4), b4(4,2)
    complex(dp) :: complex_roots(3), roots4(4)
    type(pole_placement_result) :: fit

    a = 0.0_dp
    a(1,2) = 1.0_dp
    a(2,3) = 1.0_dp
    b = 0.0_dp
    b(2,1) = 1.0_dp
    b(3,2) = 1.0_dp

    call check_poles(a,b,cmplx([-1.0_dp,-2.0_dp,-3.0_dp],0.0_dp,dp))
    call check_poles(a,b,cmplx([-1.0_dp,-1.0_dp,-3.0_dp],0.0_dp,dp))
    call check_poles(a,b,cmplx([0.0_dp,-1.0_dp,-3.0_dp],0.0_dp,dp))
    complex_roots = [cmplx(-1.0_dp,2.0_dp,dp), cmplx(-1.0_dp,-2.0_dp,dp), &
        cmplx(-4.0_dp,0.0_dp,dp)]
    call check_poles(a,b,complex_roots)
    complex_roots = [cmplx(-4.0_dp,0.0_dp,dp), cmplx(-1.0_dp,-2.0_dp,dp), &
        cmplx(-1.0_dp,2.0_dp,dp)]
    call check_poles(a,b,complex_roots)

    fit = place_poles(a,b,[-1.0_dp,-1.0_dp,-1.0_dp])
    if (fit%status /= pole_status_invalid_input) error stop 'excess multiplicity'
    fit = place_poles(a,b,[cmplx(-1.0_dp,1.0_dp,dp), &
        cmplx(-2.0_dp,-1.0_dp,dp), cmplx(-3.0_dp,0.0_dp,dp)])
    if (fit%status /= pole_status_invalid_input) error stop 'unpaired complex roots'

    ! Duplicate B columns give only one independent input, but the
    ! controllable pair still admits three distinct real assigned poles.
    bad = 0.0_dp
    bad(3,1) = 1.0_dp
    bad(3,2) = 1.0_dp
    call check_poles(a,bad,cmplx([-1.0_dp,-2.0_dp,-3.0_dp],0.0_dp,dp))
    fit = place_poles(a,bad,[-1.0_dp,-1.0_dp,-2.0_dp])
    if (fit%status /= pole_status_invalid_input) error stop 'effective input rank'

    ! A 3-by-3 rank-2 B must use the nullspace solver even though m>=n.
    b3(:,1:2) = b
    b3(:,3) = b(:,1)+b(:,2)
    call check_poles(a,b3,complex_roots)
    call check_poles(a,b3,cmplx([-1.0_dp,-2.0_dp,-3.0_dp],0.0_dp,dp))

    a4 = 0.0_dp
    a4(1,2) = 1.0_dp
    a4(3,4) = 1.0_dp
    b4 = 0.0_dp
    b4(2,1) = 1.0_dp
    b4(4,2) = 1.0_dp
    roots4 = [cmplx(-1.0_dp,2.0_dp,dp), cmplx(-1.0_dp,-2.0_dp,dp), &
        cmplx(-1.0_dp,2.0_dp,dp), cmplx(-1.0_dp,-2.0_dp,dp)]
    call check_poles(a4,b4,roots4)
    call check_poles(a4,b4,cmplx([-1.0_dp,-1.0_dp,-2.0_dp,-2.0_dp],0.0_dp,dp))
    fit = place_poles(a4,b4,[-1.0_dp,-1.0_dp,-1.0_dp,-2.0_dp])
    if (fit%status /= pole_status_invalid_input) error stop 'n4 excessive repeats'

    bad = 0.0_dp
    fit = place_poles(a,bad,[-1.0_dp,-2.0_dp,-3.0_dp])
    if (fit%status /= pole_status_uncontrollable) error stop 'zero B'
    a = 0.0_dp
    a(1,1) = 1.0_dp
    a(2,2) = 2.0_dp
    a(3,3) = 3.0_dp
    bad(1,1) = 1.0_dp
    bad(2,2) = 1.0_dp
    fit = place_poles(a,bad,[-1.0_dp,-2.0_dp,-3.0_dp])
    if (fit%status == pole_status_success) error stop 'uncontrollable mode moved'
    print *, 'signal MIMO place poles passed'
contains
    subroutine check_poles(aa,bb,roots)
        real(dp), intent(in) :: aa(:,:) !! test state matrix
        real(dp), intent(in) :: bb(:,:) !! test input matrix
        complex(dp), intent(in) :: roots(:) !! desired closed-loop poles
        type(pole_placement_result) :: result
        real(dp), allocatable :: closed(:,:)
        complex(dp), allocatable :: desired(:)
        real(dp), allocatable :: got(:), t(:,:), work(:,:)
        real(dp) :: c, scale
        integer :: n, i, j

        n = size(aa,1)
        result = place_poles(aa,bb,roots)
        if (result%status /= pole_status_success) then
            print *, 'MIMO status=', result%status, ' poles=', roots
            error stop 'MIMO failed status'
        end if
        if (size(result%gain_matrix,1) /= size(bb,2) .or. &
            size(result%gain_matrix,2) /= n) error stop 'MIMO gain shape'
        if (maxval(abs(result%requested_poles-roots)) > 0.0_dp) error stop 'requested poles'
        closed = aa-matmul(bb,result%gain_matrix)
        allocate(work(n,n),t(n,n),got(n+1),desired(n+1))
        work = 0.0_dp
        do i=1,n
            work(i,i) = 1.0_dp
        end do
        got = 0.0_dp
        got(1) = 1.0_dp
        do i=1,n
            t = matmul(closed,work)
            c = 0.0_dp
            do j=1,n
                c = c-t(j,j)
            end do
            c = c/real(i,dp)
            got(i+1) = c
            work = t
            do j=1,n
                work(j,j) = work(j,j)+c
            end do
        end do
        desired = (0.0_dp,0.0_dp)
        desired(1) = (1.0_dp,0.0_dp)
        do i=1,n
            do j=i+1,2,-1
                desired(j) = desired(j)-roots(i)*desired(j-1)
            end do
        end do
        scale = max(1.0_dp,maxval(abs(desired)))
        if (maxval(abs(cmplx(got,0.0_dp,dp)-desired)) > 1.0e-7_dp*scale) then
            print *, 'wanted:', desired
            print *, 'computed coeffs:', got
            error stop 'MIMO characteristic polynomial'
        end if
    end subroutine check_poles
end program test_signal_place_poles_mimo
