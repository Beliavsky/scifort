program test_linalg_matrix_functions
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_no_convergence, &
        linalg_status_invalid_shape
    use scifort_linalg_matrix_functions, only : sqrtm, sqrtm_real, signm
    implicit none

    real(dp) :: ar(2, 2), sr(2, 2), residual
    complex(dp) :: ac(2, 2), xc(2, 2), scomplex(2, 2), expected(2, 2)
    integer :: status, it

    ar = 0.0_dp
    ar(1, 1) = 1.0_dp
    ar(2, 2) = -4.0_dp
    call sqrtm(ar, xc, status, residual)
    call assert_ok(status == linalg_status_success .and. residual < 1.0e-12_dp, 'real sqrtm')
    call assert_ok(abs(xc(1, 1)-cmplx(1.0_dp, 0.0_dp, dp)) < 1.0e-11_dp, 'sqrtm positive')
    call assert_ok(abs(xc(2, 2)-cmplx(0.0_dp, 2.0_dp, dp)) < 1.0e-11_dp, 'sqrtm negative')
    call assert_ok(maxval(abs(matmul(xc, xc)-cmplx(ar, kind=dp))) < 1.0e-11_dp, 'sqrtm residual')

    ar = 0.0_dp
    ar(1, 1) = 4.0_dp
    ar(2, 2) = 9.0_dp
    ar(1, 2) = 3.0_dp
    call sqrtm_real(ar, sr, status, residual)
    call assert_ok(status == linalg_status_success, 'sqrtm real result')
    call assert_ok(abs(sr(1, 2)-0.6_dp) < 1.0e-11_dp, 'sqrtm upper triangular')
    call assert_ok(maxval(abs(matmul(sr, sr)-ar)) < 1.0e-11_dp, 'real sqrtm residual')

    ! Real-only output refuses a non-real principal square root.
    ar = 0.0_dp
    ar(1, 1) = 1.0_dp
    ar(2, 2) = -4.0_dp
    call sqrtm_real(ar, sr, status)
    call assert_ok(status == linalg_status_no_convergence, 'sqrtm_real nonreal result')

    ar = 0.0_dp
    ar(1, 2) = 1.0_dp
    call sqrtm(ar, xc, status)
    call assert_ok(status == linalg_status_no_convergence, 'nilpotent square root does not exist')
    ar = 0.0_dp
    call sqrtm(ar, xc, status)
    call assert_ok(status == linalg_status_success .and. all(abs(xc) == 0.0_dp), 'zero root')
    call sqrtm(ar(1:2, 1:1), xc, status)
    call assert_ok(status == linalg_status_invalid_shape, 'sqrtm nonsquare')
    call signm(ar(1:2, 1:1), sr, status)
    call assert_ok(status == linalg_status_invalid_shape, 'signm nonsquare')
    call signm(ar, sr, status)
    call assert_ok(status == linalg_status_success .and. all(abs(sr) == 0.0_dp), 'zero sign')

    ac(1, 1) = cmplx(2.0_dp, 1.0_dp, dp)
    ac(1, 2) = cmplx(0.3_dp, 0.2_dp, dp)
    ac(2, 1) = cmplx(-0.5_dp, 0.7_dp, dp)
    ac(2, 2) = cmplx(1.5_dp, -1.0_dp, dp)
    expected = ac
    ac = matmul(expected, expected)
    call sqrtm(ac, xc, status, residual)
    call assert_ok(status == linalg_status_success .and. residual < 1.0e-11_dp, 'complex sqrtm')
    call assert_ok(maxval(abs(xc-expected)) < 1.0e-10_dp, 'principal complex sqrtm')

    ar = 0.0_dp
    ar(1, 1) = -1.0_dp
    ar(2, 2) = 2.0_dp
    ar(1, 2) = 4.0_dp
    call signm(ar, sr, status, iterations=it, residual=residual)
    call assert_ok(status == linalg_status_success .and. it > 0, 'real signm')
    call assert_ok(maxval(abs(matmul(sr, sr)-identity2())) < 1.0e-10_dp, 'sign involution')
    call assert_ok(abs(sr(1, 2)-8.0_dp/3.0_dp) < 1.0e-10_dp, 'sign offdiagonal')

    ac = cmplx(ar, kind=dp)
    ac(1, 1) = cmplx(-1.0_dp, 1.0_dp, dp)
    ac(2, 2) = cmplx(2.0_dp, -3.0_dp, dp)
    call signm(ac, scomplex, status, residual=residual)
    call assert_ok(status == linalg_status_success .and. residual < 1.0e-10_dp, 'complex signm')
    call assert_ok(maxval(abs(matmul(scomplex, scomplex)-identity2())) < 1.0e-10_dp, &
        'complex sign involution')
    call signm(ac, scomplex, status, max_iter=1)
    call assert_ok(status == linalg_status_no_convergence, 'signm nonconvergence flag')
    call signm(ac, scomplex, status, tol=-1.0_dp)
    call assert_ok(status == linalg_status_invalid_shape, 'signm bad tolerance')

    print '(a)', 'matrix_functions_native PASS'
contains
    subroutine assert_ok(ok, label)
        logical, intent(in) :: ok !! test assertion
        character(*), intent(in) :: label !! failure description
        if (.not. ok) then
            print '(a)', 'FAIL: '//label
            error stop 1
        end if
    end subroutine assert_ok

    pure function identity2() result(id)
        real(dp) :: id(2, 2)
        id = 0.0_dp
        id(1, 1) = 1.0_dp
        id(2, 2) = 1.0_dp
    end function identity2
end program test_linalg_matrix_functions
