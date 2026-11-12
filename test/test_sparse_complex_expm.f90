! SPDX-License-Identifier: MIT
program test_sparse_complex_expm
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csc_complex_matrix, &
        csr_complex_from_dense, csr_complex_to_csc
    use scifort_sparse_linalg_complex_expm, only : expm_multiply_complex, &
        expm_multiply_complex_operator, expm_multiply_complex_interval, &
        expm_multiply_complex_operator_interval, complex_expm_success, complex_expm_invalid_input, &
        complex_expm_callback_failure, complex_expm_breakdown
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none
    type(csr_complex_matrix) :: a
    type(csc_complex_matrix) :: c
    complex(dp) :: dense(3,3), rhs(3,2), expected(3,2), mu, z(3), exact(3)
    complex(dp), allocatable :: y(:,:), v(:), series(:,:,:), vector_series(:,:)
    complex(dp), allocatable :: big_diag(:), big_b(:), big_y(:)
    integer :: info, i, scaling_steps, products
    real(dp), parameter :: epsilon_test = 3.0e-12_dp
    real(dp) :: nan_value

    ! Jordan block A = mu*I+N, with N**2=0: exponential is exact.
    mu = cmplx(-0.35_dp,0.65_dp,dp)
    dense = cmplx(0.0_dp,0.0_dp,dp)
    dense(1,1) = mu
    dense(2,2) = mu
    dense(3,3) = mu
    dense(1,2) = cmplx(2.0_dp,-1.0_dp,dp)
    rhs(:,1) = [cmplx(1.0_dp,0.5_dp,dp), cmplx(-0.2_dp,1.0_dp,dp), cmplx(1.0_dp,-1.0_dp,dp)]
    rhs(:,2) = [cmplx(0.1_dp,0.2_dp,dp), cmplx(2.0_dp,-0.5_dp,dp), cmplx(0.0_dp,1.0_dp,dp)]
    expected = exp(1.25_dp*mu)*rhs
    expected(1,:) = expected(1,:) + exp(1.25_dp*mu)*1.25_dp*dense(1,2)*rhs(2,:)
    a = csr_complex_from_dense(dense,info)
    call check(info == 0,'CSR assembly')
    c = csr_complex_to_csc(a,info)
    call check(info == 0,'CSC assembly')
    call expm_multiply_complex(a,rhs,y,info,t=1.25_dp,scaling_steps=scaling_steps,products=products)
    call check(info == complex_expm_success,'CSR block action')
    call check(maxval(abs(y-expected)) < epsilon_test,'CSR block result')
    call check(products > 0 .and. scaling_steps > 0,'diagnostics')
    call expm_multiply_complex(c,rhs,y,info,t=1.25_dp,tracea=3.0_dp*mu)
    call check(info == complex_expm_success,'CSC block action')
    call check(maxval(abs(y-expected)) < epsilon_test,'CSC block result')
    call expm_multiply_complex(a,rhs(:,2),v,info,t=1.25_dp)
    call check(info == complex_expm_success,'CSR vector action')
    call check(maxval(abs(v-expected(:,2))) < epsilon_test,'CSR vector result')
    call expm_multiply_complex(c,rhs(:,1),v,info,t=1.25_dp)
    call check(info == complex_expm_success,'CSC vector action')
    call check(maxval(abs(v-expected(:,1))) < epsilon_test,'CSC vector result')

    ! Matrix-free time stepping with a documented conservative norm bound.
    call expm_multiply_complex_operator(3,rhs,apply_dense,3.0_dp,y,info,t=1.25_dp, &
        tracea=3.0_dp*mu,products=products)
    call check(info == complex_expm_success,'operator block action')
    call check(maxval(abs(y-expected)) < epsilon_test,'operator block result')
    call check(products > 0,'operator product count')
    call expm_multiply_complex_operator(3,rhs(:,1),apply_dense,3.0_dp,v,info,t=-0.4_dp)
    call check(info == complex_expm_success,'operator vector action')
    exact = exp(-0.4_dp*mu)*rhs(:,1)
    exact(1) = exact(1)-0.4_dp*exp(-0.4_dp*mu)*dense(1,2)*rhs(2,1)
    call check(maxval(abs(v-exact)) < epsilon_test,'operator negative time')

    ! Time intervals must agree with exact Jordan-block actions at every sample.
    call expm_multiply_complex_interval(a,rhs,-0.9_dp,1.1_dp,5,series,info,products=products)
    call check(info == complex_expm_success,'CSR interval action')
    call check(all(shape(series) == [3,2,5]),'CSR interval shape')
    do i = 1,5
        call check(maxval(abs(series(:,:,i)-jordan_action(-0.9_dp+0.5_dp*real(i-1,dp)))) &
                   < epsilon_test,'CSR interval result')
    end do
    call expm_multiply_complex_interval(c,rhs,-0.9_dp,1.1_dp,5,series,info)
    call check(info == complex_expm_success,'CSC interval action')
    do i = 1,5
        call check(maxval(abs(series(:,:,i)-jordan_action(-0.9_dp+0.5_dp*real(i-1,dp)))) &
                   < epsilon_test,'CSC interval result')
    end do
    call expm_multiply_complex_interval(a,rhs(:,1),-0.9_dp,1.1_dp,4,vector_series,info,endpoint=.false.)
    call check(info == complex_expm_success,'CSR endpoint-exclusive interval')
    do i = 1,4
        expected = jordan_action(-0.9_dp+0.5_dp*real(i-1,dp))
        call check(maxval(abs(vector_series(:,i)-expected(:,1))) &
                   < epsilon_test,'exclusive interval result')
    end do
    call expm_multiply_complex_operator_interval(3,rhs,apply_dense,3.0_dp,-0.9_dp,1.1_dp,5, &
        series,info,tracea=3.0_dp*mu)
    call check(info == complex_expm_success,'operator interval action')
    do i = 1,5
        call check(maxval(abs(series(:,:,i)-jordan_action(-0.9_dp+0.5_dp*real(i-1,dp)))) &
                   < epsilon_test,'operator interval result')
    end do
    call expm_multiply_complex_interval(a,rhs(:,1),0.2_dp,9.0_dp,1,vector_series,info)
    call check(info == complex_expm_success .and. size(vector_series,2) == 1,'single time interval')
    expected = jordan_action(0.2_dp)
    call check(maxval(abs(vector_series(:,1)-expected(:,1))) &
               < epsilon_test,'single interval result')
    call expm_multiply_complex_interval(a,rhs(:,1),0.2_dp,9.0_dp,0,vector_series,info)
    call check(info == complex_expm_invalid_input,'invalid interval count')

    ! Empty and zero-time cases must not invoke the operator callback.
    call expm_multiply_complex_operator(3,rhs(:,1),always_fail,3.0_dp,v,info,t=0.0_dp,products=products)
    call check(info == complex_expm_success .and. products == 0,'zero-time short circuit')
    call check(maxval(abs(v-rhs(:,1))) == 0.0_dp,'zero-time identity')
    call expm_multiply_complex_operator(3,rhs(:,1),always_fail,3.0_dp,v,info)
    call check(info == complex_expm_callback_failure,'callback failure')
    call expm_multiply_complex_operator(3,rhs(:,1),bad_data,3.0_dp,v,info)
    call check(info == complex_expm_breakdown,'nonfinite callback result')
    call expm_multiply_complex_operator(3,rhs(:,1),apply_dense,-1.0_dp,v,info)
    call check(info == complex_expm_invalid_input,'invalid norm bound')
    call expm_multiply_complex(a,rhs(:,1),v,info,t=0.0_dp,products=products)
    call check(info == complex_expm_success .and. products == 0,'CSR zero time')
    call check(all(v == rhs(:,1)),'CSR zero-time identity')
    call expm_multiply_complex(a,rhs(:,1),v,info,tol=-0.1_dp)
    call check(info == complex_expm_invalid_input,'negative tolerance')
    nan_value = ieee_value(1.0_dp,ieee_quiet_nan)
    z = rhs(:,1)
    z(2) = cmplx(nan_value,0.0_dp,dp)
    call expm_multiply_complex(a,z,v,info)
    call check(info == complex_expm_invalid_input,'nonfinite right-hand side')

    ! Diagonal and zero matrix cases (including the analytical shift-only path).
    dense = cmplx(0.0_dp,0.0_dp,dp)
    do i = 1,3
        dense(i,i) = cmplx(0.2_dp*i,-0.4_dp*i,dp)
    end do
    a = csr_complex_from_dense(dense,info)
    call expm_multiply_complex(a,rhs(:,1),v,info,t=-0.8_dp)
    call check(info == complex_expm_success,'diagonal action')
    do i = 1,3
        call check(abs(v(i)-exp(-0.8_dp*dense(i,i))*rhs(i,1)) < epsilon_test,'diagonal result')
    end do
    dense = cmplx(0.0_dp,0.0_dp,dp)
    a = csr_complex_from_dense(dense,info)
    call expm_multiply_complex(a,rhs(:,1),v,info,t=-2.0_dp,products=products)
    call check(info == complex_expm_success .and. products == 0,'zero matrix no products')
    call check(all(v == rhs(:,1)),'zero matrix result')

    ! Large truly matrix-free problem: no dense or sparse matrix is allocated.
    allocate(big_diag(4096),big_b(4096))
    do i = 1,size(big_diag)
        big_diag(i) = cmplx(-0.05_dp-0.02_dp*cos(real(i,dp)), &
                            0.25_dp*sin(real(i,dp)),dp)
        big_b(i) = cmplx(cos(real(i,dp)),sin(real(i,dp)),dp)
    end do
    call expm_multiply_complex_operator(size(big_b),big_b,apply_large,0.30_dp, &
        big_y,info,t=0.7_dp,products=products)
    call check(info == complex_expm_success,'large matrix-free exponential')
    call check(maxval(abs(big_y-exp(0.7_dp*big_diag)*big_b)) < epsilon_test, &
        'large matrix-free accuracy')
    call check(products > 0,'large matrix-free matvec count')

    print '(a)', 'test_sparse_complex_expm: PASS'
contains
    function jordan_action(time) result(result)
        real(dp), intent(in) :: time
        complex(dp) :: result(3,2)
        result = exp(time*mu)*rhs
        result(1,:) = result(1,:)+exp(time*mu)*time*dense(1,2)*rhs(2,:)
    end function jordan_action

    subroutine check(ok, label)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: label
        if (.not. ok) then
            print '(a)', 'FAILED: '//label
            error stop 1
        end if
    end subroutine check

    subroutine apply_dense(x, y, status)
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: y(:)
        integer, intent(out) :: status
        y = matmul(dense,x)
        status = 0
    end subroutine apply_dense

    subroutine apply_large(x, y, status)
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: y(:)
        integer, intent(out) :: status
        y = big_diag*x
        status = 0
    end subroutine apply_large

    subroutine always_fail(x, y, status)
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: y(:)
        integer, intent(out) :: status
        y = x
        status = 19
    end subroutine always_fail

    subroutine bad_data(x, y, status)
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: y(:)
        integer, intent(out) :: status
        y = x
        y(1) = cmplx(ieee_value(1.0_dp,ieee_quiet_nan),0.0_dp,dp)
        status = 0
    end subroutine bad_data
end program test_sparse_complex_expm
