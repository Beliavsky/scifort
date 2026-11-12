! SPDX-License-Identifier: MIT
program linalg_banded_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_banded, only : solve_banded, solveh_banded, &
        cholesky_banded, cho_solve_banded, lu_factor_banded, lu_solve_banded
    implicit none
    real(dp) :: general(3,5), spd_upper(2,5), cholesky(2,5), b(5), x(5)
    real(dp) :: lu(4,5), second_rhs(5), matrix_rhs(5,2), matrix_x(5,2)
    integer :: ipiv(5)
    integer :: status

    general = 0.0_dp
    general(2,:) = 4.0_dp
    general(1,2:5) = -1.0_dp
    general(3,1:4) = -1.0_dp
    b = [3.0_dp,2.0_dp,2.0_dp,2.0_dp,3.0_dp]
    call solve_banded(general,1,1,b,x,status)
    if (status /= linalg_status_success) error stop 'general band solve'
    if (maxval(abs(x-1.0_dp)) > 1.0e-12_dp) error stop 'general band residual'
    print '(A,5F7.2)', 'solve_banded: ',x

    ! The LU band is larger by lower=1 rows to hold pivot-induced fill.
    call lu_factor_banded(general,1,1,lu,ipiv,status)
    if (status /= linalg_status_success) error stop 'band LU factor'
    second_rhs = 2.0_dp*b
    matrix_rhs(:,1)=b
    matrix_rhs(:,2)=second_rhs
    call lu_solve_banded(lu,1,1,ipiv,matrix_rhs,matrix_x,status)
    if (status /= linalg_status_success) error stop 'band LU reuse'
    if (maxval(abs(matrix_x(:,1)-1.0_dp)) > 1.0e-12_dp) error stop 'band LU first RHS'
    if (maxval(abs(matrix_x(:,2)-2.0_dp)) > 1.0e-12_dp) error stop 'band LU second RHS'
    print '(A,5F7.2)', 'lu_solve_banded reused: ',matrix_x(:,2)

    ! To solve A^T*x=b, specify trans='T'; trans='C' uses A^H.
    call lu_solve_banded(lu,1,1,ipiv,b,x,status,trans='T')
    if (status /= linalg_status_success) error stop 'band LU transpose'
    print '(A,5F7.2)', 'transpose solve: ',x

    ! For Hermitian systems, only the upper half band is stored.
    spd_upper = 0.0_dp
    spd_upper(2,:) = 4.0_dp
    spd_upper(1,2:5) = -1.0_dp
    call solveh_banded(spd_upper,b,x,status)
    if (status /= linalg_status_success) error stop 'SPD band solve'
    print '(A,5F7.2)', 'solveh_banded: ',x

    ! Factor once, then solve many different right-hand sides.
    call cholesky_banded(spd_upper,cholesky,status)
    if (status /= linalg_status_success) error stop 'band Cholesky'
    call cho_solve_banded(cholesky,b,x,status)
    if (status /= linalg_status_success) error stop 'band Cholesky solve'
    if (maxval(abs(x-1.0_dp)) > 1.0e-12_dp) error stop 'band Cholesky accuracy'
    print '(A,5F7.2)', 'cho_solve_banded: ',x
end program linalg_banded_example
