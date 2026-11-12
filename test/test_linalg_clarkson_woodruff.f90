! SPDX-License-Identifier: MIT
! Native CountSketch tests; compile with kinds, linalg and linalg_clarkson_woodruff.
program test_linalg_clarkson_woodruff
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_clarkson_woodruff, only : countsketch_plan, make_countsketch_plan, &
        clarkson_woodruff_transform, clarkson_woodruff_csr
    implicit none

    type(countsketch_plan) :: plan, p2, bad_plan
    real(dp) :: a(4,3), expected(3,3), sketched(3,3), repeated(3,3)
    real(dp) :: batch(2,4,3), sketched_batch(2,3,3)
    real(dp) :: sparse_dense(4,3), sparse_sketched(3,3), zero_output(0,3), empty(0,3)
    real(dp) :: large(4096,6), large_sketch(768,6), norm_ratio
    complex(dp) :: ca(4,3), cs(3,3), cs_sparse(3,3)
    complex(dp) :: cbatch(2,4,3), cs_batch(2,3,3)
    integer :: row_ptr(5), col_ind(7), status, i, j, t, trials
    real(dp) :: values(7), bad_values(7)
    complex(dp) :: cvalues(7)
    real(dp), parameter :: tol = 2.0e-13_dp

    a(:,1) = [1.0_dp,2.0_dp,3.0_dp,4.0_dp]
    a(:,2) = [5.0_dp,6.0_dp,7.0_dp,8.0_dp]
    a(:,3) = [0.0_dp,9.0_dp,-2.0_dp,2.0_dp]
    plan%sketch_rows = 3
    plan%bucket = [1,3,1,2]
    plan%sign = [1,-1,-1,1]
    expected(1,:) = a(1,:)-a(3,:)
    expected(2,:) = a(4,:)
    expected(3,:) = -a(2,:)
    call clarkson_woodruff_transform(a,3,sketched,status,plan=plan)
    call verify(status == linalg_status_success, 'dense real status')
    call verify(maxval(abs(sketched - expected)) < tol, 'dense explicit CountSketch')

    ca = cmplx(a,-2.0_dp*a,dp)
    call clarkson_woodruff_transform(ca,3,cs,status,plan=plan)
    call verify(status == linalg_status_success, 'complex status')
    call verify(maxval(abs(cs-cmplx(expected,-2.0_dp*expected,dp))) < tol, &
        'complex sketch sign is real not conjugation')

    batch(1,:,:) = a
    batch(2,:,:) = 2.0_dp*a
    cbatch(1,:,:) = ca
    cbatch(2,:,:) = 2.0_dp*ca
    call clarkson_woodruff_transform(batch,3,sketched_batch,status,plan=plan)
    call verify(status == linalg_status_success, 'batch status')
    call verify(maxval(abs(sketched_batch(1,:,:) - expected)) < tol, 'batch first slice')
    call verify(maxval(abs(sketched_batch(2,:,:) - 2.0_dp*expected)) < tol, 'batch second slice')
    call clarkson_woodruff_transform(cbatch,3,cs_batch,status,plan=plan)
    call verify(status == linalg_status_success, 'complex batch status')
    call verify(maxval(abs(cs_batch(1,:,:)-cs)) < tol, 'complex batch first slice')
    call verify(maxval(abs(cs_batch(2,:,:)-2.0_dp*cs)) < tol, 'complex batch second slice')

    call make_countsketch_plan(4,3,p2,status,seed=20261109_int64)
    call verify(status == linalg_status_success, 'seeded plan status')
    call clarkson_woodruff_transform(a,3,sketched,status,seed=20261109_int64)
    call verify(status == linalg_status_success, 'seeded dense status')
    call clarkson_woodruff_transform(a,3,repeated,status,plan=p2)
    call verify(maxval(abs(sketched - repeated)) < tol, 'seeded plan/direct equivalence')
    call make_countsketch_plan(4,3,bad_plan,status,seed=20261109_int64)
    call verify(all(p2%bucket == bad_plan%bucket) .and. &
        all(p2%sign == bad_plan%sign), 'identical seeds produce identical plans')
    call make_countsketch_plan(4,3,bad_plan,status,seed=-huge(1_int64)-1_int64)
    call verify(status == linalg_status_success, 'int64 minimum seed')
    call verify(all(bad_plan%bucket >= 1) .and. all(bad_plan%bucket <= 3), 'seed range')
    call verify(all(abs(bad_plan%sign) == 1), 'seed sign range')

    ! CSR contains 7 stored values and one row with two entries.
    row_ptr = [1,3,4,6,8]
    col_ind = [1,3,2,1,2,1,3]
    values = [1.0_dp,2.0_dp,5.0_dp,7.0_dp,-2.0_dp,4.0_dp,9.0_dp]
    sparse_dense = 0.0_dp
    sparse_dense(1,1) = 1.0_dp
    sparse_dense(1,3) = 2.0_dp
    sparse_dense(2,2) = 5.0_dp
    sparse_dense(3,1) = 7.0_dp
    sparse_dense(3,2) = -2.0_dp
    sparse_dense(4,1) = 4.0_dp
    sparse_dense(4,3) = 9.0_dp
    call clarkson_woodruff_csr(row_ptr,col_ind,values,3,3,sparse_sketched,status,plan=plan)
    call verify(status == linalg_status_success, 'CSR real status')
    call clarkson_woodruff_transform(sparse_dense,3,repeated,status,plan=plan)
    call verify(maxval(abs(sparse_sketched - repeated)) < tol, 'CSR real agrees with dense')
    cvalues = cmplx(values,values/7.0_dp,dp)
    call clarkson_woodruff_csr(row_ptr,col_ind,cvalues,3,3,cs_sparse,status,plan=plan)
    call verify(status == linalg_status_success, 'CSR complex status')
    call verify(maxval(abs(cs_sparse-cmplx(repeated,repeated/7.0_dp,dp))) < tol, &
        'CSR complex agrees with dense')

    ! Unsorted or duplicated column indices and an empty CSR row are valid.
    row_ptr = [1,3,3,6,8]
    col_ind = [3,3,1,1,2,2,1]
    sparse_dense = 0.0_dp
    do i = 1, 4
        do t = row_ptr(i), row_ptr(i+1)-1
            sparse_dense(i,col_ind(t)) = sparse_dense(i,col_ind(t)) + values(t)
        end do
    end do
    call clarkson_woodruff_csr(row_ptr,col_ind,values,3,3,sparse_sketched,status,plan=plan)
    call verify(status == linalg_status_success, 'duplicate CSR status')
    call clarkson_woodruff_transform(sparse_dense,3,repeated,status,plan=plan)
    call verify(maxval(abs(sparse_sketched-repeated)) < tol, &
        'duplicate column entries and empty CSR row')
    row_ptr = [1,3,4,6,8]
    col_ind = [1,3,2,1,2,1,3]

    ! Empty core and empty sketch agree with SciPy for zero source rows.
    call clarkson_woodruff_transform(empty,0,zero_output,status)
    call verify(status == linalg_status_success, 'empty sketch')
    call clarkson_woodruff_transform(a,0,zero_output,status)
    call verify(status == linalg_status_invalid_shape, 'reject zero sketch nonempty input')
    call make_countsketch_plan(-1,3,p2,status)
    call verify(status == linalg_status_invalid_shape, 'negative input row count')
    call make_countsketch_plan(2,-1,p2,status)
    call verify(status == linalg_status_invalid_shape, 'negative sketch size')

    bad_plan = plan
    bad_plan%bucket(1) = 4
    call clarkson_woodruff_transform(a,3,sketched,status,plan=bad_plan)
    call verify(status == linalg_status_invalid_shape, 'out of bounds hash')
    call verify(maxval(abs(sketched)) == 0.0_dp, 'zero result on invalid plan')
    bad_plan = plan
    bad_plan%sign(2) = 0
    call clarkson_woodruff_transform(a,3,sketched,status,plan=bad_plan)
    call verify(status == linalg_status_invalid_shape, 'invalid sign')
    call clarkson_woodruff_transform(a,3,sketched,status,plan=plan,seed=4_int64)
    call verify(status == linalg_status_invalid_shape, 'seed and plan are exclusive')
    row_ptr(2) = 0
    call clarkson_woodruff_csr(row_ptr,col_ind,values,3,3,sketched,status,plan=plan)
    call verify(status == linalg_status_invalid_shape, 'bad CSR offsets')
    row_ptr(2) = 3
    col_ind(1) = 4
    call clarkson_woodruff_csr(row_ptr,col_ind,values,3,3,sketched,status,plan=plan)
    call verify(status == linalg_status_invalid_shape, 'bad CSR column')
    col_ind(1) = 1
    bad_values = values
    bad_values(3) = ieee_value(0.0_dp,ieee_quiet_nan)
    call clarkson_woodruff_csr(row_ptr,col_ind,bad_values,3,3,sketched,status,plan=plan)
    call verify(status == linalg_status_no_convergence, 'nonfinite CSR values')
    call verify(maxval(abs(sketched)) == 0.0_dp, 'clear nonfinite CSR output')
    a(1,1) = ieee_value(0.0_dp,ieee_quiet_nan)
    call clarkson_woodruff_transform(a,3,sketched,status,plan=plan)
    call verify(status == linalg_status_no_convergence, 'nonfinite dense input')
    call verify(maxval(abs(sketched)) == 0.0_dp, 'clear nonfinite dense output')

    ! Norm preservation is probabilistic, so average over independent seeds.
    ! A generous 10 percent threshold avoids asserting a single random outcome.
    do i = 1, size(large,1)
        do j = 1, size(large,2)
            large(i,j) = sin(real(i*j,dp)*0.113_dp) + &
                cos(real(i*(j+9),dp)*0.291_dp)
        end do
    end do
    norm_ratio = 0.0_dp
    trials = 20
    do i = 1, trials
        call clarkson_woodruff_transform(large,768,large_sketch,status,seed=int(i,int64))
        call verify(status == linalg_status_success, 'large dense sketch status')
        norm_ratio = norm_ratio + sum(large_sketch**2)/sum(large**2)
    end do
    norm_ratio = norm_ratio/real(trials,dp)
    call verify(abs(norm_ratio-1.0_dp) < 0.10_dp, 'mean squared Frobenius norm')
    print '(a, i0, a, f10.6)', 'CountSketch native tests passed; mean norm ratio (', trials, ') = ', norm_ratio

contains
    subroutine verify(condition, name)
        logical, intent(in) :: condition !! assertion to check
        character(*), intent(in) :: name !! failure message
        if (.not. condition) then
            print '(a,a)', 'FAIL: ', name
            error stop 1
        end if
    end subroutine verify
end program test_linalg_clarkson_woodruff
