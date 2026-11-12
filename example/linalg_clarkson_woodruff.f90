! SPDX-License-Identifier: MIT
! Example: one CountSketch mapping applied to both a dense matrix and CSR data.
program example_linalg_clarkson_woodruff
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_linalg_clarkson_woodruff, only : countsketch_plan, &
        make_countsketch_plan, clarkson_woodruff_transform, clarkson_woodruff_csr
    implicit none
    type(countsketch_plan) :: plan
    real(dp) :: a(6,3), sketch(4,3), sketch_csr(4,3)
    integer :: ptr(7), col(18), status, i, j, idx
    real(dp) :: val(18)

    do i = 1, 6
        do j = 1, 3
            a(i,j) = sin(0.1_dp*real(i*j,dp))
        end do
    end do
    call make_countsketch_plan(6,4,plan,status,seed=20261109_int64)
    if (status /= 0) error stop 'CountSketch plan failed'
    call clarkson_woodruff_transform(a,4,sketch,status,plan=plan)
    if (status /= 0) error stop 'Dense CountSketch failed'

    do i = 1, 7
        ptr(i) = 3*(i-1)+1
    end do
    idx = 0
    do i = 1, 6
        do j = 1, 3
            idx = idx+1
            col(idx) = j
            val(idx) = a(i,j)
        end do
    end do
    call clarkson_woodruff_csr(ptr,col,val,3,4,sketch_csr,status,plan=plan)
    if (status /= 0) error stop 'Sparse CountSketch failed'
    if (maxval(abs(sketch_csr-sketch)) > 1.0e-14_dp) &
        error stop 'Dense/sparse CountSketch mismatch'
    print '(a)', 'Reusable CountSketch bucket/sign mapping:'
    do i = 1, 6
        print '(3(i0,1x))', i, plan%bucket(i), plan%sign(i)
    end do
    print '(a)', 'Sketch of 6-by-3 matrix (4-by-3 output):'
    do i = 1, 4
        print '(*(f12.6,1x))', sketch(i,:)
    end do
end program example_linalg_clarkson_woodruff
