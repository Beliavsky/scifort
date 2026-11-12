! SPDX-License-Identifier: MIT
program sparse_complex_expm_example
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csr_complex_from_dense
    use scifort_sparse_linalg_complex_expm, only : expm_multiply_complex, &
        expm_multiply_complex_interval, complex_expm_success
    implicit none
    type(csr_complex_matrix) :: a
    complex(dp) :: dense(2,2), b(2)
    complex(dp), allocatable :: y(:), time_series(:,:)
    integer :: info

    dense = cmplx(0.0_dp,0.0_dp,dp)
    dense(1,1) = cmplx(-0.1_dp,0.3_dp,dp)
    dense(2,2) = cmplx(-0.2_dp,0.0_dp,dp)
    dense(1,2) = cmplx(1.0_dp,0.5_dp,dp)
    b = [cmplx(1.0_dp,0.0_dp,dp),cmplx(0.0_dp,1.0_dp,dp)]
    a = csr_complex_from_dense(dense,info)
    if (info /= 0) error stop 'invalid CSR matrix'
    call expm_multiply_complex(a,b,y,info,t=0.5_dp)
    if (info /= complex_expm_success) error stop 'single-time exponential failed'
    print *, 'exp(0.5*A)*b = ', y
    call expm_multiply_complex_interval(a,b,-0.5_dp,0.5_dp,3,time_series,info)
    if (info /= complex_expm_success) error stop 'interval exponential failed'
    print *, 'middle time sample = ', time_series(:,2)
    if (maxval(abs(time_series(:,2)-b)) > 5.0e-13_dp) error stop 't=0 mismatch'
end program sparse_complex_expm_example
