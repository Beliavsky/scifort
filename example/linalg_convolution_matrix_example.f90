! SPDX-License-Identifier: MIT
program linalg_convolution_matrix_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_convolution_matrix, only : convolution_matrix, convolution_matrix_rows
    use scifort_linalg_matmul_convolution_matrix, only : matmul_convolution_matrix
    implicit none
    real(dp) :: coeff(3), x(5), a(5,5), y(5)
    complex(dp) :: cz(3), zx(5), zy(5)
    integer :: status

    coeff = [-1.0_dp,4.0_dp,-2.0_dp]
    x = [1.0_dp,2.0_dp,0.0_dp,-3.0_dp,0.5_dp]
    if (convolution_matrix_rows(size(coeff),size(x),'same') /= size(y)) &
        error stop 'expected output length five'
    call convolution_matrix(coeff,a,status,mode='same')
    if (status /= linalg_status_success) error stop 'convolution_matrix failed'
    call matmul_convolution_matrix(coeff,x,y,status,mode='same')
    if (status /= linalg_status_success) error stop 'matmul_convolution_matrix failed'
    if (maxval(abs(y-matmul(a,x))) > 1.0e-12_dp) error stop 'product mismatch'
    print '(a)', 'SciPy example, same mode:'
    print '(5f10.3)', y
    ! Expected: 2, 6, -1, -12.5, 8.

    cz = cmplx(coeff,[0.2_dp,0.1_dp,-0.3_dp],dp)
    zx = cmplx(x,[0.0_dp,-0.1_dp,1.0_dp,0.0_dp,0.2_dp],dp)
    call matmul_convolution_matrix(cz,zx,zy,status,'same')
    if (status /= linalg_status_success) error stop 'complex FFT product failed'
    print '(a,2es16.7)', 'First complex product (real, imag): ',real(zy(1),dp),aimag(zy(1))
end program linalg_convolution_matrix_example
