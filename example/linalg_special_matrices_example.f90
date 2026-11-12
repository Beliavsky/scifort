! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program linalg_special_matrices_example
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_helmert, only : helmert
    use scifort_linalg_hilbert, only : hilbert, invhilbert, invhilbert_exact_decimal
    implicit none
    real(dp) :: contrasts(3,4), h(4,4), h_inv(4,4)
    integer(int64) :: h_inv_exact(4,4)
    character(80) :: h16_exact(16,16)
    integer :: status, i

    call helmert(4,contrasts,status)
    if (status /= linalg_status_success) error stop 'Helmert failed'
    print '(a)', 'Helmert contrast matrix, n=4:'
    do i=1,3
        print '(4(f11.7,1x))', contrasts(i,:)
    end do

    call hilbert(4,h,status)
    if (status /= linalg_status_success) error stop 'Hilbert failed'
    call invhilbert(4,h_inv,status)
    if (status /= linalg_status_success) error stop 'inverse Hilbert failed'
    call invhilbert(4,h_inv_exact,status)
    if (status /= linalg_status_success) error stop 'exact inverse Hilbert failed'
    print '(a)', 'Hilbert(4) first row:'
    print '(4(f11.7,1x))', h(1,:)
    print '(a)', 'Inverse Hilbert(4) first row, exact int64:'
    print '(4(i8,1x))', h_inv_exact(1,:)
    print '(a,es13.5)', 'Maximum |H * inv(H) - I| (n=4): ', &
        maxval(abs(matmul(h,h_inv)-identity(4)))

    call invhilbert_exact_decimal(16,h16_exact,status)
    if (status /= linalg_status_success) error stop 'arbitrary-precision inverse Hilbert failed'
    print '(a,a)', 'Exact invhilbert(16) element (8,8): ',trim(h16_exact(8,8))
contains
    pure function identity(n) result(eye)
        integer, intent(in) :: n !! square matrix order
        real(dp) :: eye(n,n) !! identity matrix
        integer :: j
        eye = 0.0_dp
        do j=1,n
            eye(j,j) = 1.0_dp
        end do
    end function identity
end program linalg_special_matrices_example
