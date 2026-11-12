! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program linalg_pascal_example
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_pascal, only : pascal, invpascal, pascal_exact_decimal, invpascal_exact_decimal
    implicit none
    integer(int64) :: p(5,5), pinv(5,5), lower(5,5)
    real(dp) :: floating(5,5)
    character(100) :: p50(50,50), inv50(50,50)
    integer :: status, i

    call pascal(5,p,status)
    if (status /= linalg_status_success) error stop 'Pascal failed'
    call invpascal(5,pinv,status)
    if (status /= linalg_status_success) error stop 'inverse Pascal failed'
    print '(a)', 'Symmetric Pascal matrix of order 5:'
    do i=1,5
        print '(5(i6,1x))', p(i,:)
    end do
    print '(a)', 'Its exact inverse:'
    do i=1,5
        print '(5(i6,1x))', pinv(i,:)
    end do

    call pascal(5,lower,status,kind='lower')
    if (status /= linalg_status_success) error stop 'lower Pascal failed'
    print '(a)', 'Lower Pascal matrix of order 5:'
    do i=1,5
        print '(5(i6,1x))', lower(i,:)
    end do

    call invpascal(5,floating,status,kind='upper')
    if (status /= linalg_status_success) error stop 'floating inverse Pascal failed'
    print '(a)', 'Upper inverse Pascal, floating point, first row:'
    print '(5(f6.1,1x))', floating(1,:)

    call pascal_exact_decimal(50,p50,status)
    if (status /= linalg_status_success) error stop 'arbitrary-precision Pascal failed'
    call invpascal_exact_decimal(50,inv50,status)
    if (status /= linalg_status_success) error stop 'arbitrary-precision inverse Pascal failed'
    print '(a,a)', 'Pascal(50,50), exact decimal: ', trim(p50(50,50))
    print '(a,a)', 'Inverse Pascal(1,50), exact decimal: ', trim(inv50(1,50))
end program linalg_pascal_example
