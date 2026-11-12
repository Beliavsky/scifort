! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program linalg_hadamard_leslie_dft_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_hadamard, only : hadamard
    use scifort_linalg_leslie, only : leslie
    use scifort_linalg_dft, only : dft
    implicit none

    integer :: i, status
    integer :: h(4,4)
    real(dp) :: a(4,4)
    complex(dp) :: fourier(4,4), signal(4), transformed(4)

    call hadamard(4,h,status)
    if (status /= linalg_status_success) error stop 'hadamard failed'
    print '(a)', 'Hadamard matrix (n=4):'
    do i=1,4
        print '(4(i3,1x))', h(i,:)
    end do

    call leslie([0.1_dp,2.0_dp,1.0_dp,0.1_dp], &
        [0.2_dp,0.8_dp,0.7_dp],a,status)
    if (status /= linalg_status_success) error stop 'leslie failed'
    print '(a)', 'Leslie matrix:'
    do i=1,4
        print '(4(f5.2,1x))', a(i,:)
    end do

    call dft(4,fourier,status,scale='sqrtn')
    if (status /= linalg_status_success) error stop 'dft failed'
    signal = cmplx([1.0_dp,2.0_dp,3.0_dp,4.0_dp],0.0_dp,dp)
    transformed = matmul(fourier,signal)
    print '(a)', 'Unitary four-point Fourier transform:'
    do i=1,4
        print '(2(f10.5,1x))', real(transformed(i),dp),aimag(transformed(i))
    end do
end program linalg_hadamard_leslie_dft_example
