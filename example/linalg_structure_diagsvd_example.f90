! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program linalg_structure_diagsvd_example
    use scifort_kinds, only : dp
    use scifort_linalg_structure, only : bandwidth, issymmetric, ishermitian
    use scifort_linalg_diagsvd, only : diagsvd
    implicit none
    real(dp) :: a(4,4), sigmas(3), rectangular_sigma(3,5)
    complex(dp) :: h(2,2)
    integer :: lower, upper, status, i
    logical :: symmetric, hermitian

    a = 0.0_dp
    do i = 1, 4
        a(i,i) = real(i,dp)
    end do
    a(4,1) = 2.0_dp
    a(1,2) = 3.0_dp
    call bandwidth(a, lower, upper, status)
    if (status /= 0) error stop 'bandwidth failed'
    print '(a,2(1x,i0))', 'lower/upper bandwidth:', lower, upper

    h(1,1) = cmplx(1.0_dp,0.0_dp,dp)
    h(2,2) = cmplx(2.0_dp,0.0_dp,dp)
    h(1,2) = cmplx(3.0_dp,4.0_dp,dp)
    h(2,1) = conjg(h(1,2))
    call ishermitian(h, hermitian, status)
    if (status /= 0) error stop 'ishermitian failed'
    call issymmetric(h, symmetric, status)
    if (status /= 0) error stop 'issymmetric failed'
    print '(a,l1)', 'Hermitian? ', hermitian
    print '(a,l1)', 'Symmetric? ', symmetric

    sigmas = [4.0_dp, 2.0_dp, 1.0_dp]
    call diagsvd(sigmas, 3, 5, rectangular_sigma, status)
    if (status /= 0) error stop 'diagsvd failed'
    print *, 'rectangular SVD sigma matrix (3x5):'
    do i = 1, 3
        print '(5f6.1)', rectangular_sigma(i,:)
    end do
end program linalg_structure_diagsvd_example
