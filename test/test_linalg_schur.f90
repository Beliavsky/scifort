! SPDX-License-Identifier: MIT
program test_linalg_schur
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_schur, only : schur
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none

    real(dp) :: a(4,4), rotation(2,2), malformed(2,3), big(41,41)
    complex(dp) :: c(4,4), t(4,4), z(4,4), rt(2,2), rz(2,2), badt(2,2), &
        too_t(41,41), too_z(41,41), empty_a(0,0), empty_t(0,0), empty_z(0,0), &
        one(1,1), one_t(1,1), one_z(1,1)
    integer :: code, nselected

    a = 0.0_dp
    a(1,1) = -4.0_dp
    a(2,2) = 2.0_dp
    a(3,3) = -0.8_dp
    a(4,4) = 3.0_dp
    a(1,2) = 7.0_dp
    a(2,3) = -3.0_dp
    a(3,4) = 0.75_dp
    call schur(a,t,z,code,sort='lhp',sdim=nselected)
    if (code /= linalg_status_success .or. nselected /= 2) error stop 'lhp sort'
    if (any(real([t(1,1),t(2,2)],dp) >= 0.0_dp)) error stop 'lhp selection'
    if (any(real([t(3,3),t(4,4)],dp) < 0.0_dp)) error stop 'lhp unselected'
    call verify(cmplx(a,kind=dp),t,z)
    call schur(a,t,z,code,sort='rhp',sdim=nselected)
    if (code /= linalg_status_success .or. nselected /= 2) error stop 'rhp sort'
    call verify(cmplx(a,kind=dp),t,z)
    call schur(a,t,z,code,sort='iuc',sdim=nselected)
    if (code /= linalg_status_success .or. nselected /= 1) error stop 'iuc sort'
    call verify(cmplx(a,kind=dp),t,z)
    call schur(a,t,z,code,sort='ouc',sdim=nselected)
    if (code /= linalg_status_success .or. nselected /= 3) error stop 'ouc sort'
    call verify(cmplx(a,kind=dp),t,z)
    call schur(a,t,z,code,selector=positive_imag,sdim=nselected)
    if (code /= linalg_status_success .or. nselected /= 0) error stop 'user selector'

    rotation = reshape([0.0_dp,1.0_dp,-1.0_dp,0.0_dp],[2,2])
    call schur(rotation,rt,rz,code,selector=positive_imag,sdim=nselected)
    if (code /= linalg_status_success .or. nselected /= 1) error stop 'rotation selector'
    if (aimag(rt(1,1)) <= 0.0_dp) error stop 'selected rotation eigenvalue'
    call verify(cmplx(rotation,kind=dp),rt,rz)

    c = cmplx(a,kind=dp)
    c(2,1) = cmplx(0.5_dp,1.0_dp,dp)
    c(4,2) = cmplx(-0.3_dp,0.4_dp,dp)
    call schur(c,t,z,code)
    if (code /= linalg_status_success) error stop 'complex unsorted'
    call verify(c,t,z)
    call schur(c,t,z,code,sort='lhp',sdim=nselected)
    if (code /= linalg_status_success .or. nselected /= 2) error stop 'complex sorted'
    call verify(c,t,z)

    ! Edge cases and invalid inputs.
    call schur(empty_a,empty_t,empty_z,code)
    if (code /= linalg_status_success) error stop 'empty'
    one = cmplx(0.0_dp,0.0_dp,dp)
    one(1,1) = cmplx(-2.0_dp,3.0_dp,dp)
    call schur(one,one_t,one_z,code,sort='lhp',sdim=nselected)
    if (code /= linalg_status_success .or. nselected /= 1) error stop 'order 1'
    one(1,1) = cmplx(0.0_dp,0.0_dp,dp)
    call schur(one,one_t,one_z,code,sort='rhp',sdim=nselected)
    if (code /= linalg_status_success .or. nselected /= 1) error stop 'zero rhp'
    call schur(one,one_t,one_z,code,sort='ouc',sdim=nselected)
    if (code /= linalg_status_success .or. nselected /= 0) error stop 'zero ouc'
    malformed = 0.0_dp
    call schur(malformed,rt,rz,code)
    if (code /= linalg_status_invalid_shape) error stop 'nonsquare accepted'
    call schur(rotation,one_t,rz,code)
    if (code /= linalg_status_invalid_shape) error stop 'wrong output shape accepted'
    call schur(rotation,badt,rz,code,sort='bogus')
    if (code /= linalg_status_invalid_shape) error stop 'bad sort accepted'
    call schur(rotation,badt,rz,code,sort='lhp',selector=positive_imag)
    if (code /= linalg_status_invalid_shape) error stop 'multiple selectors accepted'
    big = 0.0_dp
    call schur(big,too_t,too_z,code)
    if (code /= linalg_status_invalid_shape) error stop 'order bound accepted'
    rotation(1,1) = ieee_value(0.0_dp,ieee_quiet_nan)
    call schur(rotation,rt,rz,code)
    if (code /= linalg_status_no_convergence) error stop 'nonfinite accepted'

    print '(a)', 'Schur native tests passed'

contains

    logical function positive_imag(eigenvalue)
        complex(dp), intent(in) :: eigenvalue !! candidate complex eigenvalue
        positive_imag = aimag(eigenvalue) > 0.0_dp
    end function positive_imag

    subroutine verify(matrix,t,z)
        complex(dp), intent(in) :: matrix(:, :) !! input matrix
        complex(dp), intent(in) :: t(:, :) !! upper-triangular Schur factor
        complex(dp), intent(in) :: z(:, :) !! unitary transformation
        complex(dp) :: identity(size(t,1),size(t,1))
        real(dp) :: err, unitary, triangular
        integer :: k, n

        n = size(t,1)
        identity = cmplx(0.0_dp,0.0_dp,dp)
        do k = 1, n
            identity(k,k) = cmplx(1.0_dp,0.0_dp,dp)
        end do
        err = maxval(abs(matmul(matmul(z,t),transpose(conjg(z)))-matrix))
        unitary = maxval(abs(matmul(transpose(conjg(z)),z)-identity))
        triangular = 0.0_dp
        do k = 1, n-1
            triangular = max(triangular,maxval(abs(t(k+1:n,k))))
        end do
        if (err > 1e-10_dp*max(1.0_dp,maxval(abs(matrix)))) error stop 'Schur reconstruction'
        if (unitary > 1e-10_dp) error stop 'Schur unitarity'
        if (triangular > 1e-10_dp*max(1.0_dp,maxval(abs(t)))) error stop 'Schur triangularity'
    end subroutine verify

end program test_linalg_schur
