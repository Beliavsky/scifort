! SPDX-License-Identifier: MIT
program test_linalg_real_schur
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_schur, only : schur
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none

    real(dp) :: a(4,4), t(4,4), z(4,4), r(2,2), rt(2,2), rz(2,2), near(2,2)
    real(dp) :: diagonal(3,3), dt(3,3), dz(3,3), bad(2,3), &
        one(1,1), one_t(1,1), one_z(1,1), &
        empty_a(0,0), empty_t(0,0), empty_z(0,0), big_a(41,41), big_t(41,41), big_z(41,41)
    integer :: status, sdim

    r = reshape([0.0_dp,1.0_dp,-1.0_dp,0.0_dp],[2,2])
    call schur(r,rt,rz,status)
    if (status /= linalg_status_success) error stop 'rotation real schur failed'
    call verify(r,rt,rz)
    if (abs(rt(2,1)) < 0.1_dp) error stop 'missing complex pair block'
    if (rt(1,2)*rt(2,1) >= 0.0_dp) error stop '2x2 complex block invalid'

    ! A genuinely complex pair near the real axis must not be split.
    near = reshape([1.0_dp, 1.0e-11_dp, -1.0e-11_dp, 1.0_dp], [2,2])
    call schur(near,rt,rz,status)
    if (status /= linalg_status_success) error stop 'near-real conjugate pair'
    if (abs(rt(2,1)) < 5.0e-12_dp) error stop 'near-real conjugate pair was split'
    call verify(near,rt,rz)

    a = 0.0_dp
    a(1:2,1:2) = r
    a(3,3) = -3.0_dp
    a(4,4) = 2.0_dp
    a(1,3) = 0.3_dp
    a(3,4) = 0.8_dp
    a(2,4) = -0.5_dp
    call schur(a,t,z,status)
    if (status /= linalg_status_success) error stop 'fourth-order real schur failed'
    call verify(a,t,z)
    call schur(a,t,z,status,sort='lhp',sdim=sdim)
    if (status /= linalg_status_success .or. sdim /= 1) error stop 'lhp selection failed'
    call verify(a,t,z)
    call schur(a,t,z,status,sort='iuc',sdim=sdim)
    if (status /= linalg_status_success .or. sdim /= 2) error stop 'iuc pair selection failed'
    call verify(a,t,z)
    call schur(a,t,z,status,sort='rhp',sdim=sdim)
    if (status /= linalg_status_success .or. sdim /= 3) error stop 'rhp selection failed'
    call verify(a,t,z)
    call schur(a,t,z,status,selector=positive_imag,sdim=sdim)
    if (status /= linalg_status_success .or. sdim /= 2) error stop 'pair selector failed'
    call verify(a,t,z)

    ! A skew-symmetric order-3 matrix has a true real zero eigenvalue.
    ! Its Schur eigenvalue can acquire an imaginary roundoff component;
    ! the custom selector must not count it as part of the complex pair.
    diagonal = 0.0_dp
    diagonal(1,2) = -1.0_dp
    diagonal(2,1) = 1.0_dp
    call schur(diagonal,dt,dz,status,selector=positive_imag,sdim=sdim)
    if (status /= linalg_status_success .or. sdim /= 2) error stop 'spurious selected real root'
    call verify(diagonal,dt,dz)

    diagonal = 0.0_dp
    diagonal(1,1) = 4.0_dp
    diagonal(2,2) = 4.0_dp
    diagonal(3,3) = -1.0_dp
    call schur(diagonal,dt,dz,status,sort='lhp',sdim=sdim)
    if (status /= linalg_status_success .or. sdim /= 1) error stop 'repeated real spectrum'
    call verify(diagonal,dt,dz)

    one = 0.0_dp
    call schur(one,one_t,one_z,status,sort='rhp',sdim=sdim)
    if (status /= linalg_status_success .or. sdim /= 1) error stop 'zero matrix rhp'
    call schur(empty_a,empty_t,empty_z,status)
    if (status /= linalg_status_success) error stop 'empty real Schur'
    bad = 0.0_dp
    call schur(bad,rt,rz,status)
    if (status /= linalg_status_invalid_shape) error stop 'invalid input shape accepted'
    call schur(r,one_t,rz,status)
    if (status /= linalg_status_invalid_shape) error stop 'invalid output shape accepted'
    call schur(r,rt,rz,status,sort='bogus')
    if (status /= linalg_status_invalid_shape) error stop 'invalid sort accepted'
    call schur(r,rt,rz,status,sort='lhp',selector=positive_imag)
    if (status /= linalg_status_invalid_shape) error stop 'invalid double selection accepted'
    big_a = 0.0_dp
    call schur(big_a,big_t,big_z,status)
    if (status /= linalg_status_invalid_shape) error stop 'invalid order accepted'
    r(1,1) = ieee_value(0.0_dp,ieee_quiet_nan)
    call schur(r,rt,rz,status)
    if (status /= linalg_status_no_convergence) error stop 'NaN accepted'

    print '(a)', 'Real Schur native tests passed'
contains
    logical function positive_imag(w)
        complex(dp), intent(in) :: w !! eigenvalue to select
        positive_imag = aimag(w) > 0.0_dp
    end function positive_imag

    subroutine verify(x,u,v)
        real(dp), intent(in) :: x(:, :) !! matrix to reconstruct
        real(dp), intent(in) :: u(:, :) !! real quasi-triangular factor
        real(dp), intent(in) :: v(:, :) !! orthogonal basis
        real(dp) :: identity(size(x,1),size(x,1)), scale, resid, orth
        integer :: n, k

        n = size(x,1)
        identity = 0.0_dp
        do k = 1,n
            identity(k,k) = 1.0_dp
        end do
        scale = max(1.0_dp, maxval(abs(x)))
        resid = maxval(abs(matmul(matmul(v,u),transpose(v))-x))/scale
        orth = maxval(abs(matmul(transpose(v),v)-identity))
        if (resid > 1.0e-9_dp) error stop 'real Schur reconstruction'
        if (orth > 1.0e-9_dp) error stop 'real Schur orthogonality'
        do k = 1,n-2
            if (any(abs(u(k+2:n,k)) > 1.0e-9_dp*scale)) error stop 'not upper Hessenberg'
            if (abs(u(k+1,k)) > 1.0e-9_dp*scale) then
                if (abs(u(k+2,k+1)) > 1.0e-9_dp*scale) error stop 'overlapping 2x2 blocks'
            end if
        end do
    end subroutine verify
end program test_linalg_real_schur
