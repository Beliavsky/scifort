! SPDX-License-Identifier: MIT
! Native tests for generalized eigenvalues/eigenvectors, including beta=0.
program test_linalg_generalized_eigen
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence, linalg_status_singular
    use scifort_linalg_generalized_eigen, only : eig, eigvals
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    implicit none
    real(dp) :: a(3,3), b(3,3), rotation(2,2), mass(2,2), invalid(3,2)
    complex(dp) :: w(3), alpha(3), beta(3), vr(3,3), vl(3,3)
    complex(dp) :: w2(2), u2(2,2), v2(2,2), al2(2), be2(2), cb(2,2)
    complex(dp) :: bad_w(2), short(1), emptyw(0), tiny_w(1), tiny_alpha(1), tiny_beta(1)
    complex(dp) :: emptya(0,0), emptyb(0,0), emptyvec(0,0)
    real(dp) :: residual, norme, tiny_a(1,1), tiny_b(1,1)
    integer :: status, i

    a = 0.0_dp
    b = 0.0_dp
    a(1,1) = 2.0_dp
    a(2,2) = 6.0_dp
    a(3,3) = 12.0_dp
    b(1,1) = 1.0_dp
    b(2,2) = 2.0_dp
    b(3,3) = 3.0_dp
    call eig(a,b,w,status,vr=vr,vl=vl,alpha=alpha,beta=beta,method='iterative')
    if (status /= linalg_status_success) error stop 'finite: status'
    if (maxval(abs(w-[cmplx(2.0_dp,0.0_dp,dp),cmplx(3.0_dp,0.0_dp,dp), &
        cmplx(4.0_dp,0.0_dp,dp)])) > 1.0e-10_dp) then
        ! QZ does not promise a particular eigenvalue order.
        if (.not. contains_roots(w,[2.0_dp,3.0_dp,4.0_dp])) error stop 'finite eigenvalues'
    end if
    do i=1,3
        residual = maxval(abs(beta(i)*matmul(cmplx(a,kind=dp),vr(:,i)) - &
            alpha(i)*matmul(cmplx(b,kind=dp),vr(:,i))))
        if (residual > 1.0e-9_dp) error stop 'finite right eigenvectors'
        residual = maxval(abs(conjg(beta(i))*matmul(transpose(cmplx(a,kind=dp)),vl(:,i)) - &
            conjg(alpha(i))*matmul(transpose(cmplx(b,kind=dp)),vl(:,i))))
        if (residual > 1.0e-9_dp) error stop 'finite left eigenvectors'
        norme = sqrt(sum(abs(vr(:,i))**2))
        if (abs(norme-1.0_dp) > 1.0e-12_dp) error stop 'right normalization'
    end do
    call eigvals(a,b,w,status,alpha,beta,method='inverse')
    if (status /= linalg_status_success) error stop 'inverse eigvals'
    if (.not. contains_roots(w,[2.0_dp,3.0_dp,4.0_dp])) error stop 'inverse spectrum'

    ! Both A and B singular, regular pencil with 0 and infinity.
    a(3,3) = 1.0_dp
    a(2,2) = 0.0_dp
    b(1,1) = 0.0_dp
    b(2,2) = 3.0_dp
    b(3,3) = 1.0_dp
    call eig(a,b,w,status,vr=vr,vl=vl,alpha=alpha,beta=beta)
    if (status /= linalg_status_success) error stop 'regular singular: status'
    if (.not. any(.not. ieee_is_finite(real(w,dp)))) error stop 'infinite root missing'
    if (.not. any(abs(w) < 1.0e-12_dp)) error stop 'zero root missing'
    do i=1,3
        residual = maxval(abs(beta(i)*matmul(cmplx(a,kind=dp),vr(:,i)) - &
            alpha(i)*matmul(cmplx(b,kind=dp),vr(:,i))))
        if (residual > 1.0e-8_dp) error stop 'singular pencil right residual'
    end do

    rotation = reshape([0.0_dp,1.0_dp,-1.0_dp,0.0_dp],[2,2])
    mass = 0.0_dp
    mass(1,1) = 2.0_dp
    mass(2,2) = 2.0_dp
    cb = cmplx(mass,kind=dp)
    call eig(rotation,cb,w2,status,vr=u2,vl=v2,alpha=al2,beta=be2)
    if (status /= linalg_status_success) error stop 'real-complex overload'
    if (maxval(abs(abs(w2)-0.5_dp)) > 1.0e-10_dp) error stop 'rotation roots'
    call eig(cb,rotation,w2,status,vr=u2,vl=v2)
    if (status /= linalg_status_success) error stop 'complex-real overload'
    call eig(cb,cb,w2,status,vr=u2,vl=v2)
    if (status /= linalg_status_success) error stop 'complex-complex overload'

    call eig(emptya,emptyb,emptyw,status,vr=emptyvec,vl=emptyvec)
    if (status /= linalg_status_success) error stop 'empty case'
    invalid = 0.0_dp
    call eigvals(invalid,b,w,status)
    if (status /= linalg_status_invalid_shape) error stop 'shape accepted'
    call eigvals(rotation,mass,short,status)
    if (status /= linalg_status_invalid_shape) error stop 'short output accepted'
    call eigvals(rotation,mass,bad_w,status,method='unknown')
    if (status /= linalg_status_invalid_shape) error stop 'unknown method accepted'
    mass(1,1) = ieee_value(0.0_dp,ieee_quiet_nan)
    call eigvals(rotation,mass,bad_w,status)
    if (status /= linalg_status_no_convergence) error stop 'NaN accepted'
    ! Physical eigenvalues can exceed the finite range; homogeneous
    ! coordinates still allow a meaningful result without overflow traps.
    tiny_a(1,1) = 1.0e200_dp
    tiny_b(1,1) = 1.0e-200_dp
    call eigvals(tiny_a,tiny_b,tiny_w,status,alpha=tiny_alpha,beta=tiny_beta)
    if (status /= linalg_status_success) error stop 'huge ratio: status'
    if (ieee_is_finite(real(tiny_w(1),dp))) error stop 'huge ratio: infinity'
    if (abs(tiny_alpha(1)) <= 0.0_dp .or. abs(tiny_beta(1)) <= 0.0_dp) &
        error stop 'huge ratio: homogeneous representation'
    tiny_a(1,1) = 1.0e-200_dp
    tiny_b(1,1) = 1.0e200_dp
    call eigvals(tiny_a,tiny_b,tiny_w,status)
    if (status /= linalg_status_success) error stop 'tiny ratio: status'
    if (abs(tiny_w(1)) > 1.0e-200_dp) error stop 'tiny ratio: underflow'
    a = 0.0_dp
    b = 0.0_dp
    call eigvals(a,b,w,status)
    if (status /= linalg_status_singular) error stop 'singular pencil accepted'
    print '(a)', 'generalized eig native tests passed'
contains
    logical function contains_roots(roots, targets) result(ok)
        complex(dp), intent(in) :: roots(:) !! generalized eigenvalues
        real(dp), intent(in) :: targets(:) !! expected real eigenvalues
        integer :: j
        ok = .true.
        do j=1,size(targets)
            if (.not. any(abs(roots-targets(j)) < 1.0e-10_dp)) ok = .false.
        end do
    end function contains_roots
end program test_linalg_generalized_eigen
