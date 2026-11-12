! SPDX-License-Identifier: MIT
program test_linalg_cdf2rdf
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_cdf2rdf, only : cdf2rdf
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf
    implicit none

    complex(dp), allocatable :: w(:), v(:, :), wb(:, :), vb(:, :, :)
    real(dp), allocatable :: wr(:, :), vr(:, :), a(:, :), rb(:, :, :), qb(:, :, :)
    real(dp), allocatable :: rw(:), rv(:, :), rwb(:, :), rvb(:, :, :)
    real(dp) :: nan_value, inf_value, residual, scale, alpha, beta
    integer :: n, i, j, k, status, count

    count = 0
    do n = 0, 31
        allocate(w(n),v(n,n),wr(n,n),vr(n,n),a(n,n),rw(n),rv(n,n))
        do k = 1, 7
            w = cmplx(0.0_dp,0.0_dp,dp)
            v = cmplx(0.0_dp,0.0_dp,dp)
            a = 0.0_dp
            rw = 0.0_dp
            rv = 0.0_dp
            i = 1
            do while (i <= n)
                alpha = real(i+k,dp)/7.0_dp
                if (i == n .or. mod(i+k,3) == 0) then
                    w(i) = cmplx(alpha,0.0_dp,dp)
                    v(i,i) = cmplx(1.0_dp,0.0_dp,dp)
                    a(i,i) = alpha
                    i = i+1
                else
                    beta = real(i+2*k,dp)/11.0_dp
                    if (mod(k,2) == 0) beta = -beta
                    w(i) = cmplx(alpha,beta,dp)
                    w(i+1) = conjg(w(i))
                    v(i,i) = cmplx(1.0_dp,0.0_dp,dp)
                    v(i+1,i) = cmplx(0.0_dp,1.0_dp,dp)
                    v(i,i+1) = cmplx(0.0_dp,1.0_dp,dp)*conjg(v(i,i))
                    v(i+1,i+1) = cmplx(0.0_dp,1.0_dp,dp)*conjg(v(i+1,i))
                    a(i,i) = alpha
                    a(i+1,i+1) = alpha
                    a(i,i+1) = beta
                    a(i+1,i) = -beta
                    i = i+2
                end if
            end do
            ! Deliberately rotate the second eigenvector's complex phase.
            ! Its phase must not affect the real invariant basis.
            call cdf2rdf(w,v,wr,vr,status)
            if (status /= linalg_status_success) error stop 'cdf2rdf native complex status'
            if (n > 0) then
                scale = max(1.0_dp,maxval(abs(a)))
                residual = maxval(abs(matmul(a,vr)-matmul(vr,wr)))/scale
                if (residual > 2.0e-13_dp) error stop 'cdf2rdf native reconstruction'
            end if
            count = count + 1
            do j = 1, n
                rw(j) = real(j+k,dp)/5.0_dp
                rv(j,j) = 1.0_dp
            end do
            call cdf2rdf(rw,rv,wr,vr,status)
            if (status /= linalg_status_success) error stop 'cdf2rdf native real status'
            if (n > 0) then
                if (maxval(abs(vr-rv)) > 0.0_dp) error stop 'cdf2rdf real vectors'
            end if
            call cdf2rdf(rw,cmplx(rv,0.0_dp,dp),wr,vr,status)
            if (status /= linalg_status_success) error stop 'cdf2rdf real w complex v'
            count = count + 2
        end do
        deallocate(w,v,wr,vr,a,rw,rv)
    end do

    allocate(w(3),v(3,3),wr(3,3),vr(3,3),a(3,3))
    w = [cmplx(3.0_dp,5.0_dp,dp), cmplx(2.0_dp,0.0_dp,dp), &
        cmplx(3.0_dp,-5.0_dp,dp)]
    v = cmplx(0.0_dp,0.0_dp,dp)
    v(:,1) = [cmplx(0.0_dp,0.0_dp,dp),cmplx(1.0_dp,0.0_dp,dp), &
        cmplx(0.0_dp,1.0_dp,dp)]
    v(1,2) = cmplx(1.0_dp,0.0_dp,dp)
    v(:,3) = conjg(v(:,1))
    a = 0.0_dp
    a(1,1) = 2.0_dp
    a(2,2) = 3.0_dp
    a(3,3) = 3.0_dp
    a(2,3) = 5.0_dp
    a(3,2) = -5.0_dp
    call cdf2rdf(w,v,wr,vr,status)
    if (status /= linalg_status_success) error stop 'nonadjacent pair status'
    if (maxval(abs(matmul(a,vr)-matmul(vr,wr))) > 2.0e-14_dp) &
        error stop 'nonadjacent pair reconstruction'
    if (wr(1,3) /= 5.0_dp .or. wr(3,1) /= -5.0_dp) &
        error stop 'nonadjacent pair placement'
    count = count+1

    w(3) = cmplx(7.0_dp,-5.0_dp,dp)
    call cdf2rdf(w,v,wr,vr,status)
    if (status /= linalg_status_invalid_shape) error stop 'nonconjugate pair'
    if (any(wr /= 0.0_dp) .or. any(vr /= 0.0_dp)) error stop 'nonconjugate output'
    w(3) = cmplx(3.0_dp,-5.0_dp,dp)
    w(2) = cmplx(4.0_dp,1.0_dp,dp)
    call cdf2rdf(w,v,wr,vr,status)
    if (status /= linalg_status_invalid_shape) error stop 'odd complex count'
    w(2) = cmplx(2.0_dp,0.0_dp,dp)
    w(1) = cmplx(3.0_dp,5.0_dp,dp)
    w(3) = cmplx(3.0_dp,-5.0_dp,dp)
    nan_value = ieee_value(1.0_dp,ieee_quiet_nan)
    inf_value = ieee_value(1.0_dp,ieee_positive_inf)
    w(1) = cmplx(nan_value,5.0_dp,dp)
    call cdf2rdf(w,v,wr,vr,status)
    if (status /= linalg_status_no_convergence) error stop 'eigenvalue NaN'
    w(1) = cmplx(3.0_dp,5.0_dp,dp)
    v(1,2) = cmplx(inf_value,0.0_dp,dp)
    call cdf2rdf(w,v,wr,vr,status)
    if (status /= linalg_status_no_convergence) error stop 'eigenvector infinity'
    v(1,2) = cmplx(1.0_dp,0.0_dp,dp)
    call cdf2rdf(w,v,wr(1:2,:),vr,status)
    if (status /= linalg_status_invalid_shape) error stop 'bad output shape'
    w(2) = cmplx(2.0_dp,3.0e-15_dp,dp)
    call cdf2rdf(w,v,wr,vr,status)
    if (status /= linalg_status_success) error stop 'near-real eigenvalue fallback'
    if (wr(2,2) /= 2.0_dp .or. wr(2,3) /= 0.0_dp) &
        error stop 'near-real eigenvalue placement'
    w(1) = cmplx(2.0_dp,1.0e-20_dp,dp)
    w(2) = conjg(w(1))
    w(3) = cmplx(7.0_dp,0.0_dp,dp)
    call cdf2rdf(w,v,wr,vr,status)
    if (status /= linalg_status_success) error stop 'exact tiny conjugate pair'
    if (wr(1,2) /= 1.0e-20_dp) error stop 'tiny pair must remain complex'
    count = count+7
    deallocate(w,v,wr,vr,a)

    allocate(wb(3,3),vb(3,3,3),rb(3,3,3),qb(3,3,3))
    wb = cmplx(0.0_dp,0.0_dp,dp)
    vb = cmplx(0.0_dp,0.0_dp,dp)
    do i = 1, 3
        wb(i,:) = [cmplx(1.0_dp,0.0_dp,dp),cmplx(real(i,dp),2.0_dp,dp), &
            cmplx(real(i,dp),-2.0_dp,dp)]
        vb(i,1,1) = cmplx(1.0_dp,0.0_dp,dp)
        vb(i,2,2) = cmplx(1.0_dp,0.0_dp,dp)
        vb(i,3,2) = cmplx(0.0_dp,1.0_dp,dp)
        vb(i,2,3) = cmplx(1.0_dp,0.0_dp,dp)
        vb(i,3,3) = cmplx(0.0_dp,-1.0_dp,dp)
    end do
    call cdf2rdf(wb,vb,rb,qb,status)
    if (status /= linalg_status_success) error stop 'complex batch status'
    if (maxval(abs(rb(:,2,3)-2.0_dp)) > 0.0_dp) error stop 'complex batch blocks'
    wb(2,3) = cmplx(2.0_dp,-1.0_dp,dp)
    call cdf2rdf(wb,vb,rb,qb,status)
    if (status /= linalg_status_invalid_shape) error stop 'complex batch bad pair'
    if (any(rb /= 0.0_dp) .or. any(qb /= 0.0_dp)) error stop 'batch failure cleanup'
    deallocate(wb,vb,rb,qb)

    allocate(rwb(2,3),rvb(2,3,3),rb(2,3,3),qb(2,3,3))
    rwb = 2.0_dp
    rvb = 0.0_dp
    do i = 1, 3
        rvb(:,i,i) = 1.0_dp
    end do
    call cdf2rdf(rwb,rvb,rb,qb,status)
    if (status /= linalg_status_success) error stop 'real batch status'
    if (maxval(abs(qb-rvb)) > 0.0_dp) error stop 'real batch values'
    call cdf2rdf(rwb,cmplx(rvb,0.0_dp,dp),rb,qb,status)
    if (status /= linalg_status_success) error stop 'real w complex v batch'
    count = count+4

    print '(a,i0)', 'PASS cdf2rdf native cases: ',count
end program test_linalg_cdf2rdf
