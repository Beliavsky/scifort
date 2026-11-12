! SPDX-License-Identifier: MIT
program test_linalg_rq
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_rq, only : rq, rq_r
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none
    real(dp), allocatable :: a(:,:), r(:,:), q(:,:), only_r(:,:), identity(:,:), batch_a(:,:,:), &
        batch_r(:,:,:), batch_q(:,:,:)
    complex(dp), allocatable :: z(:,:), rz(:,:), qz(:,:), z_only_r(:,:), identity_z(:,:), &
        batch_z(:,:,:), batch_rz(:,:,:), batch_qz(:,:,:)
    real(dp) :: norm_a, residual, orth_err, bound, nan_value
    integer :: m, n, k, nq, i, j, mode_id, variant, status, count, b
    character(len=8) :: mode
    real(dp), parameter :: tol = 2.0e-12_dp

    count = 0
    do m = 0, 14
        do n = 0, 14
            k = min(m,n)
            allocate(a(m,n),z(m,n),only_r(m,n),z_only_r(m,n))
            do variant = 1, 5
                do j = 1, n
                    do i = 1, m
                        a(i,j) = sin(real(9*i+5*j+variant,dp))/real(i+j+1,dp)
                        z(i,j) = cmplx(a(i,j),cos(real(i+7*j+variant,dp))/real(i+2*j,dp),dp)
                    end do
                end do
                if (variant == 2 .and. n > 0) then
                    a(:,1) = 0.0_dp
                    z(:,1) = cmplx(0.0_dp,0.0_dp,dp)
                else if (variant == 3 .and. m > 1) then
                    a(m,:) = a(1,:)
                    z(m,:) = z(1,:)
                else if (variant == 4) then
                    a = 0.0_dp
                    z = cmplx(0.0_dp,0.0_dp,dp)
                else if (variant == 5) then
                    a = 1.0e-30_dp*a
                    z = cmplx(1.0e-30_dp,0.0_dp,dp)*z
                end if
                call rq_r(a,only_r,status)
                if (status /= linalg_status_success) error stop 'RQ R-only real'
                call rq_r(z,z_only_r,status)
                if (status /= linalg_status_success) error stop 'RQ R-only complex'
                do mode_id = 1, 2
                    mode = 'full'
                    nq = n
                    if (mode_id == 2) then
                        mode = 'economic'
                        nq = k
                    end if
                    allocate(r(m,nq),q(nq,n),rz(m,nq),qz(nq,n))
                    call rq(a,r,q,status,mode=trim(mode))
                    if (status /= linalg_status_success) error stop 'RQ real'
                    call rq(z,rz,qz,status,mode=trim(mode))
                    if (status /= linalg_status_success) error stop 'RQ complex'
                    norm_a = max(1.0e-90_dp, sqrt(sum(a*a)))
                    if (m*n > 0) then
                        residual = sqrt(sum((a-matmul(r,q))**2))/norm_a
                        if (residual > tol) error stop 'RQ real reconstruction'
                        norm_a = max(1.0e-90_dp, sqrt(sum(abs(z)**2)))
                        residual = sqrt(sum(abs(z-matmul(rz,qz))**2))/norm_a
                        if (residual > tol) error stop 'RQ complex reconstruction'
                    end if
                    allocate(identity(nq,nq),identity_z(nq,nq))
                    identity = 0.0_dp
                    identity_z = cmplx(0.0_dp,0.0_dp,dp)
                    do i = 1, nq
                        identity(i,i) = 1.0_dp
                        identity_z(i,i) = cmplx(1.0_dp,0.0_dp,dp)
                    end do
                    if (nq > 0) then
                        orth_err = maxval(abs(matmul(q,transpose(q))-identity))
                        if (orth_err > tol) error stop 'RQ real orthogonality'
                        orth_err = maxval(abs(matmul(qz,transpose(conjg(qz)))-identity_z))
                        if (orth_err > tol) error stop 'RQ complex unitarity'
                    end if
                    do j = 1, nq
                        do i = 1, m
                            if (j < i+nq-m) then
                                if (r(i,j) /= 0.0_dp) error stop 'RQ real triangular'
                                if (rz(i,j) /= cmplx(0.0_dp,0.0_dp,dp)) error stop 'RQ complex triangular'
                            end if
                        end do
                    end do
                    if (mode_id == 1) then
                        if (any(abs(only_r-r) > 1.0e-12_dp*max(1.0_dp,maxval(abs(r))))) &
                            error stop 'RQ real only-R mismatch'
                        if (any(abs(z_only_r-rz) > 1.0e-12_dp*max(1.0_dp,maxval(abs(rz))))) &
                            error stop 'RQ complex only-R mismatch'
                    end if
                    count = count+2
                    deallocate(r,q,rz,qz,identity,identity_z)
                end do
            end do
            deallocate(a,z,only_r,z_only_r)
        end do
    end do

    ! A batch occupies the first (leading) index as in SciPy's batch convention.
    do m = 0, 6
        do n = 0, 6
            k = min(m,n)
            allocate(batch_a(3,m,n),batch_z(3,m,n))
            do b = 1, 3
                do j = 1, n
                    do i = 1, m
                        batch_a(b,i,j) = sin(real(b+i+7*j,dp))
                        batch_z(b,i,j) = cmplx(batch_a(b,i,j),cos(real(i+3*b+j,dp)),dp)
                    end do
                end do
            end do
            do mode_id = 1,2
                mode = 'full'
                nq = n
                if (mode_id == 2) then
                    mode = 'economic'
                    nq = k
                end if
                allocate(batch_r(3,m,nq),batch_q(3,nq,n),batch_rz(3,m,nq),batch_qz(3,nq,n))
                call rq(batch_a,batch_r,batch_q,status,mode)
                if (status /= linalg_status_success) error stop 'RQ batch real'
                call rq(batch_z,batch_rz,batch_qz,status,mode)
                if (status /= linalg_status_success) error stop 'RQ batch complex'
                do b = 1, 3
                    if (m*n > 0) then
                        bound = tol*max(1.0_dp,maxval(abs(batch_a(b,:,:))))
                        if (maxval(abs(batch_a(b,:,:)-matmul(batch_r(b,:,:),batch_q(b,:,:)))) > bound) &
                            error stop 'RQ batch real reconstruction'
                        bound = tol*max(1.0_dp,maxval(abs(batch_z(b,:,:))))
                        if (maxval(abs(batch_z(b,:,:)-matmul(batch_rz(b,:,:),batch_qz(b,:,:)))) > bound) &
                            error stop 'RQ batch complex reconstruction'
                    end if
                    count = count+2
                end do
                deallocate(batch_r,batch_q,batch_rz,batch_qz)
            end do
            allocate(batch_r(3,m,n),batch_rz(3,m,n))
            call rq_r(batch_a,batch_r,status)
            if (status /= linalg_status_success) error stop 'RQ batched R-only real'
            call rq_r(batch_z,batch_rz,status)
            if (status /= linalg_status_success) error stop 'RQ batched R-only complex'
            count = count+6
            deallocate(batch_a,batch_z,batch_r,batch_rz)
        end do
    end do

    allocate(a(3,4),r(3,4),q(4,4),z(3,4),rz(3,4),qz(4,4))
    a = 1.0_dp
    z = cmplx(1.0_dp,1.0_dp,dp)
    call rq(a,r,q,status,mode='unknown')
    if (status /= linalg_status_invalid_shape) error stop 'RQ invalid mode'
    if (any(r /= 0.0_dp) .or. any(q /= 0.0_dp)) error stop 'RQ failed output'
    call rq(a,r(:,1:3),q,status)
    if (status /= linalg_status_invalid_shape) error stop 'RQ invalid output shape'
    call rq_r(a,r(:,1:3),status)
    if (status /= linalg_status_invalid_shape) error stop 'RQ invalid R shape'
    nan_value = ieee_value(1.0_dp,ieee_quiet_nan)
    a(1,1) = nan_value
    call rq(a,r,q,status)
    if (status /= linalg_status_no_convergence) error stop 'RQ NaN real'
    call rq(a,r,q,status,check_finite=.false.)
    if (status /= linalg_status_no_convergence) error stop 'RQ unchecked NaN'
    z(1,1) = cmplx(nan_value,0.0_dp,dp)
    call rq(z,rz,qz,status)
    if (status /= linalg_status_no_convergence) error stop 'RQ NaN complex'
    deallocate(a,r,q,z,rz,qz)
    allocate(batch_a(0,3,4),batch_r(0,3,4),batch_q(0,4,4))
    call rq(batch_a,batch_r,batch_q,status)
    if (status /= linalg_status_invalid_shape) error stop 'RQ zero-sized batch'

    print '(a,i0)', 'PASS linalg RQ native checks: ',count
end program test_linalg_rq
