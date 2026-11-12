! SPDX-License-Identifier: MIT
program test_linalg_interpolative
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, &
        linalg_status_invalid_shape, linalg_status_no_convergence
    use scifort_linalg_interpolative, only : interp_decomp, reconstruct_skel_matrix, &
        reconstruct_interp_matrix, reconstruct_matrix_from_id, id_to_svd
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none
    real(dp), allocatable :: a(:, :), b(:, :), p(:, :), rec(:, :), proj(:, :)
    real(dp), allocatable :: u(:, :), sv(:), v(:, :), t(:, :)
    complex(dp), allocatable :: z(:, :), bz(:, :), pz(:, :), rz(:, :), projz(:, :)
    complex(dp), allocatable :: uz(:, :), vz(:, :), tz(:, :)
    integer, allocatable :: idx(:), idxz(:)
    integer :: m, n, k, kk, status, i, j, variant, checks
    real(dp) :: scale, err
    real(dp), parameter :: tol = 4.0e-9_dp

    checks = 0
    do m = 0, 12
        do n = 0, 14
            allocate(a(m,n),z(m,n))
            do variant = 1, 5
                do j = 1, n
                    do i = 1, m
                        a(i,j) = sin(real(i*3+j*7+variant,dp)) + 0.2_dp*cos(real(i*j,dp))
                        z(i,j) = cmplx(a(i,j), sin(real(i*11-j*4+variant,dp)), dp)
                    end do
                end do
                select case (variant)
                case (2)
                    a = 0.0_dp
                    z = cmplx(0.0_dp,0.0_dp,dp)
                case (3)
                    if (m > 0 .and. n > 0) then
                        do j = 1, n
                            a(:,j) = a(:,1)*real(j,dp)
                            z(:,j) = z(:,1)*cmplx(real(j,dp),0.3_dp*real(j-1,dp),dp)
                        end do
                    end if
                case (4)
                    a = a * 1.0e-45_dp
                    z = z * 1.0e-45_dp
                case (5)
                    do j = 1, n
                        a(:,j) = a(:,j)*real(j,dp)
                        z(:,j) = z(:,j)*real(j,dp)
                    end do
                end select

                ! Fixed rank ID, including k=0 and k=min(m,n).
                k = min(m,n)/2
                if (variant == 2) k = min(m,n)
                if (variant == 5) k = min(m,n)
                call interp_decomp(a,k,idx,proj,status)
                if (status /= linalg_status_success) then
                    print *, 'failed real rank', m,n,k,variant,status
                    error stop 'real fixed rank ID'
                end if
                call interp_decomp(z,k,idxz,projz,status)
                if (status /= linalg_status_success) then
                    print *, 'failed complex rank', m,n,k,variant,status
                    error stop 'complex fixed rank ID'
                end if
                allocate(b(m,k),p(k,n),rec(m,n),bz(m,k),pz(k,n),rz(m,n))
                call reconstruct_skel_matrix(a,k,idx,b,status)
                if (status /= linalg_status_success) error stop 'real skeleton'
                call reconstruct_interp_matrix(idx,proj,p,status)
                if (status /= linalg_status_success) error stop 'real interpolation matrix'
                call reconstruct_matrix_from_id(b,idx,proj,rec,status)
                if (status /= linalg_status_success) error stop 'real matrix reconstruction'
                if (m*n > 0) then
                    if (maxval(abs(rec-matmul(b,p))) > tol*max(1.0e-90_dp,maxval(abs(a)))) &
                        error stop 'real reconstruction identity'
                    if (maxval(abs(b-a(:,idx(1:k)))) > 0.0_dp) error stop 'real original skeleton'
                end if
                call reconstruct_skel_matrix(z,k,idxz,bz,status)
                if (status /= linalg_status_success) error stop 'complex skeleton'
                call reconstruct_interp_matrix(idxz,projz,pz,status)
                if (status /= linalg_status_success) error stop 'complex interpolation matrix'
                call reconstruct_matrix_from_id(bz,idxz,projz,rz,status)
                if (status /= linalg_status_success) error stop 'complex reconstruction'
                if (m*n > 0) then
                    if (maxval(abs(rz-matmul(bz,pz))) > tol*max(1.0e-90_dp,maxval(abs(z)))) &
                        error stop 'complex reconstruction identity'
                    if (maxval(abs(bz-z(:,idxz(1:k)))) > 0.0_dp) error stop 'complex original skeleton'
                end if
                checks = checks+8
                if (variant == 2 .or. variant == 5) then
                    if (m*n > 0) then
                        scale = max(1.0e-90_dp,maxval(abs(a)))
                        if (maxval(abs(rec-a))/scale > tol) error stop 'real full-rank recovery'
                        scale = max(1.0e-90_dp,maxval(abs(z)))
                        if (maxval(abs(rz-z))/scale > tol) error stop 'complex full-rank recovery'
                    end if
                    checks = checks+2
                end if
                if (k > 0 .and. mod(m+n,3) == 0) then
                    allocate(u(m,k),sv(k),v(n,k),t(m,n),uz(m,k),vz(n,k),tz(m,n))
                    call id_to_svd(b,idx,proj,u,sv,v,status)
                    if (status /= linalg_status_success) error stop 'real ID to SVD'
                    t = matmul(u,matmul(diagonal_real(sv),transpose(v)))
                    scale = max(1.0e-90_dp,maxval(abs(rec)))
                    err = maxval(abs(t-rec))/scale
                    if (err > tol) then
                        print *, m,n,k,err
                        error stop 'real ID-to-SVD reconstruction'
                    end if
                    call id_to_svd(bz,idxz,projz,uz,sv,vz,status)
                    if (status /= linalg_status_success) error stop 'complex ID to SVD'
                    tz = matmul(uz,matmul(diagonal_real(sv),transpose(conjg(vz))))
                    scale = max(1.0e-90_dp,maxval(abs(rz)))
                    err = maxval(abs(tz-rz))/scale
                    if (err > tol) then
                        print *, m,n,k,err
                        error stop 'complex ID-to-SVD reconstruction'
                    end if
                    checks = checks+2
                    deallocate(u,sv,v,t,uz,vz,tz)
                end if
                deallocate(idx,proj,idxz,projz,b,p,rec,bz,pz,rz)

                ! Tolerance-based variant; exact rank-one inputs remain accurate.
                call interp_decomp(a,1.0e-8_dp,kk,idx,proj,status)
                if (status /= linalg_status_success) error stop 'real tolerance ID'
                allocate(b(m,kk),rec(m,n))
                call reconstruct_skel_matrix(a,kk,idx,b,status)
                if (status /= linalg_status_success) error stop 'real tolerance skeleton'
                call reconstruct_matrix_from_id(b,idx,proj,rec,status)
                if (status /= linalg_status_success) error stop 'real tolerance reconstruction'
                if (m*n > 0) then
                    scale = max(1.0e-90_dp,sqrt(sum((a/max(1.0e-90_dp,maxval(abs(a))))**2)))
                    ! Ratio computed with the same scale to avoid overflow/underflow.
                    if (maxval(abs(a)) > 0.0_dp) then
                        err = sqrt(sum(((rec-a)/maxval(abs(a)))**2))/scale
                        if (err > 3.0e-8_dp) then
                            print *, 'real tol',m,n,variant,kk,err
                            error stop 'real relative Frobenius tolerance'
                        end if
                    end if
                end if
                deallocate(idx,proj,b,rec)
                call interp_decomp(z,1.0e-8_dp,kk,idxz,projz,status)
                if (status /= linalg_status_success) error stop 'complex tolerance ID'
                allocate(bz(m,kk),rz(m,n))
                call reconstruct_skel_matrix(z,kk,idxz,bz,status)
                if (status /= linalg_status_success) error stop 'complex tolerance skeleton'
                call reconstruct_matrix_from_id(bz,idxz,projz,rz,status)
                if (status /= linalg_status_success) error stop 'complex tolerance reconstruction'
                if (m*n > 0) then
                    if (maxval(abs(z)) > 0.0_dp) then
                        scale = maxval(abs(z))
                        err = sqrt(sum((abs((rz-z)/scale))**2))/sqrt(sum((abs(z/scale))**2))
                        if (err > 3.0e-8_dp) then
                            print *, 'complex tol',m,n,variant,kk,err
                            error stop 'complex relative Frobenius tolerance'
                        end if
                    end if
                end if
                checks = checks+6
                deallocate(idxz,projz,bz,rz)
            end do
            deallocate(a,z)
        end do
    end do

    ! Invalid options and permutations must return a status, not touch invalid memory.
    allocate(a(3,3),z(3,3),b(3,2),rec(3,3),p(2,3),bz(3,2),rz(3,3),pz(2,3))
    a = 0.0_dp
    a(1,1) = 1.0_dp
    a(2,2) = 2.0_dp
    a(3,3) = 3.0_dp
    z = cmplx(a,0.0_dp,dp)
    call interp_decomp(a,4,idx,proj,status)
    if (status /= linalg_status_invalid_shape) error stop 'rank outside range'
    call interp_decomp(a,-0.01_dp,kk,idx,proj,status)
    if (status /= linalg_status_invalid_shape) error stop 'negative tolerance'
    call interp_decomp(a,1.0_dp,kk,idx,proj,status)
    if (status /= linalg_status_invalid_shape) error stop 'invalid tolerance'
    call interp_decomp(z,2,idxz,projz,status)
    if (status /= linalg_status_success) error stop 'error test decomposition'
    idxz = [1,1,3]
    call reconstruct_skel_matrix(z,2,idxz,bz,status)
    if (status /= linalg_status_invalid_shape) error stop 'duplicate permutation for skeleton'
    call reconstruct_interp_matrix(idxz,projz,pz,status)
    if (status /= linalg_status_invalid_shape) error stop 'duplicate permutation for P'
    call reconstruct_matrix_from_id(bz,idxz,projz,rz,status)
    if (status /= linalg_status_invalid_shape) error stop 'duplicate permutation for A'
    a(1,2) = ieee_value(1.0_dp,ieee_quiet_nan)
    call interp_decomp(a,1,idx,proj,status)
    if (status /= linalg_status_no_convergence) error stop 'nonfinite real input'
    z(1,2) = cmplx(ieee_value(1.0_dp,ieee_quiet_nan),0.0_dp,dp)
    call interp_decomp(z,0.5_dp,kk,idxz,projz,status)
    if (status /= linalg_status_no_convergence) error stop 'nonfinite complex input'
    checks = checks+8
    print '(a,i0)', 'PASS interpolative native cases: ',checks

contains
    pure function diagonal_real(s) result(d)
        real(dp), intent(in) :: s(:) !! singular values to place on diagonal
        real(dp) :: d(size(s),size(s))
        integer :: j
        d = 0.0_dp
        do j = 1, size(s)
            d(j,j) = s(j)
        end do
    end function diagonal_real
end program test_linalg_interpolative
