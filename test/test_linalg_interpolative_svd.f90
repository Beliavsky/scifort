! SPDX-License-Identifier: MIT
program test_linalg_interpolative_svd
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_interpolative_svd, only : interpolative_svd
    use scifort_linalg_interpolative, only : interp_decomp, reconstruct_matrix_from_id
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none

    real(dp), allocatable :: a(:, :), u(:, :), s(:), v(:, :), proj(:, :), rec(:, :), b(:, :)
    real(dp), allocatable :: ua(:, :), sa(:), va(:, :)
    complex(dp), allocatable :: z(:, :), uz(:, :), vz(:, :), projz(:, :), recz(:, :), bz(:, :)
    complex(dp), allocatable :: uza(:, :), vza(:, :)
    integer, allocatable :: idx(:), idxz(:)
    integer :: m, n, k, kk, j, i, variant, status, checks, mode
    real(dp) :: eps, scale, rel, nan
    real(dp), parameter :: tol = 1.0e-8_dp

    checks = 0
    do m = 0, 12
        do n = 0, 13
            allocate(a(m,n),z(m,n))
            do variant = 1, 5
                do j = 1, n
                    do i = 1, m
                        a(i,j) = sin(real(7*i+11*j+variant,dp)) + 0.01_dp*cos(real(i*j,dp))
                        z(i,j) = cmplx(a(i,j),cos(real(13*i-3*j+variant,dp)),dp)
                    end do
                end do
                select case (variant)
                case (2)
                    a = 0.0_dp
                    z = cmplx(0.0_dp,0.0_dp,dp)
                case (3)
                    do j = 2, n
                        if (m > 0) then
                            a(:,j) = a(:,1)*real(j,dp)
                            z(:,j) = z(:,1)*cmplx(real(j,dp),0.5_dp*real(j-1,dp),dp)
                        end if
                    end do
                case (4)
                    a = 1.0e-40_dp*a
                    z = 1.0e-40_dp*z
                case (5)
                    do j = 1, n
                        a(:,j) = real(j,dp)*a(:,j)
                        z(:,j) = real(j,dp)*z(:,j)
                    end do
                end select
                scale = max(1.0e-90_dp, sqrt(sum(abs(a)**2)))

                do mode = 1, 3
                    select case (mode)
                    case (1)
                        k = 0
                    case (2)
                        k = min(m,n)/2
                    case (3)
                        k = min(m,n)
                    end select
                    allocate(u(m,k),s(k),v(n,k),uz(m,k),vz(n,k),rec(m,n),recz(m,n))
                    call interpolative_svd(a,k,u,s,v,status)
                    if (status /= linalg_status_success) then
                        print *, 'real fixed rank',m,n,k,variant,status
                        error stop 'real fixed-rank SVD'
                    end if
                    call interpolative_svd(z,k,uz,s,vz,status)
                    if (status /= linalg_status_success) then
                        print *, 'complex fixed rank',m,n,k,variant,status
                        error stop 'complex fixed-rank SVD'
                    end if
                    ! Compare with the existing ID, which need not match full A.
                    call interp_decomp(a,k,idx,proj,status)
                    if (status /= linalg_status_success) error stop 'real comparison ID'
                    allocate(b(m,k))
                    b = a(:,idx(1:k))
                    call reconstruct_matrix_from_id(b,idx,proj,rec,status)
                    if (status /= linalg_status_success) error stop 'real comparison reconstruct'
                    if (k > 0) then
                        if (any(s < -tol)) error stop 'negative singular value'
                        if (any(s(2:k) > s(1:k-1)+tol*maxval(s))) error stop 'singular values unordered'
                        if (maxval(abs(matmul(transpose(u),u)-identity(k))) > tol) &
                            error stop 'real U orthonormality'
                        if (maxval(abs(matmul(transpose(v),v)-identity(k))) > tol) &
                            error stop 'real V orthonormality'
                    end if
                    ! Recompute the real decomposition because the complex call overwrote s.
                    call interpolative_svd(a,k,u,s,v,status)
                    rec = rec - matmul(u,matmul(diagonal(s),transpose(v)))
                    if (sqrt(sum(abs(rec)**2))/scale > tol) error stop 'real ID-SVD reconstruction'
                    deallocate(b,idx,proj)

                    call interp_decomp(z,k,idxz,projz,status)
                    if (status /= linalg_status_success) error stop 'complex comparison ID'
                    allocate(bz(m,k))
                    bz = z(:,idxz(1:k))
                    call reconstruct_matrix_from_id(bz,idxz,projz,recz,status)
                    if (status /= linalg_status_success) error stop 'complex comparison reconstruct'
                    if (k > 0) then
                        if (maxval(abs(matmul(transpose(conjg(uz)),uz)-identity(k))) > tol) &
                            error stop 'complex U orthonormality'
                        if (maxval(abs(matmul(transpose(conjg(vz)),vz)-identity(k))) > tol) &
                            error stop 'complex V orthonormality'
                    end if
                    call interpolative_svd(z,k,uz,s,vz,status)
                    recz = recz - matmul(uz,matmul(diagonal(s),transpose(conjg(vz))))
                    rel = sqrt(sum(abs(recz)**2))/max(1.0e-90_dp,sqrt(sum(abs(z)**2)))
                    if (rel > tol) error stop 'complex ID-SVD reconstruction'
                    checks = checks+6
                    deallocate(u,s,v,uz,vz,rec,recz,bz,idxz,projz)
                end do

                do mode = 1, 3
                    select case (mode)
                    case (1)
                        eps = 0.0_dp
                    case (2)
                        eps = 1.0e-5_dp
                    case (3)
                        eps = 0.2_dp
                    end select
                    call interpolative_svd(a,eps,k,ua,sa,va,status)
                    if (status /= linalg_status_success) error stop 'real adaptive SVD'
                    if (size(ua,1) /= m .or. size(ua,2) /= k) error stop 'real adaptive U shape'
                    if (size(va,1) /= n .or. size(va,2) /= k) error stop 'real adaptive V shape'
                    if (size(sa) /= k) error stop 'real adaptive S shape'
                    if (k > 0) then
                        if (maxval(abs(matmul(transpose(ua),ua)-identity(k))) > tol) &
                            error stop 'real adaptive U orthogonality'
                        if (maxval(abs(matmul(transpose(va),va)-identity(k))) > tol) &
                            error stop 'real adaptive V orthogonality'
                    end if
                    rel = sqrt(sum(abs(a-matmul(ua,matmul(diagonal(sa),transpose(va))))**2))/scale
                    if (rel > max(tol,eps*1.03_dp)) then
                        print *,m,n,variant,k,eps,rel
                        error stop 'real adaptive tolerance'
                    end if
                    call interpolative_svd(z,eps,kk,uza,sa,vza,status)
                    if (status /= linalg_status_success) error stop 'complex adaptive SVD'
                    if (size(uza,1) /= m .or. size(uza,2) /= kk) error stop 'complex adaptive U shape'
                    if (size(vza,1) /= n .or. size(vza,2) /= kk) error stop 'complex adaptive V shape'
                    if (size(sa) /= kk) error stop 'complex adaptive S shape'
                    if (kk > 0) then
                        if (maxval(abs(matmul(transpose(conjg(uza)),uza)-identity(kk))) > tol) &
                            error stop 'complex adaptive U orthogonality'
                        if (maxval(abs(matmul(transpose(conjg(vza)),vza)-identity(kk))) > tol) &
                            error stop 'complex adaptive V orthogonality'
                    end if
                    rel = sqrt(sum(abs(z-matmul(uza,matmul(diagonal(sa),transpose(conjg(vza)))))**2)) / &
                        max(1.0e-90_dp,sqrt(sum(abs(z)**2)))
                    if (rel > max(tol,eps*1.03_dp)) error stop 'complex adaptive tolerance'
                    checks = checks+6
                end do
            end do
            deallocate(a,z)
        end do
    end do

    ! Validate shape errors, NaNs, empty outputs, invalid tolerance.
    allocate(a(4,3),u(4,2),s(2),v(3,2))
    a = 1.0_dp
    call interpolative_svd(a,4,u,s,v,status)
    if (status /= linalg_status_invalid_shape) error stop 'invalid rank'
    call interpolative_svd(a,-1,u,s,v,status)
    if (status /= linalg_status_invalid_shape) error stop 'negative rank'
    call interpolative_svd(a,2,u(:,1:1),s,v,status)
    if (status /= linalg_status_invalid_shape) error stop 'invalid U size'
    call interpolative_svd(a,2,u,s,v(:,1:1),status)
    if (status /= linalg_status_invalid_shape) error stop 'invalid V size'
    call interpolative_svd(a,1.0_dp,k,ua,sa,va,status)
    if (status /= linalg_status_invalid_shape) error stop 'invalid eps'
    if (k /= 0) error stop 'invalid eps rank'
    if (size(ua) /= 0) error stop 'invalid eps U'
    nan = ieee_value(0.0_dp,ieee_quiet_nan)
    call interpolative_svd(a,nan,k,ua,sa,va,status)
    if (status /= linalg_status_invalid_shape) error stop 'nan eps'
    a(1,1) = nan
    call interpolative_svd(a,2,u,s,v,status)
    if (status /= linalg_status_no_convergence) error stop 'nonfinite fixed rank'
    call interpolative_svd(a,0.1_dp,k,ua,sa,va,status)
    if (status /= linalg_status_no_convergence) error stop 'nonfinite tolerance'
    checks = checks+10
    print '(a,i0)', 'interpolative SVD native checks passed: ', checks

contains

    function identity(n) result(x)
        integer, intent(in) :: n !! square identity matrix order
        real(dp) :: x(n,n)
        integer :: j
        x = 0.0_dp
        do j=1,n
            x(j,j) = 1.0_dp
        end do
    end function identity

    function diagonal(x) result(d)
        real(dp), intent(in) :: x(:) !! diagonal entries
        real(dp) :: d(size(x),size(x))
        integer :: j
        d = 0.0_dp
        do j=1,size(x)
            d(j,j) = x(j)
        end do
    end function diagonal
end program test_linalg_interpolative_svd
