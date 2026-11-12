! SPDX-License-Identifier: MIT
! Native regression checks of randomized interpolative SVD.
program test_linalg_interpolative_svd_random
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_interpolative_svd_random, only : randomized_interpolative_svd
    implicit none

    real(dp), allocatable :: a(:, :), u(:, :), v(:, :), s(:), rec(:, :)
    real(dp), allocatable :: ua(:, :), va(:, :), sa(:), previous(:)
    complex(dp), allocatable :: z(:, :), uz(:, :), vz(:, :), rz(:, :)
    complex(dp), allocatable :: uaz(:, :), vaz(:, :)
    integer :: m, n, k, rank, j, i, mode, family, status, checks
    real(dp) :: scale, rel, tolerance, xnan
    integer(int64), parameter :: seed = 314159_int64
    real(dp), parameter :: tol = 5.0e-8_dp

    checks = 0
    do m = 0, 13
        do n = 0, 14
            allocate(a(m,n),z(m,n))
            do family = 1, 5
                do j = 1, n
                    do i = 1, m
                        a(i,j) = sin(real(11*i+3*j+family,dp)) + cos(real(7*i-9*j,dp))*0.14_dp
                        z(i,j) = cmplx(a(i,j),cos(real(13*i+5*j,dp)),dp)
                    end do
                end do
                select case (family)
                case (2)
                    a = 0.0_dp
                    z = cmplx(0.0_dp,0.0_dp,dp)
                case (3)
                    do j = 2, n
                        if (m == 0) exit
                        a(:,j) = a(:,1)*real(j,dp)
                        z(:,j) = z(:,1)*cmplx(real(j,dp),0.15_dp,dp)
                    end do
                case (4)
                    a = a*1.0e-40_dp
                    z = z*1.0e-40_dp
                case (5)
                    a = a*1.0e40_dp
                    z = z*1.0e40_dp
                end select
                scale = 1.0_dp
                if (m > 0 .and. n > 0) scale = max(1.0e-50_dp,maxval(abs(a)))

                do mode = 0, 2
                    k = 0
                    if (mode == 1) k = min(m,n)/2
                    if (mode == 2) k = min(m,n)
                    allocate(u(m,k),v(n,k),s(k),uz(m,k),vz(n,k),rec(m,n),rz(m,n))
                    call randomized_interpolative_svd(a,k,u,s,v,status,seed=seed)
                    if (status /= linalg_status_success) then
                        print *, 'real fixed error:',m,n,family,k,status
                        error stop 'real randomized SVD'
                    end if
                    if (k > 0) then
                        if (any(s < 0.0_dp)) error stop 'negative singular values'
                        if (any(s(2:k) > s(1:k-1)*(1.0_dp+1.0e-12_dp))) &
                            error stop 'unordered real singular values'
                        if (maxval(abs(matmul(transpose(u),u)-eye(k))) > tol) error stop 'real U not orthogonal'
                        if (maxval(abs(matmul(transpose(v),v)-eye(k))) > tol) error stop 'real V not orthogonal'
                    end if
                    if (mode == 2) then
                        rec = a-matmul(u,matmul(diagonal(s),transpose(v)))
                        rel = sqrt(sum((rec/scale)**2)) / max(1.0e-40_dp,sqrt(sum((a/scale)**2)))
                        if (rel > tol) then
                            print *, 'real full-rank residual:',m,n,family,rel
                            error stop 'real full rank'
                        end if
                    end if
                    call randomized_interpolative_svd(z,k,uz,s,vz,status,seed=seed)
                    if (status /= linalg_status_success) then
                        print *, 'complex fixed error:',m,n,family,k,status
                        error stop 'complex randomized SVD'
                    end if
                    if (k > 0) then
                        if (maxval(abs(matmul(transpose(conjg(uz)),uz)-eye(k))) > tol) &
                            error stop 'complex U not unitary'
                        if (maxval(abs(matmul(transpose(conjg(vz)),vz)-eye(k))) > tol) &
                            error stop 'complex V not unitary'
                    end if
                    if (mode == 2) then
                        rz = z-matmul(uz,matmul(diagonal(s),transpose(conjg(vz))))
                        rel = sqrt(sum(abs(rz/scale)**2)) / &
                            max(1.0e-40_dp,sqrt(sum(abs(z/scale)**2)))
                        if (rel > tol) then
                            print *, 'complex full-rank residual:',m,n,family,rel
                            error stop 'complex full rank'
                        end if
                    end if
                    checks = checks+8
                    deallocate(u,v,s,uz,vz,rec,rz)
                end do

                do mode = 1, 3
                    select case (mode)
                    case (1)
                        tolerance = 0.15_dp
                    case (2)
                        tolerance = 1.0e-5_dp
                    case default
                        tolerance = 0.0_dp
                    end select
                    call randomized_interpolative_svd(a,tolerance,rank,ua,sa,va,status,seed=seed)
                    if (status /= linalg_status_success) then
                        print *, 'real adaptive error:',m,n,family,mode,status
                        error stop 'real adaptive SVD'
                    end if
                    if (rank < 0 .or. rank > min(m,n)) error stop 'real invalid rank'
                    if (size(ua,1) /= m .or. size(ua,2) /= rank) error stop 'adaptive real U shape'
                    if (size(va,1) /= n .or. size(va,2) /= rank) error stop 'adaptive real V shape'
                    if (size(sa) /= rank) error stop 'adaptive real S shape'
                    if (rank > 0) then
                        if (maxval(abs(matmul(transpose(ua),ua)-eye(rank))) > tol) &
                            error stop 'adaptive real U orthogonality'
                        if (maxval(abs(matmul(transpose(va),va)-eye(rank))) > tol) &
                            error stop 'adaptive real V orthogonality'
                    end if
                    rec = a - matmul(ua,matmul(diagonal(sa),transpose(va)))
                    rel = sqrt(sum((rec/scale)**2)) / max(1.0e-40_dp,sqrt(sum((a/scale)**2)))
                    if (rel > max(tol,tolerance*1.00001_dp)) then
                        print *, 'real adaptive residual:',m,n,family,mode,rank,rel
                        error stop 'real adaptive tolerance'
                    end if
                    call randomized_interpolative_svd(z,tolerance,rank,uaz,sa,vaz,status,seed=seed)
                    if (status /= linalg_status_success) then
                        print *, 'complex adaptive error:',m,n,family,mode,status
                        error stop 'complex adaptive SVD'
                    end if
                    if (rank > 0) then
                        if (maxval(abs(matmul(transpose(conjg(uaz)),uaz)-eye(rank))) > tol) &
                            error stop 'adaptive complex U orthogonality'
                        if (maxval(abs(matmul(transpose(conjg(vaz)),vaz)-eye(rank))) > tol) &
                            error stop 'adaptive complex V orthogonality'
                    end if
                    rz = z - matmul(uaz,matmul(diagonal(sa),transpose(conjg(vaz))))
                    rel = sqrt(sum(abs(rz/scale)**2)) / max(1.0e-40_dp,sqrt(sum(abs(z/scale)**2)))
                    if (rel > max(tol,tolerance*1.00001_dp)) then
                        print *, 'complex adaptive residual:',m,n,family,mode,rank,rel
                        error stop 'complex adaptive tolerance'
                    end if
                    checks = checks+8
                end do
                deallocate(ua,sa,va,uaz,vaz,rec,rz)
            end do
            deallocate(a,z)
        end do
    end do

    allocate(a(5,4),u(5,2),s(2),v(4,2))
    a = 1.0_dp
    call randomized_interpolative_svd(a,5,u,s,v,status,seed=seed)
    if (status /= linalg_status_invalid_shape) error stop 'invalid rank'
    call randomized_interpolative_svd(a,2,u,s,v,status,seed=seed,oversampling=-1)
    if (status /= linalg_status_invalid_shape) error stop 'negative oversampling'
    call randomized_interpolative_svd(a,2,u(:,1:1),s,v,status,seed=seed)
    if (status /= linalg_status_invalid_shape) error stop 'incorrect output shape'
    call randomized_interpolative_svd(a,1.0_dp,rank,ua,sa,va,status,seed=seed)
    if (status /= linalg_status_invalid_shape) error stop 'invalid tolerance'
    call randomized_interpolative_svd(a,-1.0_dp,rank,ua,sa,va,status,seed=seed)
    if (status /= linalg_status_invalid_shape) error stop 'negative tolerance'
    xnan = ieee_value(0.0_dp, ieee_quiet_nan)
    a(1,1) = xnan
    call randomized_interpolative_svd(a,2,u,s,v,status,seed=seed)
    if (status /= linalg_status_no_convergence) error stop 'nonfinite matrix'
    call randomized_interpolative_svd(a,0.01_dp,rank,ua,sa,va,status,seed=seed)
    if (status /= linalg_status_no_convergence) error stop 'nonfinite adaptive matrix'
    checks = checks+7

    a = 1.0_dp
    call randomized_interpolative_svd(a,2,u,s,v,status,seed=seed)
    previous = s
    call randomized_interpolative_svd(a,2,u,s,v,status,seed=seed)
    if (any(s /= previous)) error stop 'seeded repeatability'
    checks = checks+1

    print '(a,i0)', 'randomized interpolative SVD native checks passed: ',checks

contains

    function eye(n) result(x)
        integer, intent(in) :: n !! identity order
        real(dp) :: x(n,n)
        integer :: ii
        x = 0.0_dp
        do ii = 1,n
            x(ii,ii) = 1.0_dp
        end do
    end function eye

    function diagonal(values) result(x)
        real(dp), intent(in) :: values(:) !! entries
        real(dp) :: x(size(values),size(values))
        integer :: ii
        x = 0.0_dp
        do ii = 1,size(values)
            x(ii,ii) = values(ii)
        end do
    end function diagonal
end program test_linalg_interpolative_svd_random
