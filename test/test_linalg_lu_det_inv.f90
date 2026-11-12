! SPDX-License-Identifier: MIT
program test_linalg_lu_det_inv
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_singular, linalg_status_no_convergence
    use scifort_linalg_lu, only : lu, lu_indices, lu_permuted
    use scifort_linalg_det_inv, only : det, inv
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none
    real(dp), allocatable :: a(:, :), p(:, :), l(:, :), u(:, :), pl(:, :), inverse(:, :), id(:, :)
    complex(dp), allocatable :: z(:, :), pz(:, :), lz(:, :), uz(:, :), plz(:, :), iz(:, :), idz(:, :)
    integer, allocatable :: perm(:)
    real(dp) :: d, scale, nan_value
    complex(dp) :: dz
    integer :: m, n, k, i, j, status, variant, count

    count = 0
    do m = 0, 13
        do n = 0, 13
            k = min(m,n)
            allocate(a(m,n),z(m,n),p(m,m),pz(m,m),l(m,k),u(k,n),lz(m,k),uz(k,n))
            allocate(pl(m,k),plz(m,k),perm(m))
            do variant = 1, 6
                do j = 1, n
                    do i = 1, m
                        a(i,j) = sin(real(7*i+11*j+variant,dp)) / real(i+j+1,dp)
                        z(i,j) = cmplx(a(i,j),cos(real(3*i+5*j+variant,dp)),dp)
                    end do
                end do
                if (variant == 2 .and. m > 1) then
                    a(m,:) = a(1,:)
                    z(m,:) = z(1,:)
                else if (variant == 3 .and. n > 0) then
                    a(:,1) = 0.0_dp
                    z(:,1) = cmplx(0.0_dp,0.0_dp,dp)
                else if (variant == 4) then
                    a = 0.0_dp
                    z = cmplx(0.0_dp,0.0_dp,dp)
                else if (variant == 5) then
                    a = a * 1.0e-40_dp
                    z = z * 1.0e-40_dp
                else if (variant == 6 .and. m*n > 0) then
                    a(m,n) = 1.0e30_dp
                    z(m,n) = cmplx(1.0e30_dp,1.0e30_dp,dp)
                end if
                call lu(a,p,l,u,status)
                if (status /= linalg_status_success) error stop 'lu real status'
                call lu(z,pz,lz,uz,status)
                if (status /= linalg_status_success) error stop 'lu complex status'
                call lu_indices(a,perm,pl,u,status)
                if (status /= linalg_status_success) error stop 'lu real indices status'
                if (m > 0) then
                    do i = 1, m
                        if (perm(i) < 1 .or. perm(i) > m) error stop 'lu real index range'
                        if (p(i,perm(i)) /= 1.0_dp) error stop 'lu P/indices mismatch'
                    end do
                end if
                if (m*n > 0) then
                    scale = max(1.0e-90_dp,maxval(abs(a)))
                    if (maxval(abs(a-matmul(p,matmul(pl,u)))) > 4.0e-12_dp*scale) &
                        error stop 'lu real reconstruction'
                    scale = max(1.0e-90_dp,maxval(abs(z)))
                    if (maxval(abs(z-matmul(pz,matmul(lz,uz)))) > 4.0e-12_dp*scale) &
                        error stop 'lu complex reconstruction'
                    do j = 1, k
                        if (pl(j,j) /= 1.0_dp .or. lz(j,j) /= cmplx(1.0_dp,0.0_dp,dp)) &
                            error stop 'lu unit diagonal'
                        if (j > 1) then
                            if (any(pl(1:j-1,j) /= 0.0_dp)) error stop 'lu lower structure'
                            if (any(lz(1:j-1,j) /= cmplx(0.0_dp,0.0_dp,dp))) &
                                error stop 'lu complex lower structure'
                            if (any(u(j,1:j-1) /= 0.0_dp)) error stop 'lu upper structure'
                        end if
                    end do
                end if
                call lu_permuted(a,pl,u,status)
                if (status /= linalg_status_success) error stop 'lu_permuted real'
                call lu_permuted(z,plz,uz,status)
                if (status /= linalg_status_success) error stop 'lu_permuted complex'
                if (m*n > 0) then
                    scale = max(1.0e-90_dp,maxval(abs(a)))
                    if (maxval(abs(a-matmul(pl,u))) > 4.0e-12_dp*scale) &
                        error stop 'lu PL reconstruction real'
                    scale = max(1.0e-90_dp,maxval(abs(z)))
                    if (maxval(abs(z-matmul(plz,uz))) > 4.0e-12_dp*scale) &
                        error stop 'lu PL reconstruction complex'
                end if
                count = count+2
            end do
            deallocate(a,z,p,pz,l,u,lz,uz,pl,plz,perm)
        end do
    end do

    do n = 0, 15
        allocate(a(n,n),z(n,n),inverse(n,n),iz(n,n),id(n,n),idz(n,n))
        id = 0.0_dp
        idz = cmplx(0.0_dp,0.0_dp,dp)
        do i = 1, n
            id(i,i) = 1.0_dp
            idz(i,i) = cmplx(1.0_dp,0.0_dp,dp)
        end do
        a = 0.0_dp
        z = cmplx(0.0_dp,0.0_dp,dp)
        do i = 1, n
            a(i,i) = real(i,dp)/4.0_dp
            z(i,i) = cmplx(real(i,dp)/4.0_dp,real(i,dp)/7.0_dp,dp)
        end do
        call det(a,d,status)
        if (status /= linalg_status_success) error stop 'det real diagonal'
        call det(z,dz,status)
        if (status /= linalg_status_success) error stop 'det complex diagonal'
        call inv(a,inverse,status)
        if (status /= linalg_status_success) error stop 'inv real diagonal'
        call inv(z,iz,status)
        if (status /= linalg_status_success) error stop 'inv complex diagonal'
        if (n > 0) then
            if (maxval(abs(matmul(a,inverse)-id)) > 1.0e-12_dp) error stop 'inv real identity'
            if (maxval(abs(matmul(z,iz)-idz)) > 1.0e-12_dp) error stop 'inv complex identity'
        else
            if (d /= 1.0_dp .or. dz /= cmplx(1.0_dp,0.0_dp,dp)) error stop 'det empty'
        end if
        if (n > 0) then
            a(1,:) = 0.0_dp
            z(1,:) = cmplx(0.0_dp,0.0_dp,dp)
            call det(a,d,status)
            if (status /= linalg_status_success .or. d /= 0.0_dp) error stop 'det singular real'
            call det(z,dz,status)
            if (status /= linalg_status_success .or. dz /= cmplx(0.0_dp,0.0_dp,dp)) &
                error stop 'det singular complex'
            call inv(a,inverse,status)
            if (status /= linalg_status_singular .or. any(inverse /= 0.0_dp)) &
                error stop 'inv singular real'
            call inv(z,iz,status)
            if (status /= linalg_status_singular .or. any(iz /= cmplx(0.0_dp,0.0_dp,dp))) &
                error stop 'inv singular complex'
        end if
        count = count+4
        deallocate(a,z,inverse,iz,id,idz)
    end do

    allocate(a(2,3),z(2,3),inverse(2,3),iz(2,3))
    call det(a,d,status)
    if (status /= linalg_status_invalid_shape) error stop 'det nonsquare'
    call inv(a,inverse,status)
    if (status /= linalg_status_invalid_shape) error stop 'inv nonsquare'
    deallocate(a,z,inverse,iz)
    allocate(a(2,2),p(2,2),l(2,2),u(2,2),perm(2),z(2,2),pz(2,2),lz(2,2),uz(2,2))
    a = 1.0_dp
    z = cmplx(1.0_dp,0.0_dp,dp)
    nan_value = ieee_value(1.0_dp,ieee_quiet_nan)
    a(1,1) = nan_value
    z(1,1) = cmplx(nan_value,0.0_dp,dp)
    call lu(a,p,l,u,status)
    if (status /= linalg_status_no_convergence) error stop 'lu nonfinite real'
    call lu(z,pz,lz,uz,status)
    if (status /= linalg_status_no_convergence) error stop 'lu nonfinite complex'
    call det(a,d,status)
    if (status /= linalg_status_no_convergence) error stop 'det nonfinite'
    call det(z,dz,status)
    if (status /= linalg_status_no_convergence) error stop 'det complex nonfinite'
    call lu_indices(a,perm,l(:,1:1),u,status)
    if (status /= linalg_status_invalid_shape) error stop 'lu bad shape'

    print '(a,i0)', 'PASS LU/det/inv native cases: ',count
end program test_linalg_lu_det_inv
