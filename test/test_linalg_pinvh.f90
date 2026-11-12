! SPDX-License-Identifier: MIT
program test_linalg_pinvh
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_pinvh, only : pinvh
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none

    real(dp), allocatable :: a(:,:), p(:,:), b(:,:), batch(:,:,:), batch_inv(:,:,:)
    complex(dp), allocatable :: z(:,:), zi(:,:), bz(:,:,:), bz_inv(:,:,:)
    integer, allocatable :: ranks(:)
    real(dp) :: residual, norm_a, nan_value
    integer :: n, i, j, variant, direction, status, rank, k, checks
    real(dp), parameter :: tol = 3.0e-9_dp

    checks = 0
    do n = 0, 12
        allocate(a(n,n),p(n,n),b(n,n),z(n,n),zi(n,n))
        do variant = 1, 6
            b = 0.0_dp
            z = cmplx(0.0_dp,0.0_dp,dp)
            do j = 1, n
                do i = 1, n
                    b(i,j) = sin(real(2*i+7*j+variant,dp))/real(2+i+j,dp)
                    z(i,j) = cmplx(b(i,j),cos(real(5*i+3*j+variant,dp))/real(2+i+j,dp),dp)
                end do
            end do
            a = 0.5_dp*(b+transpose(b))
            z = 0.5_dp*(z+transpose(conjg(z)))
            select case (variant)
            case (2)
                a = 0.0_dp
                z = cmplx(0.0_dp,0.0_dp,dp)
            case (3,4)
                a = 0.0_dp
                z = cmplx(0.0_dp,0.0_dp,dp)
                do i = 1, n
                    if (i == n .and. variant == 3) cycle
                    a(i,i) = (-1.0_dp)**i*real(i,dp)
                    z(i,i) = cmplx(a(i,i),0.0_dp,dp)
                end do
            case (5)
                if (n > 2) then
                    b(:,3:n) = 0.0_dp
                    z(:,3:n) = cmplx(0.0_dp,0.0_dp,dp)
                end if
                a = matmul(b,transpose(b))
                z = matmul(z,transpose(conjg(z)))
            case (6)
                a = a * 1.0e-12_dp
                z = z * 1.0e-12_dp
            end select
            do direction = 1,2
                if (direction == 1) then
                    call pinvh(a,p,status,rank,lower=.true., &
                        rtol=merge(1.0e-9_dp,real(n,dp)*epsilon(1.0_dp),variant==5))
                    if (status /= linalg_status_success) error stop 'real pinvh status'
                    call pinvh(z,zi,status,lower=.true., &
                        rtol=merge(1.0e-9_dp,real(n,dp)*epsilon(1.0_dp),variant==5))
                    if (status /= linalg_status_success) error stop 'complex pinvh status'
                else
                    call pinvh(a,p,status,rank,lower=.false., &
                        rtol=merge(1.0e-9_dp,real(n,dp)*epsilon(1.0_dp),variant==5))
                    if (status /= linalg_status_success) error stop 'real pinvh upper status'
                    call pinvh(z,zi,status,lower=.false., &
                        rtol=merge(1.0e-9_dp,real(n,dp)*epsilon(1.0_dp),variant==5))
                    if (status /= linalg_status_success) error stop 'complex pinvh upper status'
                end if
                if (n > 0) then
                    norm_a = max(1.0e-90_dp,maxval(abs(a)))
                    residual = maxval(abs(matmul(a,matmul(p,a))-a))/norm_a
                    if (residual > tol) then
                        print *, 'bad real penrose',n,variant,direction,residual,rank
                        error stop 'real penrose A P A'
                    end if
                    norm_a = max(1.0e-90_dp,maxval(abs(z)))
                    residual = maxval(abs(matmul(z,matmul(zi,z))-z))/norm_a
                    if (residual > tol) then
                        print *, 'bad complex penrose',n,variant,direction,residual
                        error stop 'complex penrose A P A'
                    end if
                    if (maxval(abs(p-transpose(p))) > tol * max(1.0e-90_dp,maxval(abs(p)))) &
                        error stop 'real inverse symmetry'
                    if (maxval(abs(zi-transpose(conjg(zi)))) > &
                        tol * max(1.0e-90_dp,maxval(abs(zi)))) error stop 'complex inverse hermiticity'
                end if
                checks = checks+2
            end do
        end do
        deallocate(a,p,b,z,zi)
    end do

    ! The criterion uses eigenvalue magnitudes, not positive values only.
    allocate(a(3,3),p(3,3),z(3,3),zi(3,3))
    a = 0.0_dp
    a(1,1) = -4.0_dp
    a(2,2) = 0.25_dp
    a(3,3) = 0.001_dp
    call pinvh(a,p,status,rank,rtol=0.1_dp)
    if (status /= linalg_status_success .or. rank /= 1) error stop 'relative rank cutoff'
    if (abs(p(1,1)+0.25_dp) > 1.0e-12_dp) error stop 'negative eigenvalue reciprocal'
    if (any(p(2:3,2:3) /= 0.0_dp)) error stop 'truncated small eigenvalues'
    call pinvh(a,p,status,rank,atol=0.01_dp,rtol=0.0_dp)
    if (status /= linalg_status_success .or. rank /= 2) error stop 'absolute rank cutoff'
    if (abs(p(2,2)-4.0_dp) > 1.0e-12_dp) error stop 'small retained eigenvalue'
    call pinvh(a,p,status,rank,atol=0.0_dp,rtol=0.0_dp)
    if (status /= linalg_status_success .or. rank /= 3) error stop 'zero tolerance rank'

    ! An ignored triangle is not used in the eigensolver; finite checks are
    ! separate and (as in SciPy) examine the entire input by default.
    a = 0.0_dp
    a(1,1) = 2.0_dp
    a(2,1) = 1.0_dp
    a(2,2) = 3.0_dp
    a(3,3) = 4.0_dp
    a(1,2) = 1000.0_dp
    call pinvh(a,p,status,lower=.true., &
                        rtol=merge(1.0e-9_dp,real(n,dp)*epsilon(1.0_dp),variant==5))
    if (status /= linalg_status_success) error stop 'lower selected triangle'
    if (maxval(abs(matmul(reshape([2.0_dp,1.0_dp,0.0_dp, &
        1.0_dp,3.0_dp,0.0_dp,0.0_dp,0.0_dp,4.0_dp],[3,3]),p)- &
        identity3())) > 1.0e-11_dp) error stop 'lower triangle ignored upper'
    a = transpose(a)
    call pinvh(a,p,status,lower=.false., &
                        rtol=merge(1.0e-9_dp,real(n,dp)*epsilon(1.0_dp),variant==5))
    if (status /= linalg_status_success) error stop 'upper selected triangle'
    if (maxval(abs(matmul(reshape([2.0_dp,1.0_dp,0.0_dp, &
        1.0_dp,3.0_dp,0.0_dp,0.0_dp,0.0_dp,4.0_dp],[3,3]),p)- &
        identity3())) > 1.0e-11_dp) error stop 'upper triangle ignored lower'

    nan_value = ieee_value(1.0_dp,ieee_quiet_nan)
    a(2,1) = nan_value
    call pinvh(a,p,status)
    if (status /= linalg_status_no_convergence) error stop 'default finite checks'
    call pinvh(a,p,status,check_finite=.false.,lower=.false.)
    if (status /= linalg_status_success) error stop 'ignored nonfinite triangle'
    if (any(.not. (p == p))) error stop 'output finite with ignored nonfinite'
    call pinvh(a,p,status,atol=-1.0_dp)
    if (status /= linalg_status_invalid_shape) error stop 'invalid atol'
    if (any(p /= 0.0_dp)) error stop 'invalid atol output not zero'
    call pinvh(a,p,status,rtol=-1.0_dp)
    if (status /= linalg_status_invalid_shape) error stop 'invalid rtol'
    call pinvh(a,p,status,max_iterations=0)
    if (status /= linalg_status_invalid_shape) error stop 'invalid iterations'
    call pinvh(a,p(1:2,:),status)
    if (status /= linalg_status_invalid_shape) error stop 'bad output shape'
    checks = checks+12

    ! Complex offdiagonal phase, and imaginary diagonal in ignored half.
    z = cmplx(0.0_dp,0.0_dp,dp)
    z(1,1) = cmplx(4.0_dp,0.0_dp,dp)
    z(2,2) = cmplx(3.0_dp,0.0_dp,dp)
    z(3,3) = cmplx(-2.0_dp,0.0_dp,dp)
    z(2,1) = cmplx(1.0_dp,-2.0_dp,dp)
    z(1,2) = conjg(z(2,1))
    call pinvh(z,zi,status,rank)
    if (status /= linalg_status_success .or. rank /= 3) error stop 'complex phase'
    if (maxval(abs(matmul(z,zi)-cmplx(identity3(),0.0_dp,dp))) > 1.0e-10_dp) &
        error stop 'complex inverse residual'
    z(1,2) = cmplx(nan_value,0.0_dp,dp)
    call pinvh(z,zi,status,check_finite=.false.)
    if (status /= linalg_status_success) error stop 'complex ignored nonfinite'
    checks = checks+3

    deallocate(a,p,z,zi)
    allocate(batch(3,2,2),batch_inv(3,2,2),ranks(3),bz(3,2,2),bz_inv(3,2,2))
    batch = 0.0_dp
    bz = cmplx(0.0_dp,0.0_dp,dp)
    do k = 1,3
        batch(k,1,1) = real(k,dp)
        batch(k,2,2) = 0.0_dp
        bz(k,1,1) = cmplx(real(k,dp),0.0_dp,dp)
    end do
    call pinvh(batch,batch_inv,status,ranks)
    if (status /= linalg_status_success .or. any(ranks /= 1)) error stop 'real batch ranks'
    call pinvh(bz,bz_inv,status,ranks)
    if (status /= linalg_status_success .or. any(ranks /= 1)) error stop 'complex batch ranks'
    if (abs(batch_inv(2,1,1)-0.5_dp) > 1.0e-13_dp) error stop 'batch inverse values'
    if (abs(bz_inv(2,1,1)-cmplx(0.5_dp,0.0_dp,dp)) > 1.0e-13_dp) &
        error stop 'complex batch inverse values'
    batch(3,2,2) = nan_value
    call pinvh(batch,batch_inv,status,ranks)
    if (status /= linalg_status_no_convergence) error stop 'batch member failure'
    if (any(batch_inv /= 0.0_dp) .or. any(ranks /= 0)) error stop 'batch cleanup'
    call pinvh(batch,batch_inv(1:2,:,:),status)
    if (status /= linalg_status_invalid_shape) error stop 'batch shape'
    call pinvh(batch,batch_inv,status,ranks(1:2))
    if (status /= linalg_status_invalid_shape) error stop 'batch rank shape'
    call pinvh(bz(1:0,:,:),bz_inv(1:0,:,:),status)
    if (status /= linalg_status_invalid_shape) error stop 'empty batch'
    checks = checks+8

    print '(a,i0)', 'PASS pinvh native cases: ',checks

contains
    pure function identity3() result(y)
        real(dp) :: y(3,3)
        integer :: t
        y = 0.0_dp
        do t = 1,3
            y(t,t) = 1.0_dp
        end do
    end function identity3
end program test_linalg_pinvh
