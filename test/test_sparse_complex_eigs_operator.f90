! SPDX-License-Identifier: MIT
! Regression coverage for matrix-free complex Arnoldi and reusable shift factors.
program test_sparse_complex_eigs_operator
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csc_complex_matrix, &
        csr_complex_from_coo, csr_complex_to_csc
    use scifort_sparse_linalg_complex_eigs, only : eigs_complex, eigs_complex_operator, &
        eigs_complex_prepare_shift, complex_eigs_shift_cache, complex_eigs_success, &
        complex_eigs_invalid_input, complex_eigs_breakdown, complex_eigs_singular_shift
    implicit none
    type(csr_complex_matrix) :: csr, changed
    type(csc_complex_matrix) :: csc
    type(complex_eigs_shift_cache) :: cache, csc_cache, unready
    complex(dp), allocatable :: diagonal(:), upper(:), values(:), vectors(:,:), reference(:), &
        start(:), old_values(:)
    real(dp), allocatable :: residuals(:)
    integer, allocatable :: rows(:), cols(:)
    integer :: n, i, info, cycles, products, failures
    complex(dp) :: shift
    character(len=2) :: target

    n = 12
    allocate(diagonal(n),upper(n-1),rows(2*n-1),cols(2*n-1),reference(n))
    do i = 1,n
        diagonal(i) = cmplx(real(i,dp),0.09_dp*real(i-6,dp),dp)
        rows(i) = i
        cols(i) = i
    end do
    do i = 1,n-1
        upper(i) = cmplx(0.10_dp,-0.025_dp,dp)
        rows(n+i) = i
        cols(n+i) = i+1
    end do
    csr = csr_complex_from_coo(n,n,rows,cols,[diagonal,upper],info)
    if (info /= 0) error stop 'matrix-free test: CSR assembly failed'
    csc = csr_complex_to_csc(csr,info)
    if (info /= 0) error stop 'matrix-free test: CSC conversion failed'

    do i = 1,6
        select case (i)
        case (1)
            target = 'LM'
            reference(:2) = [diagonal(n),diagonal(n-1)]
        case (2)
            target = 'SM'
            reference(:2) = diagonal(:2)
        case (3)
            target = 'LR'
            reference(:2) = [diagonal(n),diagonal(n-1)]
        case (4)
            target = 'SR'
            reference(:2) = diagonal(:2)
        case (5)
            target = 'LI'
            reference(:2) = [diagonal(n),diagonal(n-1)]
        case (6)
            target = 'SI'
            reference(:2) = diagonal(:2)
        end select
        call eigs_complex_operator(n,2,apply_original,values,vectors,info, &
            which=target,ncv=n,tol=1.0e-10_dp,residual_norms=residuals, &
            iterations=cycles,matvecs=products,norm_estimate=14.0_dp)
        if (info /= complex_eigs_success) error stop 'matrix-free test: Arnoldi failed'
        if (maxval(abs(values-reference(:2))) > 1.0e-7_dp) &
            error stop 'matrix-free test: incorrect unshifted spectrum'
        if (cycles /= 1 .or. products < n) error stop 'matrix-free test: bad Arnoldi statistics'
        call verify_pairs(values,vectors,residuals,1.0e-7_dp)
    end do

    shift = diagonal(6)+cmplx(0.16_dp,0.06_dp,dp)
    call eigs_complex_operator(n,2,apply_original,values,vectors,info, &
        which='LM',sigma=shift,inverse_matvec=apply_inverse, &
        ncv=n,tol=1.0e-10_dp,residual_norms=residuals)
    if (info /= 0) error stop 'matrix-free test: analytic inverse failed'
    if (maxval(abs(values-[diagonal(6),diagonal(7)])) > 1.0e-7_dp) &
        error stop 'matrix-free test: wrong interior eigenvalues'
    call verify_pairs(values,vectors,residuals,1.0e-7_dp)

    call eigs_complex_prepare_shift(csr,shift,cache,info)
    if (info /= 0) error stop 'matrix-free test: cache preparation failed'
    call eigs_complex(csr,2,values,vectors,info,ncv=n,tol=1.0e-10_dp, &
        sigma=shift,shift_cache=cache,residual_norms=residuals)
    if (info /= 0) error stop 'matrix-free test: cached CSR solve failed'
    if (maxval(abs(values-[diagonal(6),diagonal(7)])) > 1.0e-7_dp) &
        error stop 'matrix-free test: cached CSR eigenvalues differ'
    call verify_pairs(values,vectors,residuals,1.0e-7_dp)
    old_values = values
    allocate(start(n))
    start = cmplx(1.0_dp,0.21_dp,dp)
    call eigs_complex(csr,2,values,vectors,info,sigma=shift, &
        shift_cache=cache,ncv=n,v0=start,tol=1.0e-10_dp)
    if (info /= 0 .or. maxval(abs(values-old_values)) > 1.0e-7_dp) &
        error stop 'matrix-free test: repeat cached solve failed'
    call eigs_complex(csc,2,values,vectors,info,sigma=shift, &
        shift_cache=cache,ncv=n,tol=1.0e-10_dp)
    if (info /= 0 .or. maxval(abs(values-old_values)) > 1.0e-7_dp) &
        error stop 'matrix-free test: cached CSC solve failed'
    call eigs_complex_prepare_shift(csc,shift,csc_cache,info)
    if (info /= 0) error stop 'matrix-free test: CSC cache preparation failed'
    call eigs_complex(csc,2,values,vectors,info,sigma=shift, &
        shift_cache=csc_cache,ncv=n,tol=1.0e-10_dp)
    if (info /= 0) error stop 'matrix-free test: CSC-origin cache rejected'

    changed = csr
    changed%data(1) = changed%data(1) + cmplx(0.4_dp,0.1_dp,dp)
    call eigs_complex(changed,2,values,vectors,info,sigma=shift,shift_cache=cache,ncv=n)
    if (info /= complex_eigs_invalid_input .or. size(values) /= 0) &
        error stop 'matrix-free test: stale matrix cache accepted'
    call eigs_complex(csr,2,values,vectors,info,sigma=shift+0.1_dp,shift_cache=cache,ncv=n)
    if (info /= complex_eigs_invalid_input) error stop 'matrix-free test: stale shift accepted'
    call eigs_complex(csr,2,values,vectors,info,shift_cache=cache,ncv=n)
    if (info /= complex_eigs_invalid_input) error stop 'matrix-free test: missing shift accepted'
    call eigs_complex(csr,2,values,vectors,info,sigma=shift,shift_cache=unready,ncv=n)
    if (info /= complex_eigs_invalid_input) error stop 'matrix-free test: invalid cache accepted'
    call eigs_complex_prepare_shift(csr,diagonal(5),unready,info)
    if (info /= complex_eigs_singular_shift) error stop 'matrix-free test: singular cache accepted'

    ! Signature validation must not invoke any callback.
    call eigs_complex_operator(n,2,apply_original,values,vectors,info,sigma=shift,ncv=n)
    if (info /= complex_eigs_invalid_input) error stop 'matrix-free test: missing inverse accepted'
    call eigs_complex_operator(n,2,apply_original,values,vectors,info,inverse_matvec=apply_inverse,ncv=n)
    if (info /= complex_eigs_invalid_input) error stop 'matrix-free test: inverse without shift accepted'
    call eigs_complex_operator(n,2,apply_original,values,vectors,info,ncv=2)
    if (info /= complex_eigs_invalid_input) error stop 'matrix-free test: invalid ncv accepted'
    call eigs_complex_operator(n,2,apply_original,values,vectors,info,ncv=n,norm_estimate=-1.0_dp)
    if (info /= complex_eigs_invalid_input) error stop 'matrix-free test: negative norm accepted'
    call eigs_complex_operator(n,2,apply_original,values,vectors,info,ncv=n, &
        v0=cmplx(0.0_dp,0.0_dp,dp)*start)
    if (info /= complex_eigs_invalid_input) error stop 'matrix-free test: zero start accepted'
    call eigs_complex_operator(n,2,apply_original,values,vectors,info,ncv=n, &
        sigma=cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp),inverse_matvec=apply_inverse)
    if (info /= complex_eigs_invalid_input) error stop 'matrix-free test: nonfinite shift accepted'
    do failures = 1,2
        call eigs_complex_operator(n,2,apply_failure,values,vectors,info,ncv=n)
        if (info /= complex_eigs_breakdown) error stop 'matrix-free test: invalid callback success'
    end do
    failures = 3
    call eigs_complex_operator(n,2,apply_original,values,vectors,info,ncv=n, &
        sigma=shift,inverse_matvec=apply_inverse_failure)
    if (info /= complex_eigs_breakdown) error stop 'matrix-free test: inverse callback error accepted'
    failures = 0

    ! A large diagonal + superdiagonal operator is represented only by its
    ! matrix-free coefficients, not by a dense or CSR matrix.
    n = 2500
    deallocate(diagonal,upper)
    allocate(diagonal(n),upper(n-1))
    do i = 1,n
        diagonal(i) = cmplx(real(i,dp),0.004_dp*real(i,dp),dp)
    end do
    upper = cmplx(0.003_dp,0.001_dp,dp)
    shift = diagonal(1200)+cmplx(0.28_dp,-0.1_dp,dp)
    call eigs_complex_operator(n,2,apply_original,values,vectors,info,ncv=20, &
        maxiter=32,tol=5.0e-9_dp,sigma=shift,inverse_matvec=apply_inverse, &
        norm_estimate=real(n,dp),residual_norms=residuals,matvecs=products)
    if (info /= complex_eigs_success) error stop 'matrix-free test: large matrix-free Arnoldi failed'
    if (maxval(abs(values-[diagonal(1200),diagonal(1201)])) > 1.0e-3_dp) &
        error stop 'matrix-free test: large interior eigenvalues wrong'
    if (products < 20) error stop 'matrix-free test: expected inverse applications'
    call verify_pairs(values,vectors,residuals,5.0e-5_dp)

    print '(a)', 'test_sparse_complex_eigs_operator: PASS'
contains
    subroutine apply_original(x,y,status)
        complex(dp), intent(in) :: x(:) !! Vector to multiply.
        complex(dp), intent(out) :: y(:) !! A*x.
        integer, intent(out) :: status !! Zero for success.
        y = diagonal*x
        y(:n-1) = y(:n-1)+upper*x(2:n)
        status = 0
    end subroutine apply_original

    subroutine apply_inverse(x,y,status)
        complex(dp), intent(in) :: x(:) !! Right-hand side of shifted system.
        complex(dp), intent(out) :: y(:) !! Solution of (A-sigma I)*y=x.
        integer, intent(out) :: status !! Zero for success.
        integer :: j
        y(n) = x(n)/(diagonal(n)-shift)
        do j = n-1,1,-1
            y(j) = (x(j)-upper(j)*y(j+1))/(diagonal(j)-shift)
        end do
        status = 0
    end subroutine apply_inverse

    subroutine apply_failure(x,y,status)
        complex(dp), intent(in) :: x(:) !! Input vector.
        complex(dp), intent(out) :: y(:) !! Output may be undefined on failure.
        integer, intent(out) :: status !! Injected error status.
        if (failures == 1) then
            status = 17
        else
            status = 0
            y = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
        end if
    end subroutine apply_failure

    subroutine apply_inverse_failure(x,y,status)
        complex(dp), intent(in) :: x(:) !! Input vector.
        complex(dp), intent(out) :: y(:) !! Output undefined on failure.
        integer, intent(out) :: status !! Injected error status.
        status = 11
    end subroutine apply_inverse_failure

    subroutine verify_pairs(lam,vec,error,threshold)
        complex(dp), intent(in) :: lam(:) !! Returned eigenvalues.
        complex(dp), intent(in) :: vec(:,:) !! Returned eigenvectors.
        real(dp), intent(in) :: error(:) !! Reported true residual norms.
        real(dp), intent(in) :: threshold !! Absolute accuracy threshold.
        complex(dp), allocatable :: y(:)
        integer :: j,status
        if (size(lam) /= size(error) .or. size(vec,1) /= n .or. &
            size(vec,2) /= size(lam)) error stop 'matrix-free test: wrong result shape'
        allocate(y(n))
        do j = 1,size(lam)
            call apply_original(vec(:,j),y,status)
            if (status /= 0) error stop 'matrix-free test: multiplication failed'
            if (norm2(abs(y-lam(j)*vec(:,j))) > threshold) &
                error stop 'matrix-free test: invalid eigenvector residual'
            if (abs(norm2(abs(vec(:,j)))-1.0_dp) > 1.0e-8_dp) &
                error stop 'matrix-free test: unnormalized eigenvector'
            if (error(j) > threshold) error stop 'matrix-free test: false convergence'
        end do
    end subroutine verify_pairs
end program test_sparse_complex_eigs_operator
