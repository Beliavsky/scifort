! SPDX-License-Identifier: MIT
program test_linalg_svd_full
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_svd, only : svd, svdvals, diagsvd
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none
    integer :: m, n, k, stat, i, j, kind_case
    real(dp), allocatable :: d(:, :)
    real(dp) :: err, tol
    real(dp), allocatable :: ar(:, :), ur(:, :), vhr(:, :), sr(:), ss(:)
    complex(dp), allocatable :: ac(:, :), uc(:, :), vhc(:, :)

    tol = 2.0e-12_dp
    do kind_case = 1, 2
        do m = 1, 8
            do n = 1, 8
                k = min(m,n)
                allocate(ar(m,n),ur(m,m),vhr(n,n),sr(k),ss(k),d(m,n))
                allocate(ac(m,n),uc(m,m),vhc(n,n))
                do j = 1, n
                    do i = 1, m
                        ar(i,j) = real(mod(3*i+7*j,11)-5,dp)/10.0_dp
                        ac(i,j) = cmplx(ar(i,j),real(mod(5*i+j,9)-4,dp)/13.0_dp,dp)
                    end do
                end do
                if (kind_case == 2) then
                    ar = 0.0_dp
                    ac = cmplx(0.0_dp,0.0_dp,dp)
                    if (m >= 1 .and. n >= 1) then
                        ar(1,1) = 1.0_dp
                        ac(1,1) = cmplx(2.0_dp,-1.0_dp,dp)
                    end if
                end if
                call svd(ar,ur,sr,vhr,stat,full_matrices=.true.)
                if (stat /= linalg_status_success) error stop 'real full SVD status'
                err = maxval(abs(matmul(ur(:,1:k)*spread(sr,1,m),vhr(1:k,:))-ar))
                if (err > tol) error stop 'real full reconstruction'
                if (maxval(abs(matmul(transpose(ur),ur)-identity_real(m))) > tol) &
                    error stop 'real full U unitary'
                if (maxval(abs(matmul(vhr,transpose(vhr))-identity_real(n))) > tol) &
                    error stop 'real full V unitary'
                call diagsvd(sr,m,n,d,stat)
                if (stat /= linalg_status_success) error stop 'diagsvd real status'
                if (maxval(abs(matmul(matmul(ur,d),vhr)-ar)) > tol) &
                    error stop 'full SVD via diagsvd'
                call svdvals(ar,ss,stat)
                if (stat /= linalg_status_success) error stop 'real svdvals status'
                if (maxval(abs(sr-ss)) > tol) error stop 'real svdvals mismatch'
                call svd(ac,uc,sr,vhc,stat,full_matrices=.true.)
                if (stat /= linalg_status_success) error stop 'complex full SVD status'
                err = maxval(abs(matmul(uc(:,1:k)*spread(sr,1,m),vhc(1:k,:))-ac))
                if (err > tol) error stop 'complex full reconstruction'
                if (maxval(abs(matmul(transpose(conjg(uc)),uc)-identity_complex(m))) > tol) &
                    error stop 'complex full U unitary'
                if (maxval(abs(matmul(vhc,transpose(conjg(vhc)))-identity_complex(n))) > tol) &
                    error stop 'complex full V unitary'
                call svdvals(ac,ss,stat)
                if (stat /= linalg_status_success) error stop 'complex svdvals status'
                if (maxval(abs(sr-ss)) > tol) error stop 'complex svdvals mismatch'
                deallocate(ar,ur,vhr,sr,ss,ac,uc,vhc,d)
            end do
        end do
    end do

    ! Explicitly verify shape and nonfinite diagnostics, with zeroed output.
    allocate(ar(5,3),ur(5,5),vhr(3,3),sr(3),ss(3))
    ar = 0.0_dp
    call svd(ar,ur(:,1:3),sr,vhr,stat,full_matrices=.true.)
    if (stat /= linalg_status_invalid_shape) error stop 'full U invalid shape'
    call svd(ar,ur,sr,vhr(1:2,:),stat,full_matrices=.true.)
    if (stat /= linalg_status_invalid_shape) error stop 'full VH invalid shape'
    call svd(ar,ur,sr,vhr,stat)
    if (stat /= linalg_status_invalid_shape) error stop 'economy U invalid shape'
    call svd(ar,ur,sr,vhr,stat,max_sweeps=0,full_matrices=.true.)
    if (stat /= linalg_status_invalid_shape) error stop 'full invalid sweeps'
    call svdvals(ar,ss(:2),stat)
    if (stat /= linalg_status_invalid_shape) error stop 'svdvals invalid shape'
    call svdvals(ar,ss,stat,max_sweeps=0)
    if (stat /= linalg_status_invalid_shape) error stop 'svdvals invalid sweeps'
    allocate(d(5,3))
    call diagsvd(sr(:2),5,3,d,stat)
    if (stat /= linalg_status_invalid_shape) error stop 'diagsvd invalid count'
    if (any(d /= 0.0_dp)) error stop 'diagsvd error output'
    call diagsvd(sr,5,3,d,stat)
    if (stat /= linalg_status_success) error stop 'diagsvd valid status'
    if (any(d /= 0.0_dp)) error stop 'diagsvd zero inputs'
    sr(1) = ieee_value(1.0_dp,ieee_quiet_nan)
    call diagsvd(sr,5,3,d,stat)
    if (stat /= linalg_status_no_convergence) error stop 'diagsvd nonfinite'
    if (any(d /= 0.0_dp)) error stop 'diagsvd nonfinite output'
    sr = 0.0_dp
    deallocate(d)
    allocate(d(0,3))
    call diagsvd(sr(:0),0,3,d,stat)
    if (stat /= linalg_status_success) error stop 'diagsvd empty shape'
    ar(1,1) = ieee_value(1.0_dp,ieee_quiet_nan)
    call svd(ar,ur,sr,vhr,stat,full_matrices=.true.)
    if (stat /= linalg_status_no_convergence) error stop 'full nonfinite status'
    if (any(sr /= 0.0_dp) .or. any(ur /= 0.0_dp) .or. &
        any(vhr /= 0.0_dp)) error stop 'full nonfinite outputs'
    call svdvals(ar,ss,stat)
    if (stat /= linalg_status_no_convergence) error stop 'svdvals nonfinite status'
    if (any(ss /= 0.0_dp)) error stop 'svdvals nonfinite outputs'

    print '(A)', 'full real/complex SVD, svdvals, diagsvd and diagnostics passed'
contains
    function identity_real(n) result(id)
        integer, intent(in) :: n !! square dimension
        real(dp) :: id(n,n)
        integer :: ii
        id = 0.0_dp
        do ii = 1, n
            id(ii,ii) = 1.0_dp
        end do
    end function identity_real

    function identity_complex(n) result(id)
        integer, intent(in) :: n !! square dimension
        complex(dp) :: id(n,n)
        integer :: ii
        id = cmplx(0.0_dp,0.0_dp,dp)
        do ii = 1, n
            id(ii,ii) = cmplx(1.0_dp,0.0_dp,dp)
        end do
    end function identity_complex
end program test_linalg_svd_full
