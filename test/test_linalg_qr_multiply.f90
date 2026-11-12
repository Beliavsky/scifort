! SPDX-License-Identifier: MIT
program test_linalg_qr_multiply
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_qr, only : qr, qr_factor
    use scifort_linalg_qr_multiply, only : qr_multiply, qr_apply
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none
    real(dp), allocatable :: a(:, :), q(:, :), r(:, :), rout(:, :)
    real(dp), allocatable :: c(:, :), result(:, :), expected(:, :), af(:, :), tau(:)
    complex(dp), allocatable :: ac(:, :), qc(:, :), rc(:, :), rcout(:, :)
    complex(dp), allocatable :: cc(:, :), resultc(:, :), expectedc(:, :), afc(:, :), tauc(:)
    complex(dp), allocatable :: cr(:, :), cm(:, :)
    integer, allocatable :: piv(:), pqr(:)
    integer :: m, n, k, nrhs, i, j, side, conjugated, pivoted, status, count
    real(dp), parameter :: tol = 6.0e-13_dp
    character(len=5) :: mode
    logical :: use_conj

    count = 0
    do m = 1, 9
        do n = 1, 9
            k = min(m,n)
            nrhs = mod(m+n,4)+1
            allocate(a(m,n),ac(m,n),q(m,k),r(k,n),rout(k,n))
            allocate(qc(m,k),rc(k,n),rcout(k,n),af(m,n),tau(k),afc(m,n),tauc(k))
            allocate(piv(n),pqr(n))
            do j = 1,n
                do i = 1,m
                    a(i,j) = sin(real(7*i+3*j,dp))/real(i+j,dp)
                    ac(i,j) = cmplx(a(i,j),cos(real(11*i+5*j,dp))/real(i+2*j,dp),dp)
                end do
            end do
            if (mod(m+n,5) == 0) then
                a(:,1) = 0.0_dp
                ac(:,1) = cmplx(0.0_dp,0.0_dp,dp)
            end if
            do pivoted = 0, 1
                if (pivoted == 1) then
                    call qr(a,q,r,status,pqr)
                    if (status /= linalg_status_success) error stop 'real pivot qr'
                    call qr(ac,qc,rc,status,pqr)
                    if (status /= linalg_status_success) error stop 'complex pivot qr'
                end if
                do side = 1, 2
                    mode = 'left'
                    if (side == 2) mode = 'right'
                    if (side == 1) then
                        allocate(c(k,nrhs),cc(k,nrhs),result(m,nrhs),resultc(m,nrhs))
                        allocate(expected(m,nrhs),expectedc(m,nrhs))
                    else
                        allocate(c(nrhs,m),cc(nrhs,m),result(nrhs,k),resultc(nrhs,k))
                        allocate(expected(nrhs,k),expectedc(nrhs,k))
                    end if
                    do j = 1, size(c,2)
                        do i = 1, size(c,1)
                            c(i,j) = cos(real(3*i+7*j,dp))
                            cc(i,j) = cmplx(c(i,j),sin(real(2*i+9*j,dp)),dp)
                        end do
                    end do
                    do conjugated = 0, 1
                        use_conj = (conjugated == 1)
                        do i = 0, 1
                            if (i == 0) then
                                if (pivoted == 1) then
                                    call qr(a,q,r,status,pqr)
                                    call qr_multiply(a,c,result,rout,status,mode,piv)
                                else
                                    call qr(a,q,r,status)
                                    call qr_multiply(a,c,result,rout,status,mode)
                                end if
                                if (status /= linalg_status_success) error stop 'real multiply status'
                                if (pivoted == 1) then
                                    if (any(piv /= pqr)) error stop 'real permutation'
                                end if
                                if (maxval(abs(r-rout)) > tol) error stop 'real R mismatch'
                                if (side == 1) then
                                    expected = matmul(q,c)
                                else
                                    expected = matmul(c,q)
                                end if
                                if (maxval(abs(result-expected)) > tol) error stop 'real product mismatch'
                                if (pivoted == 1) then
                                    call qr_factor(a,af,tau,status,pqr)
                                else
                                    call qr_factor(a,af,tau,status)
                                end if
                                call qr_apply(af,tau,c,result,status,mode,use_conj)
                                if (status /= linalg_status_success) error stop 'real reuse status'
                                if (maxval(abs(result-expected)) > tol) error stop 'real reuse mismatch'
                            else
                                if (pivoted == 1) then
                                    call qr(ac,qc,rc,status,pqr)
                                    call qr_multiply(ac,cc,resultc,rcout,status,mode,piv,use_conj)
                                else
                                    call qr(ac,qc,rc,status)
                                    call qr_multiply(ac,cc,resultc,rcout,status,mode,conjugate=use_conj)
                                end if
                                if (status /= linalg_status_success) error stop 'complex multiply status'
                                if (pivoted == 1) then
                                    if (any(piv /= pqr)) error stop 'complex permutation'
                                end if
                                if (maxval(abs(rc-rcout)) > tol) error stop 'complex R mismatch'
                                if (side == 1) then
                                    if (use_conj) then
                                        expectedc = matmul(conjg(qc),cc)
                                    else
                                        expectedc = matmul(qc,cc)
                                    end if
                                else
                                    if (use_conj) then
                                        expectedc = matmul(cc,conjg(qc))
                                    else
                                        expectedc = matmul(cc,qc)
                                    end if
                                end if
                                if (maxval(abs(resultc-expectedc)) > tol) error stop 'complex product mismatch'
                                if (pivoted == 1) then
                                    call qr_factor(ac,afc,tauc,status,pqr)
                                else
                                    call qr_factor(ac,afc,tauc,status)
                                end if
                                call qr_apply(afc,tauc,cc,resultc,status,mode,use_conj)
                                if (status /= linalg_status_success) error stop 'complex reuse status'
                                if (maxval(abs(resultc-expectedc)) > tol) error stop 'complex reuse mismatch'
                            end if
                            count = count + 1
                        end do
                    end do
                    ! Check mixed-type multiplication against the complex path.
                    allocate(cr(size(resultc,1),size(resultc,2)),cm(size(resultc,1),size(resultc,2)))
                    call qr_multiply(a,cc,cr,rcout,status,mode)
                    if (status /= linalg_status_success) error stop 'real A complex c'
                    call qr_multiply(cmplx(a,0.0_dp,dp),cc,cm,rcout,status,mode)
                    if (status /= linalg_status_success) error stop 'mixed comparison status'
                    if (maxval(abs(cr-cm)) > tol) error stop 'mixed result mismatch'
                    call qr_multiply(ac,c,cr,rcout,status,mode)
                    if (status /= linalg_status_success) error stop 'complex A real c'
                    call qr_multiply(ac,cmplx(c,0.0_dp,dp),cm,rcout,status,mode)
                    if (status /= linalg_status_success) error stop 'mixed comparison status 2'
                    if (maxval(abs(cr-cm)) > tol) error stop 'mixed result mismatch 2'
                    deallocate(cr,cm)
                    deallocate(c,cc,result,resultc,expected,expectedc)
                end do
            end do
            deallocate(a,ac,q,qc,r,rc,rout,rcout,af,afc,tau,tauc,piv,pqr)
        end do
    end do

    block
        real(dp) :: a(4,2), q(4,2), r(2,2), res(4), c(2), r2(2,2)
        real(dp) :: v(4), rightres(2), out(3,2), c1(2,2), empty_c(0,4), empty_y(0,2)
        real(dp) :: af(4,2), tau(2)
        complex(dp) :: ac(4,2), cv(2), outc(4), rc(2,2)
        real(dp) :: bad(4,2), goodout(4,2), badc(2,2), badtau(2), short_r(1,2)
        integer :: s, p(2)
        integer :: i
        do i = 1,4
            a(i,:) = [real(i,dp), real(i*i,dp)]
        end do
        call qr(a,q,r,s)
        c = [2.0_dp,3.0_dp]
        v = [4.0_dp,5.0_dp,6.0_dp,7.0_dp]
        call qr_multiply(a,c,res,r2,s,mode='left')
        if (s /= 0 .or. maxval(abs(res-matmul(q,c))) > tol) error stop 'real left vector'
        call qr_multiply(a,v,rightres,r2,s)
        if (s /= 0 .or. maxval(abs(rightres-matmul(v,q))) > tol) error stop 'real right vector'
        ac = cmplx(a,0.0_dp,dp)
        ac(:,2) = ac(:,2)*cmplx(0.0_dp,1.0_dp,dp)
        cv = [cmplx(2.0_dp,3.0_dp,dp),cmplx(4.0_dp,-2.0_dp,dp)]
        call qr_multiply(ac,cv,outc,rc,s,mode='left',conjugate=.true.)
        if (s /= 0) error stop 'complex left vector'
        call qr_multiply(a,cv,outc,rc,s,mode='left')
        if (s /= 0) error stop 'mixed vector real A'
        call qr_multiply(ac,c,outc,rc,s,mode='left')
        if (s /= 0) error stop 'mixed vector complex A'
        call qr_factor(a,af,tau,s)
        call qr_apply(af,tau,c,res,s,mode='left')
        if (s /= 0 .or. maxval(abs(res-matmul(q,c))) > tol) error stop 'reusable vector'
        c1 = 1.0_dp
        call qr_multiply(a,c1,out,r2,s,mode='left')
        if (s /= linalg_status_invalid_shape .or. any(out /= 0.0_dp)) &
            error stop 'invalid product shape'
        call qr_multiply(a,c1,goodout,r2,s,mode='invalid')
        if (s /= linalg_status_invalid_shape) error stop 'invalid mode'
        if (any(goodout /= 0.0_dp) .or. any(r2 /= 0.0_dp)) error stop 'error cleanup'
        call qr_multiply(a,c1,goodout,r2,s,mode='left',pivots=p(1:1))
        if (s /= linalg_status_invalid_shape) error stop 'invalid permutation length'
        call qr_multiply(a,c1,goodout,short_r,s,mode='left')
        if (s /= linalg_status_invalid_shape) error stop 'invalid R shape'
        bad = a
        bad(1,1) = ieee_value(0.0_dp,ieee_quiet_nan)
        call qr_multiply(bad,c1,goodout,r2,s,mode='left')
        if (s /= linalg_status_no_convergence) error stop 'NaN input must fail'
        if (any(goodout /= 0.0_dp) .or. any(r2 /= 0.0_dp)) error stop 'NaN cleanup'
        badc = c1
        badc(2,1) = ieee_value(0.0_dp,ieee_quiet_nan)
        call qr_multiply(a,badc,goodout,r2,s,mode='left')
        if (s /= linalg_status_no_convergence) error stop 'NaN multiplier must fail'
        badtau = tau
        badtau(1) = ieee_value(0.0_dp,ieee_quiet_nan)
        call qr_apply(af,badtau,c1,goodout,s,mode='left')
        if (s /= linalg_status_no_convergence) error stop 'NaN tau must fail'
        call qr_apply(af,tau(1:1),c1,goodout,s,mode='left')
        if (s /= linalg_status_invalid_shape) error stop 'wrong tau size'
        call qr_multiply(a,empty_c,empty_y,r2,s,mode='right')
        if (s /= linalg_status_success) error stop 'zero right-hand sides'
        call qr_multiply(a(1:0,:),c1,goodout,r2,s,mode='left')
        if (s /= linalg_status_invalid_shape) error stop 'zero-row A'
        call qr_multiply(a,c1,goodout,r2,s,mode='left',conjugate=.true.)
        if (s /= linalg_status_success) error stop 'real conjugation no-op' 
    end block
    print '(a,i0)', 'qr_multiply native shape/type/pivot/mode cases passed: ',count
end program test_linalg_qr_multiply
