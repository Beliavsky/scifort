! SPDX-License-Identifier: MIT
program test_linalg_hankel
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_hankel, only : hankel
    use scifort_linalg_matmul_hankel, only : matmul_hankel
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none
    real(dp) :: c(4), r(5), h(4,5), x(5), y(4)
    real(dp) :: xx(5,2), yy(4,2), nan_value
    complex(dp) :: zc(4), zr(5), zh(4,5), zx(5), zy(4)
    complex(dp) :: zxx(5,2), zyy(4,2)
    integer :: status, n, m

    c = [1.0_dp,2.0_dp,3.0_dp,4.0_dp]
    r = [-77.0_dp,7.0_dp,7.0_dp,8.0_dp,9.0_dp]
    call hankel(c,h,status,r)
    call success(status,'SciPy documented rectangular Hankel example')
    if (any(h(1,:) /= [1.0_dp,2.0_dp,3.0_dp,4.0_dp,7.0_dp]) .or. &
        any(h(4,:) /= [4.0_dp,7.0_dp,7.0_dp,8.0_dp,9.0_dp])) &
        error stop 'Hankel matrix contents'
    x = [2.0_dp,-1.0_dp,0.5_dp,1.0_dp,3.0_dp]
    call matmul_hankel(c,x,y,status,r)
    call success(status,'real vector multiplication')
    if (maxval(abs(y-matmul(h,x))) > 3.0e-13_dp) error stop 'real vector product'
    xx(:,1) = x
    xx(:,2) = 1.5_dp*x-2.0_dp
    call matmul_hankel(c,xx,yy,status,r)
    call success(status,'real matrix multiplication')
    if (maxval(abs(yy-matmul(h,xx))) > 5.0e-13_dp) error stop 'real matrix product'

    zc = cmplx(c,[0.3_dp,-0.2_dp,0.0_dp,1.0_dp],dp)
    zr = cmplx(r,[-100.0_dp,0.1_dp,0.5_dp,-0.8_dp,0.9_dp],dp)
    zx = cmplx(x,[0.1_dp,2.0_dp,0.0_dp,-1.0_dp,0.3_dp],dp)
    call hankel(zc,zh,status,zr)
    call success(status,'complex Hankel')
    call matmul_hankel(zc,zx,zy,status,zr)
    call success(status,'complex product')
    if (maxval(abs(zy-matmul(zh,zx))) > 3.0e-12_dp) error stop 'complex matrix product'
    zxx(:,1) = zx
    zxx(:,2) = 2.0_dp*zx+cmplx(0.5_dp,-0.1_dp,dp)
    call matmul_hankel(zc,zxx,zyy,status,zr)
    call success(status,'complex matrix product')
    if (maxval(abs(zyy-matmul(zh,zxx))) > 3.0e-12_dp) error stop 'complex RHS matrices'
    call matmul_hankel(c,zx,zy,status,r)
    call success(status,'real coefficients and complex vector')
    if (maxval(abs(zy-matmul(cmplx(h,kind=dp),zx))) > 2.0e-12_dp) &
        error stop 'mixed real/complex vector'
    call matmul_hankel(zc,x,zy,status,zr)
    call success(status,'complex coefficients and real vector')
    if (maxval(abs(zy-matmul(zh,cmplx(x,kind=dp)))) > 2.0e-12_dp) &
        error stop 'mixed complex/real vector'
    call matmul_hankel(c,zxx,zyy,status,r)
    call success(status,'real coefficients and complex matrix')
    if (maxval(abs(zyy-matmul(cmplx(h,kind=dp),zxx))) > 2.0e-12_dp) &
        error stop 'mixed real/complex matrix'
    call matmul_hankel(zc,xx,zyy,status,zr)
    call success(status,'complex coefficients and real matrix')
    if (maxval(abs(zyy-matmul(zh,cmplx(xx,kind=dp)))) > 2.0e-12_dp) &
        error stop 'mixed complex/real matrix'

    ! Complex matrix construction with a real or complex last row.
    call hankel(c,zh,status,zr)
    call success(status,'real first column with complex last row')
    if (any(zh(:,1) /= cmplx(c,kind=dp)) .or. &
        any(zh(4,2:) /= zr(2:))) error stop 'mixed complex row constructor'
    call hankel(zc,zh,status,r)
    call success(status,'complex first column with real last row')
    if (any(zh(:,1) /= zc) .or. &
        any(zh(4,2:) /= cmplx(r(2:),kind=dp))) error stop 'mixed real row constructor'

    block
        real(dp) :: c3(3), h3(3,3), xv(3), yv(3)
        complex(dp) :: cz3(3), hz3(3,3), xz3(3), yz3(3)
        c3 = [1.0_dp,17.0_dp,99.0_dp]
        call hankel(c3,h3,status)
        call success(status,'default row')
        if (any(h3(:,1) /= c3) .or. any(h3(3,2:) /= 0.0_dp)) &
            error stop 'default row zero convention'
        xv = [1.0_dp,2.0_dp,3.0_dp]
        call matmul_hankel(c3,xv,yv,status)
        call success(status,'real default last row')
        if (maxval(abs(yv-matmul(h3,xv))) > 3.0e-13_dp) &
            error stop 'real default row multiplication'
        cz3 = cmplx(c3,[3.0_dp,-1.0_dp,0.5_dp],dp)
        xz3 = cmplx(xv,[1.0_dp,-1.0_dp,1.0_dp],dp)
        call hankel(cz3,hz3,status)
        call success(status,'complex default row')
        if (any(hz3(3,2:) /= cmplx(0.0_dp,0.0_dp,dp))) &
            error stop 'complex default row convention'
        call matmul_hankel(cz3,xz3,yz3,status)
        call success(status,'complex default last row')
        if (maxval(abs(yz3-matmul(hz3,xz3))) > 5.0e-13_dp) &
            error stop 'complex default row multiplication'
    end block

    ! Shape errors and deterministic zero outputs.
    block
        real(dp) :: bad(2), badmat(2,5), shortx(3), bady(4,2)
        bad = 99.0_dp
        call matmul_hankel(c,x,bad,status,r)
        if (status /= linalg_status_invalid_shape) error stop 'short output accepted'
        if (any(bad /= 0.0_dp)) error stop 'short output not zeroed'
        shortx = 0.0_dp
        call matmul_hankel(c,shortx,y,status,r)
        if (status /= linalg_status_invalid_shape) error stop 'short multiplier accepted'
        if (any(y /= 0.0_dp)) error stop 'error output not zeroed'
        call hankel(c,badmat,status,r)
        if (status /= linalg_status_invalid_shape) error stop 'bad Hankel shape accepted'
        call matmul_hankel(c,xx,bady,status,r)
        if (status /= linalg_status_success) error stop 'matrix shape valid'
    end block

    nan_value = ieee_value(0.0_dp,ieee_quiet_nan)
    block
        real(dp) :: badrow(5), coeff(4)
        badrow = r
        badrow(1) = nan_value
        call hankel(c,h,status,badrow)
        if (status /= linalg_status_no_convergence) error stop 'nan row accepted by constructor'
        call matmul_hankel(c,x,y,status,badrow)
        if (status /= linalg_status_no_convergence) error stop 'nan row accepted by product'
        if (any(y /= 0.0_dp)) error stop 'nan row produced nonzero result'
        coeff = c
        coeff(2) = nan_value
        call matmul_hankel(coeff,x,y,status,r)
        if (status /= linalg_status_no_convergence) error stop 'nan coefficients accepted'
        ! r(1) is ignored by Hankel indexing. With checks disabled, even
        ! a nonfinite ignored entry must not contaminate the computation.
        call matmul_hankel(c,x,y,status,badrow,check_finite=.false.)
        call success(status,'ignored NaN r(1) with finite check disabled')
        call hankel(c,h,status,r)
        call success(status,'fresh dense Hankel for ignored-row test')
        if (maxval(abs(y-matmul(h,x))) > 4.0e-13_dp) &
            error stop 'ignored NaN row contaminated the result'
        block
            real(dp) :: badx(5)
            badx = x
            badx(3) = nan_value
            call matmul_hankel(c,badx,y,status,r)
            if (status /= linalg_status_no_convergence) error stop 'nan RHS accepted'
            if (any(y /= 0.0_dp)) error stop 'nan RHS output not zeroed'
        end block
    end block

    block
        real(dp) :: e0(0), empty_out(0), r2(2), y2(2), h0(0,2), h20(2,0), xx0(0,3), yy0(2,3)
        complex(dp) :: ez0(0), zyy0(2,3)
        r2 = [0.0_dp,1.0_dp]
        call hankel(e0,h0,status,r2)
        call success(status,'zero rows in Hankel')
        call hankel(r2,h20,status,e0)
        call success(status,'zero columns in Hankel')
        call matmul_hankel(r2,e0,y2,status,e0)
        call success(status,'zero column product')
        if (any(y2 /= 0.0_dp)) error stop 'empty multiplier gave nonzero output'
        call matmul_hankel(e0,r2,empty_out,status,r2)
        call success(status,'zero row product')
        call matmul_hankel(r2,xx0,yy0,status,e0)
        call success(status,'zero RHS elements')
        if (any(yy0 /= 0.0_dp)) error stop 'zero RHS product not zero'
        call matmul_hankel(cmplx(r2,kind=dp),ez0,zyy0(:,1),status,ez0)
        call success(status,'zero complex multiplier')
    end block

    do n = 1, 130
        m = 1 + mod(19*n,67)
        call exercise_real(n,m)
        call exercise_complex(n,m)
    end do
    call exercise_scale(1.0e-80_dp)
    call exercise_scale(1.0e-20_dp)
    call exercise_scale(1.0_dp)
    call exercise_scale(1.0e20_dp)
    call exercise_scale(1.0e80_dp)
    ! Large case guards against replacing the structured algorithm by dense multiplication.
    block
        real(dp) :: cc(4096), rr(4096), xxlarge(4096), yylarge(4096), expected(4096)
        integer :: i
        cc = 0.0_dp
        cc(1:2) = [2.0_dp,3.0_dp]
        rr = 0.0_dp
        rr(1) = -100.0_dp
        do i=1,4096
            xxlarge(i) = sin(real(i,dp)*0.02_dp)
        end do
        expected = 0.0_dp
        expected(1) = 2.0_dp*xxlarge(1)+3.0_dp*xxlarge(2)
        expected(2) = 3.0_dp*xxlarge(1)
        call matmul_hankel(cc,xxlarge,yylarge,status,rr)
        call success(status,'4096-element structured FFT Hankel')
        if (maxval(abs(yylarge-expected)) > 1.0e-10_dp) &
            error stop 'large structured Hankel product'
    end block
    print '(a)', 'test_linalg_hankel: PASS'

contains

    subroutine success(s,label)
        integer, intent(in) :: s !! result status
        character(*), intent(in) :: label !! test description
        if (s /= linalg_status_success) then
            print '(a,1x,a,1x,i0)', 'FAILED',label,s
            error stop 'Hankel test failed'
        end if
    end subroutine success

    subroutine exercise_real(nn,mm)
        integer, intent(in) :: nn !! row count
        integer, intent(in) :: mm !! column count
        real(dp) :: cc(nn), rr(mm), aa(nn,mm), xxloc(mm,3), yyloc(nn,3)
        real(dp) :: expected(nn,3)
        integer :: i, j, s
        do i=1,nn
            cc(i) = sin(real(i+3*nn,dp))
        end do
        do i=1,mm
            rr(i) = cos(real(i+17*mm,dp))
        end do
        do j=1,3
            do i=1,mm
                xxloc(i,j) = sin(real(i*5+j*13,dp))
            end do
        end do
        call hankel(cc,aa,s,rr)
        call success(s,'random rectangular real construction')
        expected = matmul(aa,xxloc)
        call matmul_hankel(cc,xxloc,yyloc,s,rr)
        call success(s,'random rectangular real product')
        if (maxval(abs(yyloc-expected)) > 2.0e-11_dp*max(1.0_dp,maxval(abs(expected)))) &
            error stop 'random rectangular real Hankel'
        call matmul_hankel(cc,xxloc(:,1),yyloc(:,1),s,rr)
        call success(s,'random rectangular real vector')
        if (maxval(abs(yyloc(:,1)-expected(:,1))) > &
            2.0e-11_dp*max(1.0_dp,maxval(abs(expected(:,1))))) &
            error stop 'random rectangular real vector Hankel'
    end subroutine exercise_real

    subroutine exercise_complex(nn,mm)
        integer, intent(in) :: nn !! row count
        integer, intent(in) :: mm !! column count
        complex(dp) :: cc(nn), rr(mm), aa(nn,mm), xxloc(mm,2), yyloc(nn,2)
        complex(dp) :: expected(nn,2)
        integer :: i, j, s
        do i=1,nn
            cc(i) = cmplx(sin(real(11*i+nn,dp)),cos(real(3*i+nn,dp)),dp)
        end do
        do i=1,mm
            rr(i) = cmplx(cos(real(19*i+mm,dp)),sin(real(5*i+mm,dp)),dp)
        end do
        do j=1,2
            do i=1,mm
                xxloc(i,j) = cmplx(sin(real(13*i+j,dp)),cos(real(i+17*j,dp)),dp)
            end do
        end do
        call hankel(cc,aa,s,rr)
        call success(s,'random rectangular complex construction')
        expected = matmul(aa,xxloc)
        call matmul_hankel(cc,xxloc,yyloc,s,rr)
        call success(s,'random rectangular complex product')
        if (maxval(abs(yyloc-expected)) > 3.0e-11_dp*max(1.0_dp,maxval(abs(expected)))) &
            error stop 'random rectangular complex Hankel'
    end subroutine exercise_complex

    subroutine exercise_scale(factor)
        real(dp), intent(in) :: factor !! coefficient scaling
        real(dp) :: cc(4), rr(6), aa(4,6), xxloc(6), yyloc(4)
        integer :: s
        cc = factor*[0.2_dp,3.0_dp,2.0_dp,-1.0_dp]
        rr = factor*[40.0_dp,0.2_dp,0.4_dp,-1.0_dp,0.3_dp,1.0_dp]
        xxloc = [1.0_dp,2.0_dp,-1.0_dp,4.0_dp,5.0_dp,-2.0_dp]
        call hankel(cc,aa,s,rr)
        call success(s,'scaled Hankel construction')
        call matmul_hankel(cc,xxloc,yyloc,s,rr)
        call success(s,'scaled Hankel multiplication')
        if (maxval(abs(yyloc-matmul(aa,xxloc))) > &
            4.0e-12_dp*maxval(abs(matmul(aa,xxloc)))) &
            error stop 'scaled Hankel product'
    end subroutine exercise_scale
end program test_linalg_hankel
