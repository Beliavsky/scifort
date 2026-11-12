! SPDX-License-Identifier: MIT
program test_linalg_circulant
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_singular, linalg_status_no_convergence
    use scifort_linalg_circulant, only : circulant, solve_circulant
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none
    integer :: code, n, k
    real(dp) :: c(3), b(3), x(3), a(3,3), expected(3)
    real(dp) :: nan_value, scalar_x(1)
    real(dp) :: c0(0), b0(0), x0(0), a0(0,0)
    real(dp) :: empty_b(3,0), empty_x(3,0)
    real(dp) :: bc(3,2), xc(3,2), xr(2,3)
    complex(dp) :: zc(3), zb(3), zx(3), za(3,3)
    complex(dp) :: zbm(3,2), zxm(3,2), ztm(2,3)

    call circulant([1.0_dp,2.0_dp,3.0_dp],a,code)
    call success(code,'constructor')
    if (any(a /= reshape([1.0_dp,2.0_dp,3.0_dp, &
                           3.0_dp,1.0_dp,2.0_dp, &
                           2.0_dp,3.0_dp,1.0_dp],[3,3]))) &
        error stop 'circulant matrix orientation'

    c = [2.0_dp,2.0_dp,4.0_dp]
    b = [1.0_dp,2.0_dp,3.0_dp]
    expected = [0.75_dp,-0.25_dp,0.25_dp]
    call solve_circulant(c,b,x,code)
    call success(code,'SciPy documented example')
    if (maxval(abs(x-expected)) > 1.e-13_dp) error stop 'documented solution mismatch'

    ! Exact singular case, SciPy pseudoinverse example.
    block
        real(dp) :: cs(4), bs(4), xs(4), expected_s(4)
        cs = [1.0_dp,1.0_dp,0.0_dp,0.0_dp]
        bs = [1.0_dp,2.0_dp,3.0_dp,4.0_dp]
        expected_s = [0.25_dp,1.25_dp,2.25_dp,1.25_dp]
        xs = -99.0_dp
        call solve_circulant(cs,bs,xs,code)
        if (code /= linalg_status_singular .or. any(xs /= 0.0_dp)) &
            error stop 'default singular policy'
        call solve_circulant(cs,bs,xs,code,singular='lstsq')
        call success(code,'singular least squares')
        if (maxval(abs(xs-expected_s)) > 3.e-13_dp) &
            error stop 'minimum norm solution mismatch'
    end block

    c = [2.0_dp,0.2_dp,-0.1_dp]
    bc(:,1) = [1.0_dp,3.0_dp,-2.0_dp]
    bc(:,2) = [2.0_dp,0.0_dp,5.0_dp]
    call circulant(c,a,code)
    call success(code,'constructor for matrix RHS')
    call solve_circulant(c,matmul(a,bc),xc,code)
    call success(code,'matrix RHS')
    if (maxval(abs(xc-bc)) > 1.e-13_dp) error stop 'matrix RHS solution'

    ! Positive and negative axis positions: input n by nrhs and transposed.
    call solve_circulant(c,transpose(matmul(a,bc)),xr,code,baxis=-1,outaxis=-1)
    call success(code,'last-axis input and output')
    if (maxval(abs(xr-transpose(bc))) > 1.e-13_dp) error stop 'transposed axis solution'
    call solve_circulant(c,transpose(matmul(a,bc)),xc,code,baxis=2,outaxis=1)
    call success(code,'last-axis input first-axis output')
    if (maxval(abs(xc-bc)) > 1.e-13_dp) error stop 'axis reordering'
    call solve_circulant(c,matmul(a,bc),xr,code,baxis=1,outaxis=2)
    call success(code,'first-axis input last-axis output')
    if (maxval(abs(xr-transpose(bc))) > 1.e-13_dp) error stop 'axis output transpose'

    call solve_circulant(c,bc,xc,code,baxis=0)
    if (code /= linalg_status_invalid_shape .or. any(xc /= 0.0_dp)) &
        error stop 'invalid axis accepted'
    call solve_circulant(c,bc,xr,code,baxis=1,outaxis=1)
    if (code /= linalg_status_invalid_shape .or. any(xr /= 0.0_dp)) &
        error stop 'invalid output shape accepted'
    call solve_circulant(c,b,x,code,singular='unknown')
    if (code /= linalg_status_invalid_shape .or. any(x /= 0.0_dp)) &
        error stop 'invalid singular policy accepted'
    call solve_circulant(c,b,x,code,tol=-0.1_dp)
    if (code /= linalg_status_invalid_shape .or. any(x /= 0.0_dp)) &
        error stop 'negative tolerance accepted'

    nan_value = ieee_value(0.0_dp,ieee_quiet_nan)
    c = [2.0_dp,nan_value,0.1_dp]
    call solve_circulant(c,b,x,code)
    if (code /= linalg_status_no_convergence .or. any(x /= 0.0_dp)) &
        error stop 'nonfinite input accepted'
    call circulant(c,a,code)
    if (code /= linalg_status_no_convergence .or. any(a /= 0.0_dp)) &
        error stop 'nonfinite constructor accepted'

    call circulant(c0,a0,code)
    call success(code,'empty matrix')
    call solve_circulant(c0,b0,x0,code)
    call success(code,'empty vector solve')
    call solve_circulant([1.0_dp,0.1_dp,0.1_dp],empty_b,empty_x,code)
    call success(code,'zero RHS columns')
    zc = [cmplx(3.0_dp,0.2_dp,dp), cmplx(-0.1_dp,0.4_dp,dp), &
          cmplx(0.2_dp,-0.3_dp,dp)]
    call circulant(zc,za,code)
    call success(code,'complex constructor')
    if (maxval(abs(za(:,1)-zc)) > 0.0_dp) error stop 'complex first column'
    if (abs(za(1,2)-zc(3)) > 0.0_dp) error stop 'complex nonconjugate winding'
    zb = matmul(za,[cmplx(1.0_dp,0.5_dp,dp),cmplx(-0.5_dp,0.7_dp,dp), &
                  cmplx(0.9_dp,-0.3_dp,dp)])
    call solve_circulant(zc,zb,zx,code)
    call success(code,'complex vector solve')
    if (maxval(abs(matmul(za,zx)-zb)) > 2.e-13_dp) error stop 'complex vector residual'

    call solve_circulant(zc,real(zb,dp),zx,code)
    call success(code,'complex coefficients real RHS')
    if (maxval(abs(matmul(za,zx)-real(zb,dp))) > 2.e-13_dp) &
        error stop 'complex-real residual'
    call solve_circulant(real(zc,dp),zb,zx,code)
    call success(code,'real coefficients complex RHS')
    block
        real(dp) :: ar(3,3)
        call circulant(real(zc,dp),ar,code)
        call success(code,'real circulant for mixed solve')
        if (maxval(abs(matmul(cmplx(ar,kind=dp),zx)-zb)) > 2.e-13_dp) &
            error stop 'real-complex residual'
    end block

    zbm(:,1) = zb
    zbm(:,2) = cmplx(bc(:,2),kind=dp)
    call solve_circulant(zc,zbm,zxm,code)
    call success(code,'complex multiple RHS')
    if (maxval(abs(matmul(za,zxm)-zbm)) > 2.e-13_dp) &
        error stop 'complex multiple RHS residual'
    call solve_circulant(zc,real(zbm,dp),zxm,code)
    call success(code,'complex-real matrix RHS')
    if (maxval(abs(matmul(za,zxm)-real(zbm,dp))) > 2.e-13_dp) &
        error stop 'complex-real matrix residual'
    call solve_circulant(real(zc,dp),zbm,ztm,code,outaxis=2)
    call success(code,'real-complex matrix last outaxis')
    block
        real(dp) :: ar(3,3)
        call circulant(real(zc,dp),ar,code)
        call success(code,'real matrix for mixed solve')
        if (maxval(abs(matmul(cmplx(ar,kind=dp),transpose(ztm))-zbm)) > 2.e-13_dp) &
            error stop 'real-complex matrix residual'
    end block

    do n = 1, 140
        call exercise_real(n)
        call exercise_complex(n)
    end do
    do k = -80, 80, 80
        call exercise_scaled(10.0_dp**k)
    end do

    print '(a)', 'test_linalg_circulant: PASS'

contains
    subroutine success(s,label)
        integer, intent(in) :: s !! actual status
        character(*), intent(in) :: label !! test description
        if (s /= linalg_status_success) then
            print '(a,1x,a,1x,i0)', 'FAILED',label,s
            error stop 'circulant operation failed'
        end if
    end subroutine success

    subroutine exercise_real(nn)
        integer, intent(in) :: nn !! circulant order
        real(dp) :: cc(nn), aa(nn,nn), yy(nn,3), bb(nn,3), xx(nn,3)
        real(dp) :: vx(nn), sx(nn)
        integer :: i, j, status

        cc(1) = real(nn+3,dp)
        do i=2,nn
            cc(i) = 0.23_dp*sin(real(3*i+11*nn,dp))
        end do
        call circulant(cc,aa,status)
        call success(status,'randomized real constructor')
        do j=1,3
            do i=1,nn
                yy(i,j) = cos(real(7*i+13*j+nn,dp))
            end do
        end do
        bb = matmul(aa,yy)
        call solve_circulant(cc,bb,xx,status)
        call success(status,'randomized real solve')
        if (maxval(abs(xx-yy)) > 1.3e-11_dp) error stop 'randomized real solution'
        call solve_circulant(cc,bb(:,1),vx,status)
        call success(status,'randomized real vector')
        if (maxval(abs(vx-yy(:,1))) > 1.3e-11_dp) error stop 'randomized real vector solution'
        sx = yy(:,1)
        if (maxval(abs(matmul(aa,sx)-bb(:,1))) > 1.3e-11_dp) &
            error stop 'randomized test preparation'
    end subroutine exercise_real

    subroutine exercise_complex(nn)
        integer, intent(in) :: nn !! circulant order
        complex(dp) :: cc(nn), aa(nn,nn), yy(nn,2), bb(nn,2), xx(nn,2), vv(nn)
        integer :: i, j, status

        cc(1) = cmplx(real(nn+3,dp),0.2_dp,dp)
        do i=2,nn
            cc(i) = cmplx(0.1_dp*sin(real(i+nn,dp)),0.13_dp*cos(real(2*i+nn,dp)),dp)
        end do
        call circulant(cc,aa,status)
        call success(status,'randomized complex constructor')
        do j=1,2
            do i=1,nn
                yy(i,j) = cmplx(cos(real(i+j,dp)),sin(real(3*i-j,dp)),dp)
            end do
        end do
        bb = matmul(aa,yy)
        call solve_circulant(cc,bb,xx,status)
        call success(status,'randomized complex solve')
        if (maxval(abs(xx-yy)) > 1.6e-11_dp) error stop 'randomized complex solution'
        call solve_circulant(cc,bb(:,1),vv,status)
        call success(status,'randomized complex vector')
        if (maxval(abs(vv-yy(:,1))) > 1.6e-11_dp) error stop 'randomized complex vector solution'
    end subroutine exercise_complex

    subroutine exercise_scaled(factor)
        real(dp), intent(in) :: factor !! positive scale factor
        real(dp) :: cc(7), aa(7,7), yy(7), bb(7), xx(7)
        integer :: status

        cc = factor*[4.0_dp,0.2_dp,-0.1_dp,0.3_dp,0.02_dp,-0.05_dp,-0.2_dp]
        yy = [1.0_dp,0.2_dp,-0.1_dp,0.4_dp,-0.7_dp,0.9_dp,1.2_dp]
        call circulant(cc,aa,status)
        call success(status,'scaled constructor')
        bb = matmul(aa,yy)
        call solve_circulant(cc,bb,xx,status)
        call success(status,'scaled circulant solve')
        if (maxval(abs(xx-yy)) > 2.e-12_dp) error stop 'scaled solution'
    end subroutine exercise_scaled
end program test_linalg_circulant
