! SPDX-License-Identifier: MIT
program test_linalg_matmul_toeplitz
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_toeplitz, only : toeplitz
    use scifort_linalg_matmul_toeplitz, only : matmul_toeplitz
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none
    real(dp) :: column(3)
    real(dp) :: first_row(4)
    real(dp) :: dense(3,4)
    real(dp) :: x(4)
    real(dp) :: y(3)
    real(dp) :: rhs(4,2)
    real(dp) :: products(3,2)
    real(dp) :: nan_value
    complex(dp) :: zcolumn(3)
    complex(dp) :: zrow(4)
    complex(dp) :: zdense(3,4)
    complex(dp) :: zx(4)
    complex(dp) :: zy(3)
    complex(dp) :: zrhs(4,2)
    complex(dp) :: zproducts(3,2)
    integer :: status
    integer :: n
    integer :: m

    ! Example from the SciPy toeplitz documentation: a 3x4 matrix.
    column = [1.0_dp,2.0_dp,3.0_dp]
    first_row = [-999.0_dp,4.0_dp,5.0_dp,6.0_dp]
    call toeplitz(column,dense,status,first_row)
    call success(status,'real rectangular construction')
    x = [1.0_dp,2.0_dp,-1.0_dp,3.0_dp]
    call matmul_toeplitz(column,x,y,status,first_row)
    call success(status,'real rectangular vector')
    if (maxval(abs(y-matmul(dense,x))) > 3.0e-13_dp) &
        error stop 'real rectangular product wrong'
    rhs(:,1) = x
    rhs(:,2) = [-2.0_dp,0.2_dp,1.0_dp,4.0_dp]
    call matmul_toeplitz(column,rhs,products,status,first_row)
    call success(status,'real rectangular matrix')
    if (maxval(abs(products-matmul(dense,rhs))) > 3.0e-13_dp) &
        error stop 'real matrix multiplication wrong'

    zcolumn = [cmplx(3.0_dp,4.0_dp,dp), cmplx(0.5_dp,-0.2_dp,dp), &
        cmplx(-1.0_dp,0.7_dp,dp)]
    zrow = [cmplx(-100.0_dp,21.0_dp,dp), cmplx(1.0_dp,3.0_dp,dp), &
        cmplx(-0.9_dp,0.6_dp,dp),cmplx(0.3_dp,0.1_dp,dp)]
    zx = [cmplx(-1.0_dp,0.1_dp,dp), cmplx(0.5_dp,1.0_dp,dp), &
        cmplx(2.0_dp,-0.5_dp,dp),cmplx(0.3_dp,0.9_dp,dp)]
    call toeplitz(zcolumn,zdense,status,zrow)
    call success(status,'complex rectangular construction')
    call matmul_toeplitz(zcolumn,zx,zy,status,zrow)
    call success(status,'complex rectangular vector')
    if (maxval(abs(zy-matmul(zdense,zx))) > 4.0e-13_dp) &
        error stop 'complex rectangular product wrong'
    zrhs(:,1) = zx
    zrhs(:,2) = 2.0_dp*zx+cmplx(0.1_dp,-0.6_dp,dp)
    call matmul_toeplitz(zcolumn,zrhs,zproducts,status,zrow)
    call success(status,'complex rectangular matrix')
    if (maxval(abs(zproducts-matmul(zdense,zrhs))) > 8.0e-13_dp) &
        error stop 'complex matrix multiplication wrong'

    ! Default complex first row is conjugate(column), not column itself.
    block
        complex(dp) :: square(3,3)
        complex(dp) :: xs(3)
        complex(dp) :: ys(3)
        call toeplitz(zcolumn,square,status)
        call success(status,'Hermitian constructor')
        xs = [cmplx(0.5_dp,0.0_dp,dp),cmplx(-2.0_dp,0.4_dp,dp), &
            cmplx(1.0_dp,-3.0_dp,dp)]
        call matmul_toeplitz(zcolumn,xs,ys,status)
        call success(status,'default conjugate row')
        if (maxval(abs(ys-matmul(square,xs))) > 5.0e-13_dp) &
            error stop 'default complex first row is not conjugated'
    end block

    call matmul_toeplitz(column,zx,zy,status,first_row)
    call success(status,'real coefficients and complex vector')
    if (maxval(abs(zy-matmul(cmplx(dense,kind=dp),zx))) > 4.0e-13_dp) &
        error stop 'real/complex vector product wrong'
    call matmul_toeplitz(column,zrhs,zproducts,status,first_row)
    call success(status,'real coefficients and complex matrix')
    if (maxval(abs(zproducts-matmul(cmplx(dense,kind=dp),zrhs))) > 8.0e-13_dp) &
        error stop 'real/complex matrix product wrong'
    call matmul_toeplitz(zcolumn,x,zy,status,zrow)
    call success(status,'complex coefficients and real vector')
    if (maxval(abs(zy-matmul(zdense,x))) > 4.0e-13_dp) &
        error stop 'complex/real vector product wrong'
    call matmul_toeplitz(zcolumn,rhs,zproducts,status,zrow)
    call success(status,'complex coefficients and real matrix')
    if (maxval(abs(zproducts-matmul(zdense,rhs))) > 8.0e-13_dp) &
        error stop 'complex/real matrix product wrong'

    ! Short and long input vectors must produce a status, not out-of-bounds.
    call matmul_toeplitz(column,x(:3),y,status,first_row)
    if (status /= linalg_status_invalid_shape .or. any(y /= 0.0_dp)) &
        error stop 'incorrect multiplier length accepted'
    call matmul_toeplitz(column,x,y(:2),status,first_row)
    if (status /= linalg_status_invalid_shape .or. any(y(:2) /= 0.0_dp)) &
        error stop 'incorrect output length accepted'
    call matmul_toeplitz(column,rhs,products(:,1:1),status,first_row)
    if (status /= linalg_status_invalid_shape .or. any(products(:,1:1) /= 0.0_dp)) &
        error stop 'incorrect output RHS count accepted'
    nan_value = ieee_value(0.0_dp,ieee_quiet_nan)
    x(2) = nan_value
    call matmul_toeplitz(column,x,y,status,first_row,check_finite=.true.)
    if (status /= linalg_status_no_convergence .or. any(y /= 0.0_dp)) &
        error stop 'nonfinite multiplier accepted'
    x(2) = 2.0_dp
    first_row(1) = nan_value
    call matmul_toeplitz(column,x,y,status,first_row)
    if (status /= linalg_status_no_convergence .or. any(y /= 0.0_dp)) &
        error stop 'nonfinite ignored row entry accepted with input checking'
    first_row(1) = 1.0_dp

    block
        real(dp) :: c0(0), r0(0), x0(0), y0(0)
        real(dp) :: out20(2), zeros(0), zero_rhs(2,0), zero_product(3,0)
        call matmul_toeplitz(c0,x0,y0,status)
        call success(status,'empty default Toeplitz')
        call matmul_toeplitz([1.0_dp,2.0_dp],zeros,out20,status,r0)
        call success(status,'nonzero output length with empty multiplier')
        if (any(out20 /= 0.0_dp)) error stop 'empty input product nonzero'
        call matmul_toeplitz(column,zero_rhs,zero_product,status,r=[1.0_dp,2.0_dp])
        call success(status,'zero right-hand side columns')
        call matmul_toeplitz(c0,[1.0_dp,2.0_dp],y0,status,r=[1.0_dp,3.0_dp])
        call success(status,'zero row count with nonzero column count')
    end block

    ! Exercise all FFT backends with both rectangular and square matrices.
    do n = 1, 128
        m = 1+mod(17*n,65)
        call exercise_real(n,m)
        call exercise_complex(n,m)
    end do
    do n = -80, 80, 80
        call exercise_scale(10.0_dp**n)
    end do
    ! A moderately large prime-length embedding would make accidental
    ! dense O(n^2) construction expensive and checks Bluestein transform.
    block
        real(dp) :: big_c(4096), big_r(4096), big_x(4096), big_y(4096)
        real(dp) :: expected(4096)
        integer :: i
        big_c = 0.0_dp
        big_c(1:2) = [2.0_dp,-1.0_dp]
        big_r = 0.0_dp
        big_r(1) = 2.0_dp
        do i=1,4096
            big_x(i) = cos(real(i,dp)*0.1_dp)
        end do
        expected = 2.0_dp*big_x
        expected(2:) = expected(2:)-big_x(:4095)
        call matmul_toeplitz(big_c,big_x,big_y,status,r=big_r)
        call success(status,'4096-element FFT multiplication')
        if (maxval(abs(big_y-expected)) > 1.0e-10_dp) &
            error stop 'large structured Toeplitz product'
    end block
    print '(a)', 'test_linalg_matmul_toeplitz: PASS'

contains

    subroutine success(s, label)
        integer, intent(in) :: s !! status returned by SciFort
        character(*), intent(in) :: label !! test description
        if (s /= linalg_status_success) then
            print '(a,1x,a,1x,i0)', 'FAILED',label,s
            error stop 'Toeplitz multiplication failed'
        end if
    end subroutine success

    subroutine exercise_real(nn,mm)
        integer, intent(in) :: nn !! number of rows
        integer, intent(in) :: mm !! number of columns
        real(dp) :: cc(nn), rr(mm), aa(nn,mm), xx(mm,3), yy(nn,3)
        real(dp) :: expected(nn,3)
        integer :: i, j, s

        do i=1,nn
            cc(i) = 0.5_dp*sin(real(17*i+nn,dp))
        end do
        do i=1,mm
            rr(i) = cos(real(11*i+mm,dp))
        end do
        do j=1,3
            do i=1,mm
                xx(i,j) = sin(real(i+13*j+nn,dp))
            end do
        end do
        call toeplitz(cc,aa,s,rr)
        call success(s,'randomized real Toeplitz')
        expected = matmul(aa,xx)
        call matmul_toeplitz(cc,xx,yy,s,rr)
        call success(s,'randomized rectangular real')
        if (maxval(abs(yy-expected)) > 2.0e-11_dp*max(1.0_dp,maxval(abs(expected)))) &
            error stop 'randomized real matrix product'
        call matmul_toeplitz(cc,xx(:,1),yy(:,1),s,rr)
        call success(s,'randomized real vector')
        if (maxval(abs(yy(:,1)-expected(:,1))) > 2.0e-11_dp*max(1.0_dp,maxval(abs(expected(:,1))))) &
            error stop 'randomized real vector product'
    end subroutine exercise_real

    subroutine exercise_complex(nn,mm)
        integer, intent(in) :: nn !! number of rows
        integer, intent(in) :: mm !! number of columns
        complex(dp) :: cc(nn), rr(mm), aa(nn,mm), xx(mm,2), yy(nn,2)
        complex(dp) :: expected(nn,2)
        integer :: i, j, s

        do i=1,nn
            cc(i) = cmplx(sin(real(i+nn,dp)),cos(real(2*i+nn,dp)),dp)
        end do
        do i=1,mm
            rr(i) = cmplx(cos(real(3*i+mm,dp)),sin(real(11*i+mm,dp)),dp)
        end do
        do j=1,2
            do i=1,mm
                xx(i,j) = cmplx(cos(real(i+13*j+nn,dp)),sin(real(17*i+j,dp)),dp)
            end do
        end do
        call toeplitz(cc,aa,s,rr)
        call success(s,'randomized complex Toeplitz')
        expected = matmul(aa,xx)
        call matmul_toeplitz(cc,xx,yy,s,rr)
        call success(s,'randomized rectangular complex')
        if (maxval(abs(yy-expected)) > 3.0e-11_dp*max(1.0_dp,maxval(abs(expected)))) &
            error stop 'randomized complex matrix product'
    end subroutine exercise_complex

    subroutine exercise_scale(factor)
        real(dp), intent(in) :: factor !! scale applied to the Toeplitz coefficients
        real(dp) :: cc(4), rr(5), aa(4,5), xx(5), yy(4)
        integer :: s

        cc = factor*[4.0_dp,0.2_dp,-0.1_dp,0.3_dp]
        rr = factor*[10.0_dp,0.1_dp,0.3_dp,-0.2_dp,0.07_dp]
        xx = [1.0_dp,-0.2_dp,0.4_dp,0.5_dp,1.2_dp]
        call toeplitz(cc,aa,s,rr)
        call success(s,'scaled Toeplitz matrix')
        call matmul_toeplitz(cc,xx,yy,s,rr)
        call success(s,'scaled Toeplitz product')
        if (maxval(abs(yy-matmul(aa,xx))) > &
            3.0e-12_dp*maxval(abs(matmul(aa,xx)))) error stop 'scaled Toeplitz product'
    end subroutine exercise_scale
end program test_linalg_matmul_toeplitz
