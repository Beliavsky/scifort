! SPDX-License-Identifier: MIT
program test_linalg_convolution_matrix
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_convolution_matrix, only : convolution_matrix, convolution_matrix_rows
    use scifort_linalg_matmul_convolution_matrix, only : matmul_convolution_matrix
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_nan
    implicit none
    real(dp) :: a(3), matrix(5,3), v(3), y(5), nan_value
    real(dp) :: matrix_small(2,3), y_bad(4)
    complex(dp) :: z(3), zmat(5,3), zv(3), zy(5)
    integer :: s, i, m, n, mode_index
    character(5), parameter :: modes(3) = ['full ', 'same ', 'valid']

    ! SciPy's documented convolution matrix example.
    a = [-1.0_dp,4.0_dp,-2.0_dp]
    call convolution_matrix(a,matrix,s,mode='full')
    call success(s,'full documented example')
    if (any(matrix(:,1) /= [-1.0_dp,4.0_dp,-2.0_dp,0.0_dp,0.0_dp])) &
        error stop 'full first column'
    if (any(matrix(:,3) /= [0.0_dp,0.0_dp,-1.0_dp,4.0_dp,-2.0_dp])) &
        error stop 'full final column'
    v = [1.0_dp,2.0_dp,3.0_dp]
    call matmul_convolution_matrix(a,v,y,s)
    call success(s,'default full multiplication')
    if (maxval(abs(y-matmul(matrix,v))) > 1.0e-12_dp) error stop 'full product'

    ! Same mode, including even-length shift, and m > n (short input).
    block
        real(dp) :: c(4), mm(4,2), xx(2), yy(4)
        c = [1.0_dp,2.0_dp,3.0_dp,4.0_dp]
        xx = [5.0_dp,7.0_dp]
        call convolution_matrix(c,mm,s,'same')
        call success(s,'same mode m > n')
        if (any(mm(:,1) /= c) .or. any(mm(:,2) /= [0.0_dp,1.0_dp,2.0_dp,3.0_dp])) &
            error stop 'same m > n offset'
        call matmul_convolution_matrix(c,xx,yy,s,'same')
        call success(s,'same mode m > n multiplication')
        if (maxval(abs(yy-matmul(mm,xx))) > 1.0e-12_dp) error stop 'same m > n product'
    end block

    z = cmplx(a, [1.0_dp,-2.0_dp,0.5_dp],dp)
    zv = cmplx(v,[0.5_dp,2.0_dp,-0.75_dp],dp)
    call convolution_matrix(z,zmat,s)
    call success(s,'complex constructor')
    if (maxval(abs(zmat(:,1)-cmplx(matrix(:,1),kind=dp) - &
        cmplx(0.0_dp,[1.0_dp,-2.0_dp,0.5_dp,0.0_dp,0.0_dp],dp))) > 1.0e-12_dp) &
        error stop 'complex convolution placement'
    call matmul_convolution_matrix(z,zv,zy,s)
    call success(s,'complex product')
    if (maxval(abs(zy-matmul(zmat,zv))) > 1.0e-12_dp) error stop 'complex product mismatch'

    if (convolution_matrix_rows(0,4,status=s) /= 0 .or. s /= linalg_status_invalid_shape) &
        error stop 'empty coefficient vector accepted'
    if (convolution_matrix_rows(3,0,status=s) /= 0 .or. s /= linalg_status_invalid_shape) &
        error stop 'zero columns accepted'
    if (convolution_matrix_rows(3,3,'FULL',s) /= 0 .or. s /= linalg_status_invalid_shape) &
        error stop 'bad mode accepted'
    matrix_small = 8.0_dp
    call convolution_matrix(a,matrix_small,s)
    if (s /= linalg_status_invalid_shape .or. any(matrix_small /= 0.0_dp)) &
        error stop 'wrong output matrix shape'
    y_bad = 9.0_dp
    call matmul_convolution_matrix(a,v,y_bad,s)
    if (s /= linalg_status_invalid_shape .or. any(y_bad /= 0.0_dp)) &
        error stop 'wrong output product shape'
    call matmul_convolution_matrix(a,v,y,s,'FULL')
    if (s /= linalg_status_invalid_shape .or. any(y /= 0.0_dp)) &
        error stop 'invalid mode product'
    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    a(2) = nan_value
    call convolution_matrix(a,matrix,s)
    call success(s,'constructor preserves SciPy NaNs')
    if (.not. ieee_is_nan(matrix(2,1))) error stop 'constructor did not preserve NaN'
    call matmul_convolution_matrix(a,v,y,s)
    if (s /= linalg_status_no_convergence .or. any(y /= 0.0_dp)) &
        error stop 'nonfinite FFT product accepted'

    block
        real(dp) :: empty_a(0), invalid_matrix(2,3), empty_x(0), invalid_y(0)
        real(dp) :: batched(2,3), bad_batch(2,4,3), empty_result(5,0), no_rhs(3,0)
        integer :: int_a(2,3), int_batch(2,5,3)
        empty_a = 0.0_dp
        call convolution_matrix(empty_a,invalid_matrix,s)
        if (s /= linalg_status_invalid_shape .or. any(invalid_matrix /= 0.0_dp)) &
            error stop 'empty kernel construction'
        call matmul_convolution_matrix(empty_a,v,y,s)
        if (s /= linalg_status_invalid_shape .or. any(y /= 0.0_dp)) &
            error stop 'empty kernel multiplication'
        call matmul_convolution_matrix([1.0_dp],empty_x,invalid_y,s)
        if (s /= linalg_status_invalid_shape) error stop 'zero-length input accepted'
        batched = 1.0_dp
        call convolution_matrix(batched,bad_batch,s)
        if (s /= linalg_status_invalid_shape .or. any(bad_batch /= 0.0_dp)) &
            error stop 'invalid batched shape accepted'
        int_a(1,:) = [1,2,3]
        int_a(2,:) = [4,5,6]
        call convolution_matrix(int_a,int_batch,s)
        call success(s,'integer batch constructor')
        if (any(int_batch(1,1:3,1) /= int_a(1,:)) .or. &
            any(int_batch(2,1:3,1) /= int_a(2,:))) &
            error stop 'integer batch contents'
        call matmul_convolution_matrix([1.0_dp,2.0_dp,3.0_dp],no_rhs,empty_result,s)
        call success(s,'zero RHS columns')
    end block

    do m=1,38
        do n=1,35
            do mode_index=1,3
                call case_real(m,n,trim(modes(mode_index)))
                call case_complex(m,n,trim(modes(mode_index)))
            end do
        end do
    end do
    call case_scale(1.0e-80_dp)
    call case_scale(1.0e-20_dp)
    call case_scale(1.0_dp)
    call case_scale(1.0e20_dp)
    call case_scale(1.0e80_dp)

    ! FFT should handle thousands of inputs without forming a dense matrix.
    block
        real(dp) :: coeff(4), x(4096), got(4099), expected(4099)
        coeff = [2.0_dp,3.0_dp,0.0_dp,0.0_dp]
        do i=1,4096
            x(i) = cos(real(i,dp)*0.05_dp)
        end do
        expected = 0.0_dp
        do i=1,4096
            expected(i) = expected(i) + 2.0_dp*x(i)
            expected(i+1) = expected(i+1) + 3.0_dp*x(i)
        end do
        call matmul_convolution_matrix(coeff,x,got,s)
        call success(s,'large FFT-based convolution')
        if (maxval(abs(got-expected)) > 3.0e-11_dp) error stop 'large FFT convolution product'
    end block
    print '(a)', 'test_linalg_convolution_matrix: PASS'
contains
    subroutine success(status,label)
        integer, intent(in) :: status !! linalg status
        character(*), intent(in) :: label !! diagnostic label
        if (status /= linalg_status_success) then
            print '(a,1x,a,1x,i0)', 'FAILED',label,status
            error stop 'convolution_matrix test failed'
        end if
    end subroutine success

    subroutine case_real(m,n,mode)
        integer, intent(in) :: m !! coefficient count
        integer, intent(in) :: n !! input length
        character(*), intent(in) :: mode !! convolution mode
        real(dp) :: c(m), x(n), x2(n,3)
        real(dp), allocatable :: aa(:,:), yy(:), y2(:,:), batched(:,:,:)
        real(dp) :: bc(2,m)
        integer, allocatable :: ic(:), imatrix(:,:)
        integer :: k, i, j, q, s
        real(dp) :: norm

        k = convolution_matrix_rows(m,n,mode,s)
        call success(s,'real case geometry')
        allocate(aa(k,n),yy(k),y2(k,3),batched(2,k,n))
        do i=1,m
            c(i) = cos(real(i*7+m,dp))
        end do
        do i=1,n
            x(i) = sin(real(i*3+n,dp))
        end do
        do j=1,3
            do i=1,n
                x2(i,j) = x(i)*real(j,dp)+0.2_dp*j
            end do
        end do
        call convolution_matrix(c,aa,s,mode)
        call success(s,'real constructor')
        do j=1,n
            do i=1,k
                q = i-j+1
                select case(mode)
                case('same')
                    q = q+(min(m,n)-1)/2
                case('valid')
                    q = q+min(m,n)-1
                end select
                if (q >= 1 .and. q <= m) then
                    if (aa(i,j) /= c(q)) error stop 'real convolution location'
                else
                    if (aa(i,j) /= 0.0_dp) error stop 'unexpected nonzero convolution entry'
                end if
            end do
        end do
        call matmul_convolution_matrix(c,x,yy,s,mode)
        call success(s,'real convolution vector')
        norm = max(1.0_dp,maxval(abs(matmul(aa,x))))
        if (maxval(abs(yy-matmul(aa,x))) > 4.0e-11_dp*norm) &
            error stop 'real convolution vector product'
        call matmul_convolution_matrix(c,x2,y2,s,mode)
        call success(s,'real convolution matrix')
        norm = max(1.0_dp,maxval(abs(matmul(aa,x2))))
        if (maxval(abs(y2-matmul(aa,x2))) > 4.0e-11_dp*norm) &
            error stop 'real convolution matrix product'
        bc(1,:) = c
        bc(2,:) = c*2.0_dp
        call convolution_matrix(bc,batched,s,mode)
        call success(s,'real batch convolution')
        if (any(batched(1,:,:) /= aa) .or. &
            any(batched(2,:,:) /= aa*2.0_dp)) error stop 'real batch constructor'
        if (m <= 5 .and. n <= 5) then
            allocate(ic(m),imatrix(k,n))
            ic = [(i, i=1,m)]
            call convolution_matrix(ic,imatrix,s,mode)
            call success(s,'integer convolution')
            if (any(imatrix(:,1) /= int_projection(ic,k,n,mode))) &
                error stop 'integer convolution entries'
        end if
    end subroutine case_real

    function int_projection(ic,k,n,mode) result(c0)
        integer, intent(in) :: ic(:) !! integer coefficients
        integer, intent(in) :: k !! output rows
        integer, intent(in) :: n !! input length
        character(*), intent(in) :: mode !! mode
        integer :: c0(k), shift, i, loc
        shift = 0
        if (mode=='same') shift=(min(size(ic),n)-1)/2
        if (mode=='valid') shift=min(size(ic),n)-1
        c0 = 0
        do i=1,k
            loc=i+shift
            if (loc<=size(ic)) c0(i)=ic(loc)
        end do
    end function int_projection

    subroutine case_complex(m,n,mode)
        integer, intent(in) :: m !! coefficient count
        integer, intent(in) :: n !! input length
        character(*), intent(in) :: mode !! convolution mode
        complex(dp) :: c(m), x(n), x2(n,2)
        real(dp) :: cr(m), xr(n), xr2(n,2)
        complex(dp), allocatable :: aa(:,:), yy(:), y2(:,:), batched(:,:,:)
        complex(dp) :: bc(2,m)
        integer :: k, i, j, s
        real(dp) :: norm

        k=convolution_matrix_rows(m,n,mode,s)
        call success(s,'complex geometry')
        allocate(aa(k,n),yy(k),y2(k,2),batched(2,k,n))
        do i=1,m
            c(i) = cmplx(cos(real(i*13+m,dp)),sin(real(i*5+2*m,dp)),dp)
        end do
        do i=1,n
            x(i) = cmplx(cos(real(i*3+n,dp)),sin(real(i*7+n,dp)),dp)
        end do
        cr = real(c,dp)
        xr = real(x,dp)
        do j=1,2
            x2(:,j)=x*real(j,dp)+cmplx(0.5_dp,-0.2_dp,dp)
            xr2(:,j)=xr*real(j,dp)-0.1_dp
        end do
        call convolution_matrix(c,aa,s,mode)
        call success(s,'complex constructor')
        call matmul_convolution_matrix(c,x,yy,s,mode)
        call success(s,'complex convolution vector')
        norm = max(1.0_dp,maxval(abs(matmul(aa,x))))
        if (maxval(abs(yy-matmul(aa,x))) > 4.0e-11_dp*norm) &
            error stop 'complex vector product'
        call matmul_convolution_matrix(c,x2,y2,s,mode)
        call success(s,'complex convolution matrix')
        norm = max(1.0_dp,maxval(abs(matmul(aa,x2))))
        if (maxval(abs(y2-matmul(aa,x2))) > 4.0e-11_dp*norm) &
            error stop 'complex matrix product'
        call matmul_convolution_matrix(cr,x,yy,s,mode)
        call success(s,'real coefficients complex input')
        norm = max(1.0_dp,maxval(abs(matmul(cmplx(real(aa,dp),kind=dp),x))))
        if (maxval(abs(yy-matmul(cmplx(real(aa,dp),kind=dp),x))) > 4.0e-11_dp*norm) &
            error stop 'real coefficients complex vector'
        call matmul_convolution_matrix(c,xr,yy,s,mode)
        call success(s,'complex coefficients real input')
        norm = max(1.0_dp,maxval(abs(matmul(aa,cmplx(xr,kind=dp)))))
        if (maxval(abs(yy-matmul(aa,cmplx(xr,kind=dp)))) > 4.0e-11_dp*norm) &
            error stop 'complex coefficients real vector'
        call matmul_convolution_matrix(cr,x2,y2,s,mode)
        call success(s,'real coefficients complex RHS matrix')
        norm = max(1.0_dp,maxval(abs(matmul(cmplx(real(aa,dp),kind=dp),x2))))
        if (maxval(abs(y2-matmul(cmplx(real(aa,dp),kind=dp),x2))) > 4.0e-11_dp*norm) &
            error stop 'real coefficients complex RHS matrix'
        call matmul_convolution_matrix(c,xr2,y2,s,mode)
        call success(s,'complex coefficients real RHS matrix')
        norm = max(1.0_dp,maxval(abs(matmul(aa,cmplx(xr2,kind=dp)))))
        if (maxval(abs(y2-matmul(aa,cmplx(xr2,kind=dp)))) > 4.0e-11_dp*norm) &
            error stop 'complex coefficients real RHS matrix'
        bc(1,:) = c
        bc(2,:) = 2.0_dp*c
        call convolution_matrix(bc,batched,s,mode)
        call success(s,'complex batch constructor')
        if (any(batched(1,:,:) /= aa) .or. &
            any(batched(2,:,:) /= aa*2.0_dp)) error stop 'complex batch constructor mismatch'
    end subroutine case_complex

    subroutine case_scale(factor)
        real(dp), intent(in) :: factor !! uniform coefficient scale
        real(dp) :: c(3), x(4), y(6), mat(6,4)
        integer :: s
        c = [1.0_dp,-2.0_dp,3.0_dp]*factor
        x = [3.0_dp,0.5_dp,-1.0_dp,0.25_dp]
        call convolution_matrix(c,mat,s)
        call success(s,'scaled constructor')
        call matmul_convolution_matrix(c,x,y,s)
        call success(s,'scaled product')
        if (maxval(abs(y-matmul(mat,x)))/factor > 3.0e-11_dp) &
            error stop 'scaled product mismatch'
    end subroutine case_scale
end program test_linalg_convolution_matrix
