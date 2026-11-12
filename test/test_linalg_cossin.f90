! SPDX-License-Identifier: MIT
program test_linalg_cossin
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_cossin, only : cossin, cossin_separate
    implicit none
    integer :: m, p, q, k, stat, ncases, i, kind_case, flip_case
    complex(dp), allocatable :: xc(:, :), uc(:, :), vhc(:, :), ac(:, :)
    complex(dp), allocatable :: u1c(:, :), u2c(:, :), v1c(:, :), v2c(:, :)
    real(dp), allocatable :: xr(:, :), ur(:, :), vtr(:, :), cs(:, :), theta(:)
    real(dp) :: error, max_error
    logical :: swapped

    ncases = 0
    max_error = 0.0_dp
    do kind_case = 1, 2
        do m = 2, 13
            allocate(xc(m,m),uc(m,m),vhc(m,m),ac(m,m),xr(m,m),ur(m,m),vtr(m,m),cs(m,m))
            do k = 1, 3
                call generate_unitary(xc,kind_case==2,k)
                if (k == 2) then
                    ! Multiplying by a second unitary randomizes the phase gauge.
                    call generate_unitary(ac,kind_case==2,k+15)
                    xc = matmul(xc,ac)
                end if
                if (k == 3) then
                    xc = cmplx(0.0_dp,0.0_dp,dp)
                    do i = 1, m
                        xc(i,i) = cmplx(1.0_dp,0.0_dp,dp)
                    end do
                end if
                xr = real(xc,dp)
                do p = 1, m-1
                    do q = 1, m-1
                        do flip_case = 0, 1
                            swapped = flip_case == 1
                            if (kind_case == 2) then
                                call cossin(xc,p,q,uc,cs,vhc,stat,swap_sign=swapped)
                                call check(stat == linalg_status_success,'complex factorization status')
                                error = maxval(abs(matmul(matmul(uc,cmplx(cs,0.0_dp,dp)),vhc)-xc))
                                call check(error < 3.0e-10_dp,'complex reconstruction')
                                call check(unitary_error(uc) < 3.0e-10_dp,'left complex unitary')
                                call check(unitary_error(vhc) < 3.0e-10_dp,'right complex unitary')
                            else
                                call cossin(xr,p,q,ur,cs,vtr,stat,swap_sign=swapped)
                                call check(stat == linalg_status_success,'real factorization status')
                                error = maxval(abs(matmul(matmul(ur,cs),vtr)-xr))
                                call check(error < 3.0e-10_dp,'real reconstruction')
                            end if
                            max_error = max(max_error,error)
                            ncases = ncases+1
                        end do
                    end do
                end do
            end do
            deallocate(xc,uc,vhc,ac,xr,ur,vtr,cs)
        end do
    end do

    m = 7; p = 3; q = 5
    allocate(xc(m,m),uc(m,m),vhc(m,m),cs(m,m))
    call generate_unitary(xc,.true.,144)
    allocate(u1c(p,p),u2c(m-p,m-p),v1c(q,q),v2c(m-q,m-q),theta(min(p,m-p,q,m-q)))
    call cossin_separate(xc,p,q,u1c,u2c,theta,v1c,v2c,stat)
    call check(stat == linalg_status_success,'separate status')
    call check(all(theta >= 0.0_dp) .and. all(theta <= acos(-1.0_dp)/2), 'angles in [0,pi/2]')
    call cossin(xc,p,q,uc,cs,vhc,stat,compute_u=.false.,compute_vh=.false.)
    call check(stat == linalg_status_success .and. all(uc == (0.0_dp,0.0_dp)), 'disabled U')
    call check(all(vhc == (0.0_dp,0.0_dp)), 'disabled Vh')
    call cossin(xc,p,q,uc(:0,:0),cs,vhc(:0,:0),stat,compute_u=.false.,compute_vh=.false.)
    call check(stat == linalg_status_success, 'zero-size disabled factors')
    call cossin(xc,0,q,uc,cs,vhc,stat)
    call check(stat == linalg_status_invalid_shape, 'reject p=0')
    call cossin(xc,p,m,uc,cs,vhc,stat)
    call check(stat == linalg_status_invalid_shape, 'reject q=m')
    call cossin(xc,p,q,uc(:m-1,:m-1),cs,vhc,stat)
    call check(stat == linalg_status_invalid_shape, 'reject undersized U')
    call cossin_separate(xc,p,q,u1c(:p-1,:p-1),u2c,theta,v1c,v2c,stat)
    call check(stat == linalg_status_invalid_shape, 'reject undersized U1')
    xc(1,1) = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
    call cossin(xc,p,q,uc,cs,vhc,stat)
    call check(stat == linalg_status_no_convergence, 'reject nonfinite X')
    xc(1,1) = 9.0_dp
    call cossin(xc,p,q,uc,cs,vhc,stat)
    call check(stat == linalg_status_no_convergence, 'reject nonunitary X')
    print '(a,i0,a,es11.3)', 'PASS cossin native cases=',ncases,' max reconstruction error=',max_error

contains

    subroutine check(ok, message)
        logical, intent(in) :: ok !! test condition
        character(len=*), intent(in) :: message !! failure description
        if (.not. ok) then
            print '(2a)', 'FAIL: ',message
            error stop 1
        end if
    end subroutine check

    real(dp) function unitary_error(a) result(err)
        complex(dp), intent(in) :: a(:, :) !! square matrix
        complex(dp) :: prod(size(a,1),size(a,1))
        integer :: j
        prod = matmul(transpose(conjg(a)),a)
        do j=1,size(a,1)
            prod(j,j) = prod(j,j)-1.0_dp
        end do
        err = maxval(abs(prod))
    end function unitary_error

    subroutine generate_unitary(a, make_complex, seed)
        complex(dp), intent(out) :: a(:, :) !! square unitary matrix
        logical, intent(in) :: make_complex !! include complex coefficients
        integer, intent(in) :: seed !! deterministic variant index
        complex(dp) :: z(size(a,1)), val
        real(dp) :: norm
        integer :: i, j, pass, n
        n=size(a,1)
        a = cmplx(0.0_dp,0.0_dp,dp)
        do j=1,n
            do i=1,n
                val=cmplx(sin(real(i*i*97+j*j*53+i*j*19+seed*31,dp)),0.0_dp,dp)
                if (make_complex) val=val+cmplx(0.0_dp,cos(real(i*i*71+j*j*47+i*j*29+seed*19,dp)),dp)
                z(i)=val
            end do
            do pass=1,2
                do i=1,j-1
                    z=z-dot_product(a(:,i),z)*a(:,i)
                end do
            end do
            norm=sqrt(sum(abs(z)**2))
            call check(norm > 1.0e-10_dp,'unitary test matrix creation')
            a(:,j)=z/norm
        end do
    end subroutine generate_unitary
end program test_linalg_cossin
