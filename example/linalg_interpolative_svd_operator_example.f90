! SPDX-License-Identifier: MIT
! Demonstrate an SVD of a 120-by-180 operator without constructing the matrix.
program linalg_interpolative_svd_operator_example
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_interpolative_svd_operator, only : randomized_interpolative_svd_operator
    implicit none
    integer, parameter :: m = 120, n = 180, r = 3
    real(dp) :: left(m,r), right(n,r), x(n), y(m), approx(m)
    real(dp) :: u(m,r), s(r), v(n,r), err, estimated
    real(dp), allocatable :: ut(:, :), st(:), vt(:, :)
    complex(dp) :: leftc(m,r), rightc(n,r), xc(n), yc(m), approxc(m)
    complex(dp) :: uc(m,r), vc(n,r)
    complex(dp), allocatable :: utc(:, :), vtc(:, :)
    integer :: i, j, k, status

    do j = 1, r
        do i = 1, m
            left(i,j) = sin(real(i*j,dp)*0.13_dp) + 0.02_dp*real(j,dp)
            leftc(i,j) = cmplx(left(i,j),cos(real(i+2*j,dp)*0.17_dp),dp)
        end do
        do i = 1, n
            right(i,j) = cos(real(i*j,dp)*0.19_dp) + 0.03_dp*real(j,dp)
            rightc(i,j) = cmplx(right(i,j),sin(real(2*i+j,dp)*0.09_dp),dp)
        end do
    end do
    x = [(sin(real(i,dp)*0.25_dp),i=1,n)]
    xc = cmplx(x,[(cos(real(i,dp)*0.18_dp),i=1,n)],dp)

    call randomized_interpolative_svd_operator(m,n,r,mv,rmv,u,s,v,status,seed=773_int64)
    if (status /= linalg_status_success) error stop 'real fixed-rank operator SVD failed'
    call mv(x,y,status)
    approx = matmul(u,s*matmul(transpose(v),x))
    err = sqrt(sum((y-approx)**2))/sqrt(sum(y**2))
    print '(a,es12.4)', 'real fixed-rank relative application residual: ',err
    if (err > 1.0e-9_dp) error stop 'real SVD application residual'

    call randomized_interpolative_svd_operator(m,n,1.0e-6_dp,mv,rmv,k, &
        ut,st,vt,status,seed=773_int64,estimated_error=estimated)
    if (status /= linalg_status_success .or. k > r) error stop 'real adaptive operator SVD failed'
    print '(a,i0,a,es12.4)', 'real adaptive rank: ',k, '; estimated residual: ',estimated

    call randomized_interpolative_svd_operator(m,n,r,mvc,rmvc,uc,s,vc,status,seed=773_int64)
    if (status /= linalg_status_success) error stop 'complex fixed-rank operator SVD failed'
    call mvc(xc,yc,status)
    approxc = matmul(uc,s*matmul(transpose(conjg(vc)),xc))
    err = sqrt(sum(abs(yc-approxc)**2))/sqrt(sum(abs(yc)**2))
    print '(a,es12.4)', 'complex fixed-rank relative application residual: ',err
    if (err > 1.0e-9_dp) error stop 'complex SVD application residual'

    call randomized_interpolative_svd_operator(m,n,1.0e-6_dp,mvc,rmvc,k, &
        utc,st,vtc,status,seed=773_int64,estimated_error=estimated)
    if (status /= linalg_status_success .or. k > r) error stop 'complex adaptive operator SVD failed'
    print '(a,i0,a,es12.4)', 'complex adaptive rank: ',k, '; estimated residual: ',estimated
contains
    subroutine mv(vec,result,istat)
        real(dp), intent(in) :: vec(:) !! n-vector
        real(dp), intent(out) :: result(:) !! m-vector A*vec
        integer, intent(out) :: istat !! callback status
        result = matmul(left,matmul(transpose(right),vec))
        istat = 0
    end subroutine mv
    subroutine rmv(vec,result,istat)
        real(dp), intent(in) :: vec(:) !! m-vector
        real(dp), intent(out) :: result(:) !! n-vector A^T*vec
        integer, intent(out) :: istat !! callback status
        result = matmul(right,matmul(transpose(left),vec))
        istat = 0
    end subroutine rmv
    subroutine mvc(vec,result,istat)
        complex(dp), intent(in) :: vec(:) !! n-vector
        complex(dp), intent(out) :: result(:) !! m-vector A*vec
        integer, intent(out) :: istat !! callback status
        result = matmul(leftc,matmul(transpose(conjg(rightc)),vec))
        istat = 0
    end subroutine mvc
    subroutine rmvc(vec,result,istat)
        complex(dp), intent(in) :: vec(:) !! m-vector
        complex(dp), intent(out) :: result(:) !! n-vector A^H*vec
        integer, intent(out) :: istat !! callback status
        result = matmul(rightc,matmul(transpose(conjg(leftc)),vec))
        istat = 0
    end subroutine rmvc
end program linalg_interpolative_svd_operator_example
