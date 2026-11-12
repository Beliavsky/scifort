! SPDX-License-Identifier: MIT
! A rank-three matrix-free operator: no dense m-by-n matrix is created.
program linalg_interpolative_operator_example
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_interpolative_operator, only : randomized_interp_decomp_operator, &
        randomized_interp_decomp_operator_complex
    implicit none
    integer, parameter :: m = 120, n = 180, r = 3
    real(dp) :: u(m,r), v(n,r), x(n), y(m), y_id(m), projected(r), err, eps
    complex(dp) :: uc(m,r), vc(n,r), xc(n), yc(m), yc_id(m)
    real(dp), allocatable :: proj(:, :), b(:, :)
    complex(dp), allocatable :: projc(:, :), bc(:, :)
    integer, allocatable :: idx(:), idxc(:)
    integer :: i, j, k, status
    do j = 1, r
        do i = 1, m
            u(i,j) = sin(real(i*j,dp)*0.13_dp) + real(j,dp)*0.02_dp
            uc(i,j) = cmplx(u(i,j),cos(real(i+2*j,dp)*0.17_dp),dp)
        end do
        do i = 1, n
            v(i,j) = cos(real(i*j,dp)*0.19_dp) + real(j,dp)*0.03_dp
            vc(i,j) = cmplx(v(i,j),sin(real(2*i+j,dp)*0.09_dp),dp)
        end do
    end do
    x = [(cos(real(i,dp)*0.31_dp),i=1,n)]
    xc = cmplx(x,[(sin(real(i,dp)*0.25_dp),i=1,n)],dp)

    call randomized_interp_decomp_operator(m,n,r,mv,rmv,idx,proj,status, &
        seed=77631_int64,oversampling=7,skeleton=b)
    if (status /= linalg_status_success) error stop 'fixed-rank real ID failed'
    projected = x(idx(1:r)) + matmul(proj,x(idx(r+1:n)))
    y_id = matmul(b,projected)
    call mv(x,y,status)
    err = sqrt(sum((y-y_id)**2))/sqrt(sum(y**2))
    print '(a,i0,a,es12.4)', 'real fixed rank = ',r,'; relative probe error = ',err
    if (err > 1.0e-9_dp) error stop 'real ID reconstruction error'

    eps = 1.0e-6_dp
    call randomized_interp_decomp_operator(m,n,eps,mv,rmv,k,idx,proj,status, &
        seed=77631_int64,skeleton=b,estimated_error=err)
    if (status /= linalg_status_success) error stop 'adaptive real ID failed'
    print '(a,i0,a,es12.4)', 'real adaptive rank = ',k,'; estimated error = ',err
    if (err > eps .or. k > r) error stop 'real adaptive error or rank'

    call randomized_interp_decomp_operator_complex(m,n,r,mvc,rmvc,idxc,projc,status, &
        seed=77631_int64,oversampling=7,skeleton=bc)
    if (status /= linalg_status_success) error stop 'fixed-rank complex ID failed'
    yc_id = matmul(bc,xc(idxc(1:r))+matmul(projc,xc(idxc(r+1:n))))
    call mvc(xc,yc,status)
    err = sqrt(sum(abs(yc-yc_id)**2))/sqrt(sum(abs(yc)**2))
    print '(a,i0,a,es12.4)', 'complex fixed rank = ',r,'; relative probe error = ',err
    if (err > 1.0e-9_dp) error stop 'complex ID reconstruction error'

    call randomized_interp_decomp_operator_complex(m,n,eps,mvc,rmvc,k,idxc,projc,status, &
        seed=77631_int64,skeleton=bc,estimated_error=err)
    if (status /= linalg_status_success) error stop 'adaptive complex ID failed'
    print '(a,i0,a,es12.4)', 'complex adaptive rank = ',k,'; estimated error = ',err
    if (err > eps .or. k > r) error stop 'complex adaptive error or rank'

contains
    subroutine mv(vec,result,istat)
        real(dp), intent(in) :: vec(:) !! n-vector
        real(dp), intent(out) :: result(:) !! m-vector A*vec
        integer, intent(out) :: istat !! 0 on success
        result = matmul(u,matmul(transpose(v),vec))
        istat = 0
    end subroutine mv
    subroutine rmv(vec,result,istat)
        real(dp), intent(in) :: vec(:) !! m-vector
        real(dp), intent(out) :: result(:) !! n-vector A^T*vec
        integer, intent(out) :: istat !! 0 on success
        result = matmul(v,matmul(transpose(u),vec))
        istat = 0
    end subroutine rmv
    subroutine mvc(vec,result,istat)
        complex(dp), intent(in) :: vec(:) !! n-vector
        complex(dp), intent(out) :: result(:) !! m-vector A*vec
        integer, intent(out) :: istat !! 0 on success
        result = matmul(uc,matmul(transpose(conjg(vc)),vec))
        istat = 0
    end subroutine mvc
    subroutine rmvc(vec,result,istat)
        complex(dp), intent(in) :: vec(:) !! m-vector
        complex(dp), intent(out) :: result(:) !! n-vector A^H*vec
        integer, intent(out) :: istat !! 0 on success
        result = matmul(vc,matmul(transpose(conjg(uc)),vec))
        istat = 0
    end subroutine rmvc
end program linalg_interpolative_operator_example
