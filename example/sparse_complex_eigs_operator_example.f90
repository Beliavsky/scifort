! SPDX-License-Identifier: MIT
! Matrix-free shift-invert for a five-dimensional diagonal complex operator.
program sparse_complex_eigs_operator_example
    use scifort_kinds, only : dp
    use scifort_sparse_linalg_complex_eigs, only : eigs_complex_operator
    implicit none
    complex(dp) :: d(5), sigma
    complex(dp), allocatable :: lambda(:), vectors(:,:)
    real(dp), allocatable :: residuals(:)
    integer :: i, info
    d = [(cmplx(real(i,dp),0.1_dp*real(i,dp),dp),i=1,5)]
    sigma = d(3)+cmplx(0.2_dp,0.0_dp,dp)
    call eigs_complex_operator(5,1,apply,lambda,vectors,info, &
        sigma=sigma,inverse_matvec=solve_shift,ncv=5, &
        tol=1.0e-10_dp,residual_norms=residuals)
    if (info /= 0) error stop 'matrix-free eigenvalue search failed'
    print '(a,2f12.7,a,es12.3)', 'lambda = ', real(lambda(1),dp), aimag(lambda(1)), &
        ', residual = ', residuals(1)
contains
    subroutine apply(x,y,status)
        complex(dp), intent(in) :: x(:) !! Input vector.
        complex(dp), intent(out) :: y(:) !! Original diagonal operator action.
        integer, intent(out) :: status !! Zero on success.
        y = d*x
        status = 0
    end subroutine apply
    subroutine solve_shift(x,y,status)
        complex(dp), intent(in) :: x(:) !! Shifted-system RHS.
        complex(dp), intent(out) :: y(:) !! Shifted-system solution.
        integer, intent(out) :: status !! Zero on success.
        y = x/(d-sigma)
        status = 0
    end subroutine solve_shift
end program sparse_complex_eigs_operator_example
