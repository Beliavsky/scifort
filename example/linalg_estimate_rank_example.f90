! SPDX-License-Identifier: MIT
! Build from the project root with gfortran and the prerequisite modules.
program linalg_estimate_rank_example
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_estimate_rank, only : estimate_rank, estimate_rank_operator, &
        estimate_rank_operator_complex
    implicit none
    real(dp) :: a(6,5)
    complex(dp) :: c(6,5)
    integer :: i, j, rank, status

    a = 0.0_dp
    a(:,1) = [1.0_dp,2.0_dp,0.0_dp,3.0_dp,1.0_dp,-1.0_dp]
    a(:,2) = [0.0_dp,1.0_dp,1.0_dp,-1.0_dp,2.0_dp,3.0_dp]
    a(:,3) = a(:,1) + a(:,2)
    a(:,4) = 2.0_dp*a(:,1) - a(:,2)
    a(:,5) = 0.5_dp*a(:,2)
    do j = 1, 5
        do i = 1, 6
            c(i,j) = a(i,j)*cmplx(cos(0.12_dp*real(i+j,dp)), &
                sin(0.12_dp*real(i+j,dp)),dp)
        end do
    end do

    call estimate_rank(a,1.0e-9_dp,rank,status,seed=11_int64)
    if (status /= linalg_status_success) error stop 'real dense estimator failed'
    print '(a,i0)', 'Real dense rank: ',rank
    call estimate_rank_operator(6,5,mv_real,1.0e-9_dp,rank,status,seed=11_int64)
    if (status /= linalg_status_success) error stop 'real operator estimator failed'
    print '(a,i0)', 'Real operator rank: ',rank
    call estimate_rank(c,1.0e-9_dp,rank,status,seed=11_int64)
    if (status /= linalg_status_success) error stop 'complex dense estimator failed'
    print '(a,i0)', 'Complex dense rank: ',rank
    call estimate_rank_operator_complex(6,5,mv_complex,1.0e-9_dp,rank,status,seed=11_int64)
    if (status /= linalg_status_success) error stop 'complex operator estimator failed'
    print '(a,i0)', 'Complex operator rank: ',rank

contains
    subroutine mv_real(x,y,istat)
        real(dp), intent(in) :: x(:) !! right-hand vector
        real(dp), intent(out) :: y(:) !! product A*x
        integer, intent(out) :: istat !! zero on success
        y = matmul(a,x)
        istat = 0
    end subroutine mv_real
    subroutine mv_complex(x,y,istat)
        complex(dp), intent(in) :: x(:) !! right-hand vector
        complex(dp), intent(out) :: y(:) !! product C*x
        integer, intent(out) :: istat !! zero on success
        y = matmul(c,x)
        istat = 0
    end subroutine mv_complex
end program linalg_estimate_rank_example
