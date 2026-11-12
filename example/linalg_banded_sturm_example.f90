! SPDX-License-Identifier: MIT
program linalg_banded_sturm_example
    use scifort_kinds, only : dp
    use scifort_linalg_banded_eigen, only : eigvals_banded, eig_banded
    use scifort_linalg, only : linalg_status_success
    implicit none

    integer, parameter :: n=8, kd=2
    real(dp) :: ab(kd+1,n), w(3), allw(n)
    complex(dp) :: ac(kd+1,n)
    integer :: i, status, found

    ab = 0.0_dp
    ab(1,:) = 4.0_dp
    ab(2,1:n-1) = -1.0_dp
    ab(3,1:n-2) = -0.25_dp

    ! Only three eigenvalues are requested. Banded Sturm bisection avoids
    ! the n-by-n dense matrix used by the full eigenvector method.
    call eigvals_banded(ab,w,status,lower=.true., &
        subset_by_index=[3,5],num_selected=found)
    if (status /= linalg_status_success .or. found /= 3) error stop 'real selection failed'
    print '(A)', 'Selected real symmetric band eigenvalues:'
    print '(3F16.10)', w

    call eig_banded(ab,allw,status,lower=.true.)
    if (status /= linalg_status_success) error stop 'dense reference failed'
    if (maxval(abs(allw(3:5)-w)) > 1.0e-10_dp) error stop 'selection mismatch'

    ac = cmplx(ab,0.0_dp,dp)
    do i=1,n-1
        ac(2,i) = cmplx(0.0_dp,-1.0_dp,dp)
    end do
    call eigvals_banded(ac,w,status,lower=.true., &
        subset_by_index=[3,5],num_selected=found)
    if (status /= linalg_status_success .or. found /= 3) error stop 'complex selection failed'
    print '(A)', 'Selected complex Hermitian band eigenvalues:'
    print '(3F16.10)', w
end program linalg_banded_sturm_example
