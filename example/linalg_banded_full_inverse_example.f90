! SPDX-License-Identifier: MIT
program linalg_banded_full_inverse_example
    use scifort_kinds, only : dp
    use scifort_linalg_banded_eigen, only : eig_banded
    use scifort_linalg, only : linalg_status_success
    implicit none

    integer, parameter :: n=24
    real(dp) :: ab(3,n), eigenvalues(n), eigenvectors(n,n)
    integer :: j, status, found

    ! A real symmetric pentadiagonal matrix in SciPy lower-packed format.
    ab=0.0_dp
    ab(1,:)=4.0_dp
    ab(2,1:n-1)=-0.5_dp
    ab(3,1:n-2)=-0.15_dp
    do j=1,n
        ab(1,j)=ab(1,j)+0.001_dp*real(j,dp)
    end do

    ! Full eigenspectrum without a dense *input* matrix expansion.
    ! A full set of eigenvectors must still contain n*n elements.
    call eig_banded(ab,eigenvalues,status,eigenvectors,lower=.true., &
        method='inverse',num_selected=found)
    if (status /= linalg_status_success .or. found /= n) error stop 1
    print '(a,i0)', 'Computed eigenpairs: ',found
    print '(a,2(1x,es18.10))', 'Extreme eigenvalues:',eigenvalues(1),eigenvalues(n)
end program linalg_banded_full_inverse_example
