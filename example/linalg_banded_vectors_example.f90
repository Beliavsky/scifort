! SPDX-License-Identifier: MIT
program linalg_banded_vectors_example
    use scifort_kinds, only : dp
    use scifort_linalg_banded_eigen, only : eig_banded
    implicit none
    integer, parameter :: n=2000, nkeep=3
    real(dp) :: ab(2,n), w(nkeep), vectors(n,nkeep)
    real(dp) :: residual, current
    integer :: status, selected, i, j

    ab=0.0_dp
    ab(1,:)=2.0_dp
    ab(2,1:n-1)=-1.0_dp
    call eig_banded(ab,w,status,vectors,lower=.true., &
        subset_by_index=[n/2,n/2+nkeep-1],num_selected=selected,method='inverse')
    if (status /= 0 .or. selected /= nkeep) error stop 'banded eigenvectors failed'
    print '(a,i0,a)', 'Computed ',selected,' selected eigenpairs without a dense input matrix:'
    do j=1,nkeep
        residual=0.0_dp
        do i=1,n
            current=(2.0_dp-w(j))*vectors(i,j)
            if (i>1) current=current-vectors(i-1,j)
            if (i<n) current=current-vectors(i+1,j)
            residual=max(residual,abs(current))
        end do
        print '(i0,2es22.12)', j,w(j),residual
    end do
end program linalg_banded_vectors_example
