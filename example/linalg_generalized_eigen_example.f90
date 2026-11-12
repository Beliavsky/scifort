! SPDX-License-Identifier: MIT
! Homogeneous generalized eigenvalues, including an infinite root.
program linalg_generalized_eigen_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_generalized_eigen, only : eig, eigvals
    implicit none
    real(dp) :: a(3,3), b(3,3)
    complex(dp) :: w(3), alpha(3), beta(3), vr(3,3), vl(3,3)
    real(dp) :: residual, worst
    integer :: status, j

    a = 0.0_dp
    b = 0.0_dp
    a(1,1) = 3.0_dp
    a(2,2) = 5.0_dp
    a(3,3) = 7.0_dp
    b(1,1) = 2.0_dp
    b(2,2) = 0.0_dp ! one infinite eigenvalue
    b(3,3) = 4.0_dp

    call eig(a,b,w,status,vr=vr,vl=vl,alpha=alpha,beta=beta)
    if (status /= linalg_status_success) error stop 'generalized eig failed'
    worst = 0.0_dp
    do j=1,3
        residual = maxval(abs(beta(j)*matmul(cmplx(a,kind=dp),vr(:,j)) - &
            alpha(j)*matmul(cmplx(b,kind=dp),vr(:,j))))
        worst = max(worst,residual)
        if (abs(beta(j)) < tiny(1.0_dp)) then
            print '(a,i0,a)', 'Eigenvalue ',j, ': infinity'
        else
            print '(a,i0,a,2es14.5)', 'Eigenvalue ',j, ': ',real(w(j),dp),aimag(w(j))
        end if
    end do
    print '(a,es12.3)', 'Largest right-pencil residual: ',worst
    if (worst > 1.0e-11_dp) error stop 'generalized eig residual failed'
    call eigvals(a,b,w,status,alpha=alpha,beta=beta)
    if (status /= linalg_status_success) error stop 'generalized eigvals failed'
end program linalg_generalized_eigen_example
