! SPDX-License-Identifier: MIT
! Example of the public SciFort bounded simplicial global optimizer.
program optimize_shgo_example
    use scifort_kinds, only : dp
    use scifort_optimize, only : shgo, shgo_result
    implicit none
    type(shgo_result) :: result
    real(dp) :: bounds(2,2)
    integer :: j

    bounds(:,1) = -6.0_dp
    bounds(:,2) = 6.0_dp
    result = shgo(himmelblau, bounds, iters=3, maxfun=15000, &
        method='nelder-mead', local_maxiter=300)
    if (.not. result%success) error stop 'SHGO failed'
    if (result%fun > 1.0e-4_dp) error stop 'Himmelblau minimum not found'
    print '(a,2(1x,f12.7))', 'global minimizer:', result%x
    print '(a,es14.6)', 'objective: ', result%fun
    print '(a,i0,a,i0)', 'evaluations: ', result%nfev, ', completed levels: ', result%nit
    print '(a,i0)', 'distinct local minima found: ', result%nminima
    do j = 1, result%nminima
        print '(a,i0,a,es12.4,a,2f11.5)', 'minimum ', j, ': f=', &
            result%funl(j), ', x=', result%xl(:,j)
    end do
contains
    function himmelblau(x) result(f)
        real(dp), intent(in) :: x(:) !! parameter pair
        real(dp) :: f !! Himmelblau's nonconvex quartic
        f = (x(1)**2+x(2)-11.0_dp)**2 + (x(1)+x(2)**2-7.0_dp)**2
    end function himmelblau
end program optimize_shgo_example
