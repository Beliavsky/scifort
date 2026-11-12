! SPDX-License-Identifier: MIT
program test_polyroots_real_stress
    use scifort_kinds, only: dp
    use scifort_optional_polyroots, only: polyroots_segment_solver
    implicit none
    real(dp), allocatable :: coeff(:), roots(:)
    integer :: degree, variant, j, ncases
    ncases = 0
    do degree = 5, 25, 4
        do variant = 1, 20
            coeff = [1.0_dp]
            call multiply(coeff,[-0.1_dp-0.001_dp*real(variant,dp), 1.0_dp])
            call multiply(coeff,[-0.7_dp-0.003_dp*real(variant,dp), 1.0_dp])
            call multiply(coeff,[-1.5_dp+0.002_dp*real(variant,dp), 1.0_dp])
            do j = 1, (degree-3)/2
                call multiply(coeff,[0.5_dp+0.02_dp*real(j,dp)+0.001_dp*real(variant,dp),0.0_dp,1.0_dp])
            end do
            roots = polyroots_segment_solver(coeff(size(coeff):1:-1), 0.0_dp, 2.0_dp)
            if (size(roots) /= 3) then
                print *, 'FAILED root count', degree, variant, size(roots)
                error stop 1
            end if
            if (maxval(abs(roots - [0.1_dp+0.001_dp*real(variant,dp), &
                    0.7_dp+0.003_dp*real(variant,dp), 1.5_dp-0.002_dp*real(variant,dp)])) > 1.0e-6_dp) then
                print *, 'FAILED roots', degree, variant, roots
                error stop 2
            end if
            ncases = ncases + 1
        end do
    end do
    print '(a,i0,a)', 'Real upstream polynomial stress test PASS (', ncases, ' cases)'
contains
    subroutine multiply(c, factor)
        real(dp), allocatable, intent(inout) :: c(:) !! Polynomial coefficients in ascending powers
        real(dp), intent(in) :: factor(:) !! Factor coefficients in ascending powers
        real(dp), allocatable :: next(:)
        integer :: a,b
        allocate(next(size(c)+size(factor)-1))
        next=0.0_dp
        do a=1,size(c)
            do b=1,size(factor)
                next(a+b-1)=next(a+b-1)+c(a)*factor(b)
            end do
        end do
        call move_alloc(next,c)
    end subroutine multiply
end program test_polyroots_real_stress
