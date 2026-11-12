! SPDX-License-Identifier: MIT
! Benchmarks SciFort's default native root isolation against pinned rpoly.
! Requires the real upstream package: no benchmarks are fabricated.
program benchmark_polyroots
    use scifort_kinds, only : dp
    use scifort_interpolate, only : ppoly, ppoly_type
    use scifort_optional_polyroots, only : polyroots_segment_solver
    implicit none
    type(ppoly_type) :: pp
    real(dp), allocatable :: coefficients(:), native(:), upstream(:)
    real(dp) :: start, native_seconds, upstream_seconds
    integer :: degree, i
    print '(a)', 'degree, repetitions, native_seconds, upstream_seconds, max_root_difference'
    do degree = 5, 17, 4
        call make_polynomial(degree,coefficients)
        pp=ppoly(reshape(coefficients,[degree+1,1]),[0.0_dp,2.0_dp], &
            extrapolate=.false.)
        native=pp%roots()
        upstream=pp%roots(root_solver=polyroots_segment_solver)
        if (size(native)/=size(upstream)) then
            print '(a,i0)', 'root-count mismatch for polynomial degree ',degree
            cycle
        end if
        call cpu_time(start)
        do i=1,100
            native=pp%roots()
        end do
        call cpu_time(native_seconds)
        native_seconds=native_seconds-start
        call cpu_time(start)
        do i=1,100
            upstream=pp%roots(root_solver=polyroots_segment_solver)
        end do
        call cpu_time(upstream_seconds)
        upstream_seconds=upstream_seconds-start
        print '(i0,a,i0,a,es12.4,a,es12.4,a,es12.4)', &
            degree,',',100,',',native_seconds,',',upstream_seconds,',', &
            maxval(abs(native-upstream))
    end do
contains
    subroutine make_polynomial(degree,coeff)
        integer, intent(in) :: degree !! odd polynomial degree (5,9,13,17)
        real(dp), allocatable, intent(out) :: coeff(:) !! descending power coefficients
        integer :: j
        coeff=[1.0_dp]
        call multiply(coeff,[-0.25_dp,1.0_dp])
        call multiply(coeff,[-0.8_dp,1.0_dp])
        call multiply(coeff,[-1.6_dp,1.0_dp])
        do j=1,(degree-3)/2
            call multiply(coeff,[0.3_dp+real(j,dp)*0.1_dp,0.0_dp,1.0_dp])
        end do
        coeff=coeff(size(coeff):1:-1)
    end subroutine make_polynomial

    subroutine multiply(coeff,factor)
        real(dp), allocatable, intent(inout) :: coeff(:) !! ascending polynomial coefficients
        real(dp), intent(in) :: factor(:) !! ascending coefficients of factor
        real(dp), allocatable :: tmp(:)
        integer :: a,b
        allocate(tmp(size(coeff)+size(factor)-1))
        tmp=0.0_dp
        do a=1,size(coeff)
            do b=1,size(factor)
                tmp(a+b-1)=tmp(a+b-1)+coeff(a)*factor(b)
            end do
        end do
        call move_alloc(tmp,coeff)
    end subroutine multiply
end program benchmark_polyroots
