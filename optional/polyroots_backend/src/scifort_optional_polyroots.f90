! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! SciFort-owned adapter to Jacob Williams's polyroots-fortran v1.2.0.
! No third-party code is copied here. Compile this source ONLY when the
! optional polyroots-fortran dependency is present (see sibling fpm.toml).
! The adapter calls the standalone real Jenkins-Traub rpoly routine. Unlike
! other polyroots-fortran methods, rpoly does not use BLAS/LAPACK internally;
! nevertheless the upstream fpm.toml currently links both libraries.
module scifort_optional_polyroots
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_interpolate_piecewise, only : piecewise_native_real_roots
    use polyroots_module, only : rpoly, polyroots_module_rk
    implicit none
    private
    public :: polyroots_segment_solver
contains
    function polyroots_segment_solver(coeff,lower,upper) result(roots)
        real(dp), intent(in) :: coeff(:) !! descending local power coefficients
        real(dp), intent(in) :: lower !! inclusive lower local interval bound
        real(dp), intent(in) :: upper !! inclusive upper local interval bound
        real(dp), allocatable :: roots(:)
        real(dp), allocatable :: zr(:), zi(:), candidates(:)
        real(dp) :: root, residual, bound, imagtol, threshold, temp
        integer :: first, degree, istat, i, j, n
        integer, parameter :: min_upstream_degree = 4
        allocate(roots(0))
        if (dp /= polyroots_module_rk) error stop 'polyroots adapter requires REAL64'
        if (size(coeff)<2 .or. lower>upper) return
        if (.not. all(ieee_is_finite(coeff))) return
        first=1
        do while (first<size(coeff))
            if (coeff(first)/=0.0_dp) exit
            first=first+1
        end do
        degree=size(coeff)-first
        if (degree < min_upstream_degree) then
            roots=piecewise_native_real_roots(coeff,lower,upper)
            return
        end if
        allocate(zr(degree),zi(degree),candidates(degree))
        call rpoly(coeff(first:),degree,zr,zi,istat)
        if (istat/=0) then
            roots=piecewise_native_real_roots(coeff,lower,upper)
            return
        end if
        n=0
        do i=1,degree
            if (.not. ieee_is_finite(zr(i)) .or. .not. ieee_is_finite(zi(i))) cycle
            root=zr(i)
            imagtol=3.0e-9_dp*max(1.0_dp,abs(root))
            if (abs(zi(i))>imagtol) cycle
            threshold=64.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(lower),abs(upper))
            if (root<lower-threshold .or. root>upper+threshold) cycle
            root=min(upper,max(lower,root))
            ! Validate the candidate against the real polynomial. A small
            ! imaginary part alone is not sufficient evidence of a real root.
            residual=coeff(first)
            bound=abs(residual)
            do j=first+1,size(coeff)
                residual=residual*root+coeff(j)
                bound=bound*abs(root)+abs(coeff(j))
            end do
            if (.not. ieee_is_finite(residual) .or. .not. ieee_is_finite(bound)) cycle
            if (abs(residual)>1.0e-9_dp*bound) cycle
            if (n>0) then
                if (any(abs(candidates(:n)-root)<1.0e-8_dp*max(1.0_dp,abs(root)))) cycle
            end if
            n=n+1
            candidates(n)=root
        end do
        roots=candidates(:n)
        ! rpoly does not promise sorted output. PPoly's segment contract does.
        do i=2,n
            temp=roots(i)
            j=i-1
            do while (j>=1)
                if (roots(j)<=temp) exit
                roots(j+1)=roots(j)
                j=j-1
            end do
            roots(j+1)=temp
        end do
    end function polyroots_segment_solver
end module scifort_optional_polyroots
