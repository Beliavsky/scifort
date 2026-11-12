! SPDX-License-Identifier: MIT
! Check SciFort B-spline to PPoly conversion against SciPy.
program ppoly_from_spline_driver
    use scifort_kinds, only : dp
    use scifort_interpolate, only : bspline, bspline_array_interpolator_type, &
        ppoly_from_spline, ppoly_type
    implicit none
    type(bspline_array_interpolator_type) :: spline
    type(ppoly_type) :: poly
    real(dp), allocatable :: t(:), c(:,:), pc(:,:,:), q(:), vals(:,:)
    integer :: cases, i, k, nc, nch, nq, j, code

    read(*,*) cases
    do i = 1,cases
        read(*,*) k,nc,nch,nq
        allocate(t(nc+k+1),c(nc,nch),q(nq))
        read(*,*) t
        read(*,*) c
        read(*,*) q
        spline = bspline(t,c,k,status=code)
        if (code/=0) error stop 'from_spline driver: invalid spline'
        poly = ppoly_from_spline(spline,status=code)
        if (code/=0 .or. .not. poly%valid()) error stop 'from_spline driver: conversion failed'
        pc = poly%coefficients()
        vals = poly%evaluate(q)
        write(*,'(a,1x,i0)') 'CASE',i
        write(*,'(a,1x,i0,1x,i0,1x,i0)') 'COEFF',size(pc,1),size(pc,2),size(pc,3)
        write(*,'(*(1x,es26.17e3))') pc
        write(*,'(a,1x,i0,1x,i0)') 'VALUES',size(vals,1),size(vals,2)
        write(*,'(*(1x,es26.17e3))') vals
        deallocate(t,c,q)
    end do
end program ppoly_from_spline_driver
