! SPDX-License-Identifier: MIT
! Test cached cubic patches and indexed lookup on a moderately large mesh.
program test_interpolate_clough_tocher_cache
    use scifort_kinds, only: dp
    use scifort_interpolate_clough_tocher, only: clough_tocher_2d_interpolator, &
        clough_tocher_2d_interpolator_type, complex_clough_tocher_2d_interpolator_type, &
        ct_status_success, ct_status_nonconverged
    implicit none
    integer, parameter :: n=210
    real(dp) :: sites(n,2), values(n,3), p(2), expected(3)
    real(dp), allocatable :: queries(:,:), answers(:,:), one(:)
    complex(dp) :: cv(n,3)
    complex(dp), allocatable :: canswers(:,:)
    integer, allocatable :: triangles(:,:)
    type(clough_tocher_2d_interpolator_type) :: interp, stalled
    type(complex_clough_tocher_2d_interpolator_type) :: cinterp
    integer :: i,j,k,t,nquery,status

    do i=1,n
        sites(i,1)=100.0_dp+0.02_dp*modulo(sin(12.9898_dp*real(i,dp))*43758.5453_dp,1.0_dp)
        sites(i,2)=-80.0_dp+35.0_dp*modulo(sin(78.233_dp*real(i,dp))*12475.231_dp,1.0_dp)
    end do
    do i=1,n
        values(i,:)=affine(sites(i,:))
    end do
    cv=cmplx(values,2.0_dp*values,dp)
    interp=clough_tocher_2d_interpolator(sites,values,rescale=.true., &
        tol=1.0e-12_dp,maxiter=1200,fill_value=-99.0_dp,status=status)
    if (status/=ct_status_success .or. .not.interp%valid()) error stop 'cached large real construction'
    cinterp=clough_tocher_2d_interpolator(sites,cv,rescale=.true., &
        tol=1.0e-12_dp,maxiter=1200,fill_value=cmplx(-99.0_dp,5.0_dp,dp),status=status)
    if (status/=ct_status_success .or. .not.cinterp%valid()) error stop 'cached large complex construction'
    triangles=interp%simplices()
    nquery=n+4*size(triangles,2)+2
    allocate(queries(nquery,2))
    queries(:n,:)=sites
    k=n
    do t=1,size(triangles,2)
        p=0.0_dp
        do j=1,3
            p=p+sites(triangles(j,t),:)/3.0_dp
        end do
        k=k+1
        queries(k,:)=p
        do j=1,3
            k=k+1
            queries(k,:)=0.5_dp*(sites(triangles(j,t),:)+sites(triangles(mod(j,3)+1,t),:))
        end do
    end do
    queries(nquery-1,:)=[101.0_dp,0.0_dp]
    queries(nquery,:)=[90.0_dp,0.0_dp]
    answers=interp%evaluate(queries,status=status)
    if (status/=ct_status_success) error stop 'cached large real batch status'
    canswers=cinterp%evaluate(queries,status=status)
    if (status/=ct_status_success) error stop 'cached large complex batch status'
    do i=1,nquery-2
        expected=affine(queries(i,:))
        if (maxval(abs(answers(i,:)-expected))>1.0e-9_dp) then
            print *, 'bad affine', i, queries(i,:), answers(i,:), expected, maxval(abs(answers(i,:)-expected))
            error stop 'cached near edge affine value'
        end if
        if (maxval(abs(canswers(i,:)-cmplx(expected,2.0_dp*expected,dp)))>2.0e-9_dp) &
            error stop 'cached near edge complex value'
    end do
    if (any(abs(answers(nquery-1:nquery,:)+99.0_dp)>1.0e-12_dp)) error stop 'cached fill real'
    if (any(abs(canswers(nquery-1:nquery,:)-cmplx(-99.0_dp,5.0_dp,dp))>1.0e-12_dp)) &
        error stop 'cached fill complex'
    one=interp%evaluate(queries(n+1,:))
    if (maxval(abs(one-answers(n+1,:)))>1.0e-12_dp) error stop 'cached scalar/batch discrepancy'
    stalled=clough_tocher_2d_interpolator(sites,values(:,1), &
        tol=1.0e-16_dp,maxiter=1,status=status)
    if (status==ct_status_nonconverged) then
        if (.not.stalled%valid() .or. stalled%converged()) error stop 'cached nonconverged core'
        one=stalled%evaluate(sites(1,:))
        if (size(one)/=1) error stop 'cached nonconverged evaluation'
    else if (status/=ct_status_success) then
        error stop 'cached nonconverged construction status'
    end if
    print '(a,i0,a,i0,a)', 'Clough-Tocher cached/indexed batch PASS (', &
        nquery,' queries; ',size(triangles,2),' triangles)'
contains
    pure function affine(x) result(y)
        real(dp), intent(in) :: x(2) !! input physical coordinates
        real(dp) :: y(3)
        y=[2.0_dp+3.0_dp*x(1)-4.0_dp*x(2), &
           -1.0_dp+0.03_dp*x(1)+0.6_dp*x(2), &
           1.0_dp-1.0_dp*x(1)+0.05_dp*x(2)]
    end function affine
end program test_interpolate_clough_tocher_cache
