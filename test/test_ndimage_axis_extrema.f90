! SPDX-License-Identifier: MIT
! Focused native regression for axis-only uniform/minimum/maximum filters.
program test_ndimage_axis_extrema
    use scifort_kinds, only : dp
    use scifort_ndimage, only : uniform_filter1d_nd, minimum_filter1d_nd, &
        maximum_filter1d_nd, minimum_filter_nd, maximum_filter_nd, &
        correlate1d_nd, ndimage_status_invalid_input, ndimage_status_too_expensive
    implicit none
    real(dp) :: x(5,6),y(5,6),ref(5,6),t(10,12)
    complex(dp) :: z(5,6),w(5,6),zref(5,6)
    real(dp) :: rank15(3,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: out15(3,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    integer :: i,j,s
    do j=1,6
        do i=1,5
            x(i,j)=cos(0.3_dp*i+0.5_dp*j)+real(i*j,dp)/14.0_dp
        end do
    end do
    z=cmplx(x,2*x,dp)
    call uniform_filter1d_nd(x,4,2,y,status=s,mode='grid-wrap',origin=-1)
    call correlate1d_nd(x,[.25_dp,.25_dp,.25_dp,.25_dp],2,ref, &
        status=i,mode='grid-wrap',origin=-1)
    if (s/=0 .or. i/=0 .or. maxval(abs(y-ref))>1.e-14_dp) &
        error stop 'even uniform 1d'
    call uniform_filter1d_nd(z,4,2,w,status=s,mode='grid-constant', &
        cval=cmplx(-1.0_dp,3.0_dp,dp),origin=1)
    call correlate1d_nd(z,[cmplx(.25_dp,0.0_dp,dp),cmplx(.25_dp,0.0_dp,dp), &
        cmplx(.25_dp,0.0_dp,dp),cmplx(.25_dp,0.0_dp,dp)],2,zref, &
        status=i,mode='grid-constant',cval=cmplx(-1.0_dp,3.0_dp,dp),origin=1)
    if (s/=0 .or. i/=0 .or. maxval(abs(w-zref))>1.e-14_dp) &
        error stop 'complex uniform 1d'
    call minimum_filter1d_nd(x,3,1,y,status=s,mode='nearest')
    call minimum_filter_nd(x,[3,1],ref,status=i,mode='nearest')
    if (s/=0 .or. i/=0 .or. any(y/=ref)) error stop 'axis minimum baseline'
    call maximum_filter1d_nd(x,3,2,y,status=s,mode='nearest')
    call maximum_filter_nd(x,[1,3],ref,status=i,mode='nearest')
    if (s/=0 .or. i/=0 .or. any(y/=ref)) error stop 'axis maximum baseline'
    call minimum_filter1d_nd(x,2,1,t(1:10:2,1:12:2),status=s,origin=0)
    call minimum_filter1d_nd(x,2,1,ref,status=i,origin=0)
    if (s/=0 .or. i/=0 .or. any(t(1:10:2,1:12:2)/=ref)) &
        error stop 'strided output 1d'
    call maximum_filter1d_nd(x,0,1,y,status=s)
    if (s/=ndimage_status_invalid_input) error stop 'zero-width window'
    call minimum_filter1d_nd(x,3,0,y,status=s)
    if (s/=ndimage_status_invalid_input) error stop 'zero-based axis'
    call maximum_filter1d_nd(x,4,1,y,status=s,origin=2)
    if (s/=ndimage_status_invalid_input) error stop 'invalid window origin'
    call uniform_filter1d_nd(x,4,2,y,status=s,mode='bad')
    if (s/=ndimage_status_invalid_input) error stop 'invalid window mode'
    call uniform_filter1d_nd(x,50001,2,y,status=s)
    if (s/=ndimage_status_too_expensive) error stop 'uniform resource guard'
    call minimum_filter1d_nd(x,50001,2,y,status=s)
    if (s/=ndimage_status_too_expensive) error stop 'minimum resource guard'
    rank15=0.0_dp
    rank15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=1.0_dp
    call maximum_filter1d_nd(rank15,3,1,out15,status=s)
    if (s/=0 .or. any(abs(out15-1.0_dp)>1.e-14_dp)) &
        error stop 'rank15 maxima'
    call uniform_filter1d_nd(rank15,1,15,out15,status=s)
    if (s/=0 .or. any(abs(out15-rank15)>1.e-14_dp)) &
        error stop 'rank15 singleton axis'
    print '(a)', 'PASS: 1D uniform/minimum/maximum filters'
end program test_ndimage_axis_extrema
