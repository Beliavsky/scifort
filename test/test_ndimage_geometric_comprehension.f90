! SPDX-License-Identifier: MIT
! Tests geometric callback coordinates and indexed label callback grouping.
program test_ndimage_geometric_comprehension
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_ndimage, only : geometric_transform_nd, map_coordinates_nd, &
        labeled_comprehension_nd, ndimage_status_too_expensive
    implicit none
    real(dp) :: image(3,4),out(5,4),direct(5,4),coords(2,20),expected(20)
    real(dp) :: stretched(10,8),reduced(4),costs(6),out15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: real15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1),large(270000)
    real(dp) :: nan_value
    complex(dp) :: z(3,4),zout(5,4),zexpected(20)
    integer :: labels(3,4),label15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1),requests(6)
    integer :: i,j,k,n,status,callback_count,map_count
    real(dp), parameter :: tol=1.0e-10_dp

    do j=1,4
        do i=1,3
            image(i,j)=sin(0.2_dp*i+0.7_dp*j)+real(i*j,dp)*0.08_dp
            labels(i,j)=modulo(i+2*j,5)-2
        end do
    end do
    z=cmplx(image,0.5_dp*image+0.4_dp,dp)
    map_count=0
    call geometric_transform_nd(image,map2,out,order=1,status=status,mode='nearest')
    call check(status==0,'geometric real call')
    call check(map_count==size(out),'mapping called once per destination')
    do j=1,4
        do i=1,5
            k=(j-1)*5+i
            call map2([i,j],coords(:,k))
        end do
    end do
    call map_coordinates_nd(image,coords,expected,order=1,status=status,mode='nearest')
    call check(status==0,'reference interpolation')
    direct=reshape(expected,shape(direct))
    call check(maxval(abs(out-direct))<tol,'geometric points and order')
    call geometric_transform_nd(z,map2,zout,order=3,status=status,mode='grid-constant', &
        cval=cmplx(-0.7_dp,0.3_dp,dp))
    call check(status==0,'complex geometry')
    call map_coordinates_nd(z,coords,zexpected,order=3,status=status, &
        mode='grid-constant',prefilter=.true.,cval=cmplx(-0.7_dp,0.3_dp,dp))
    call check(status==0,'complex interpolation reference')
    call check(maxval(abs(zout-reshape(zexpected,shape(zout))))<tol,'complex geometry values')
    stretched=-999.0_dp
    call geometric_transform_nd(image,map2,stretched(1:10:2,1:8:2), &
        order=0,status=status,mode='nearest')
    call check(status==0,'noncontiguous geometry')
    call map_coordinates_nd(image,coords,expected,order=0,status=status,mode='nearest')
    call check(maxval(abs(stretched(1:10:2,1:8:2)-reshape(expected,[5,4])))<tol, &
        'noncontiguous geometric output')
    call check(all(stretched(2:10:2,:)==-999.0_dp),'untouched strided output')

    call geometric_transform_nd(image,map2,out,status=status,order=-1)
    call check(status/=0,'invalid interpolation order')
    call geometric_transform_nd(image,map2,out,status=status,mode='invalid')
    call check(status/=0,'invalid interpolation mode')
    call geometric_transform_nd(image,bad_mapping,out,status=status,order=1)
    call check(status/=0,'reject NaN mapping coordinates')
    call geometric_transform_nd(real15,map15,large,status=status,order=0)
    call check(status==ndimage_status_too_expensive,'coordinate resource limit')

    real15=reshape([3.0_dp,7.0_dp],shape(real15))
    call geometric_transform_nd(real15,map15,out15,order=1,status=status)
    call check(status==0,'rank-15 input')
    call check(abs(out15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)-3.0_dp)<tol, &
        'rank-15 first pixel')
    call check(abs(out15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)-7.0_dp)<tol, &
        'rank-15 second pixel')

    requests=[-1,2,2,0,99,-1]
    callback_count=0
    call labeled_comprehension_nd(image,labels,requests,ordered_reduce,-8.0_dp,costs,status)
    call check(status==0,'labeled comprehension reduction')
    call check(callback_count==count([(count(labels==requests(k))>0,k=1,6)]), &
        'callback invoked per requested label including duplicates')
    do k=1,size(requests)
        n=0
        expected=0.0_dp
        do i=1,3
            do j=1,4
                if (labels(i,j)/=requests(k)) cycle
                n=n+1
                expected(n)=image(i,j)
            end do
        end do
        if (n==0) then
            call check(costs(k)==-8.0_dp,'missing label default')
        else
            call check(abs(costs(k)-ordered_reduce(expected(:n)))<tol, &
                'grouped C-order callback values')
        end if
    end do
    call check(abs(costs(1)-costs(6))<tol,'repeated label values')
    call check(abs(costs(2)-costs(3))<tol,'repeated label calls')
    call labeled_comprehension_nd(image,labels,[integer::],ordered_reduce,0.0_dp, &
        reduced(:0),status)
    call check(status==0,'empty requested index')
    call labeled_comprehension_nd(image,labels,requests,ordered_reduce,0.0_dp, &
        reduced,status)
    call check(status/=0,'invalid output shape')
    nan_value=ieee_value(0.0_dp,ieee_quiet_nan)
    image(1,1)=nan_value
    call labeled_comprehension_nd(image,labels,requests,ordered_reduce,0.0_dp, &
        costs,status)
    call check(status/=0,'nonfinite image rejected')

    label15=33
    call labeled_comprehension_nd(real15,label15,[33],ordered_reduce,0.0_dp, &
        costs(:1),status)
    call check(status==0,'rank-15 label callback')
    call check(abs(costs(1)-17.0_dp)<tol,'rank-15 C-order samples')
    print '(a)', 'PASS: geometric transforms and batched label callbacks'
contains
    subroutine check(ok,message)
        logical, intent(in) :: ok !! assertion condition
        character(len=*), intent(in) :: message !! assertion description
        if (.not.ok) then
            print '(2a)', 'FAIL: ',message
            error stop 1
        end if
    end subroutine check

    subroutine map2(output_coords,input_coords)
        integer, intent(in) :: output_coords(:) !! two-dimensional output coordinate
        real(dp), intent(out) :: input_coords(:) !! source interpolation location
        map_count=map_count+1
        if (map_count<=size(out)) then
            call check(all(output_coords==[1+(map_count-1)/4, &
                1+modulo(map_count-1,4)]),'mapping invocation C-order')
        end if
        input_coords(1)=real(output_coords(1),dp)+0.35_dp*real(output_coords(2)-3,dp)
        input_coords(2)=real(output_coords(2),dp)-0.2_dp*real(output_coords(1)-2,dp)
    end subroutine map2

    subroutine map15(output_coords,input_coords)
        integer, intent(in) :: output_coords(:) !! rank-15 output point
        real(dp), intent(out) :: input_coords(:) !! rank-15 source point
        input_coords=real(output_coords,dp)
    end subroutine map15

    subroutine bad_mapping(output_coords,input_coords)
        integer, intent(in) :: output_coords(:) !! unused coordinates
        real(dp), intent(out) :: input_coords(:) !! invalid coordinate
        input_coords=ieee_value(0.0_dp,ieee_quiet_nan)
    end subroutine bad_mapping

    function ordered_reduce(samples) result(value)
        real(dp), intent(in) :: samples(:) !! C-ordered samples
        real(dp) :: value !! order-sensitive weighted sum
        integer :: a
        value=0.0_dp
        do a=1,size(samples)
            value=value+real(a,dp)*samples(a)
        end do
        callback_count=callback_count+1
    end function ordered_reduce
end program test_ndimage_geometric_comprehension
