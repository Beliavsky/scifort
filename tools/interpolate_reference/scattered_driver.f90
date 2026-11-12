! SPDX-License-Identifier: MIT
! Stream protocol for independent NumPy/SciPy scattered-data comparisons.
program scattered_driver
    use scifort_kinds, only: dp
    use scifort_interpolate_scattered, only: nearest_nd_interpolator, linear_nd_interpolator, griddata, &
        nearest_nd_interpolator_type, complex_nearest_nd_interpolator_type, &
        linear_nd_interpolator_type, complex_linear_nd_interpolator_type
    implicit none
    integer :: method,n,d,c,m,rescale,complex_data,i,j,code
    real(dp) :: fr,fi,a,b
    real(dp), allocatable :: x(:,:),q(:,:),yr(:,:),out_r(:,:)
    complex(dp), allocatable :: z(:,:),out_c(:,:)
    type(nearest_nd_interpolator_type) :: nr
    type(complex_nearest_nd_interpolator_type) :: nc
    type(linear_nd_interpolator_type) :: lr
    type(complex_linear_nd_interpolator_type) :: lc
    read(*,*) method,n,d,c,m,rescale,complex_data
    read(*,*) fr,fi
    allocate(x(n,d),q(m,d),z(n,c),yr(n,c))
    do i=1,n
        read(*,*) x(i,:)
    end do
    do i=1,n
        do j=1,c
            read(*,*) a,b
            z(i,j)=cmplx(a,b,dp)
        end do
    end do
    yr=real(z,dp)
    do i=1,m
        read(*,*) q(i,:)
    end do
    if (method==2 .or. method==3) then
        if (complex_data==1) then
            if (method==2) then
                out_c=griddata(x,z,q,method='nearest',rescale=rescale==1,status=code)
            else
                out_c=griddata(x,z,q,method='linear',rescale=rescale==1, &
                    fill_value=cmplx(fr,fi,dp),status=code)
            end if
        else
            if (method==2) then
                out_r=griddata(x,yr,q,method='nearest',rescale=rescale==1,status=code)
            else
                out_r=griddata(x,yr,q,method='linear',rescale=rescale==1, &
                    fill_value=fr,status=code)
            end if
        end if
    else if (method==0) then
        if (complex_data==1) then
            nc=nearest_nd_interpolator(x,z,rescale=rescale==1,status=code)
            if (code==0) out_c=nc%evaluate(q,status=code)
        else
            nr=nearest_nd_interpolator(x,yr,rescale=rescale==1,status=code)
            if (code==0) out_r=nr%evaluate(q,status=code)
        end if
    else
        if (complex_data==1) then
            lc=linear_nd_interpolator(x,z,rescale=rescale==1, &
                fill_value=cmplx(fr,fi,dp),status=code)
            if (code==0) out_c=lc%evaluate(q,status=code)
        else
            lr=linear_nd_interpolator(x,yr,rescale=rescale==1,fill_value=fr,status=code)
            if (code==0) out_r=lr%evaluate(q,status=code)
        end if
    end if
    print '(i0)',code
    if (code/=0) stop
    do i=1,m
        do j=1,c
            if (complex_data==1) then
                print '(2(es26.17,1x))',real(out_c(i,j),dp),aimag(out_c(i,j))
            else
                print '(2(es26.17,1x))',out_r(i,j),0.0_dp
            end if
        end do
    end do
end program scattered_driver
