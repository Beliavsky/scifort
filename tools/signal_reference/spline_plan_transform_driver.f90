! SPDX-License-Identifier: MIT
! Stream affine samples for cross-language SciPy comparison.
program spline_plan_transform_driver
    use scifort_kinds, only : dp
    use scifort_signal_spline_nd, only : spline_interpolation_plan, &
        prepare_spline_plan_nd, evaluate_spline_plan_affine_nd, &
        evaluate_spline_plan_geometric_nd
    implicit none
    character(len=16), parameter :: modes(8)=[character(len=16) :: &
        'mirror','constant','grid-constant','nearest', &
        'reflect','grid-mirror','wrap','grid-wrap']
    real(dp) :: source(19,17),real_result(12,11),matrix(2,2),offset(2)
    complex(dp) :: z(19,17),complex_result(12,11)
    type(spline_interpolation_plan) :: plan
    integer :: i,j,degree,mi,code

    do j=1,17
        do i=1,19
            source(i,j)=sin(0.17_dp*i+0.31_dp*j)+0.004_dp*i*j
        end do
    end do
    z=cmplx(source,-0.37_dp*source+0.29_dp,dp)
    matrix=reshape([1.109_dp,0.187_dp,-0.123_dp,0.943_dp],[2,2])
    offset=[-2.4131_dp,1.2373_dp]
    do mi=1,size(modes)
        do degree=0,5
            call prepare_spline_plan_nd(source,plan,order=degree, &
                mode=trim(modes(mi)),cval=1.125_dp,status=code)
            if (code/=0) error stop 'real prepare failed'
            call evaluate_spline_plan_affine_nd(plan,matrix,offset,real_result,code,batch_size=13)
            if (code/=0) error stop 'real affine failed'
            do j=1,11
                do i=1,12
                    write(*,'(5(i0,1x),2(es25.16e3,1x))') mi,degree,0,i,j, &
                        real_result(i,j),0.0_dp
                end do
            end do
            call prepare_spline_plan_nd(z,plan,order=degree, &
                mode=trim(modes(mi)),cval=cmplx(1.125_dp,-0.35_dp,dp),status=code)
            if (code/=0) error stop 'complex prepare failed'
            call evaluate_spline_plan_geometric_nd(plan,affine_mapping,complex_result,code,batch_size=37)
            if (code/=0) error stop 'complex geometric failed'
            do j=1,11
                do i=1,12
                    write(*,'(5(i0,1x),2(es25.16e3,1x))') mi,degree,1,i,j, &
                        real(complex_result(i,j),dp),aimag(complex_result(i,j))
                end do
            end do
        end do
    end do
contains
    subroutine affine_mapping(out_coord,in_coord)
        integer, intent(in) :: out_coord(:) !! one-based output pixel
        real(dp), intent(out) :: in_coord(:) !! one-based source position
        in_coord=1.0_dp+offset+matmul(matrix,real(out_coord-1,dp))
    end subroutine affine_mapping
end program spline_plan_transform_driver
