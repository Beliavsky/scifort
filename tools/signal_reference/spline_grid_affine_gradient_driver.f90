! SPDX-License-Identifier: MIT
! Machine-readable SciPy comparison of grid gradients and affine Jacobians.
program spline_grid_affine_gradient_driver
    use scifort_kinds, only : dp
    use scifort_signal_spline_nd, only : spline_interpolation_plan, &
        prepare_spline_plan_nd, evaluate_spline_plan_affine_with_gradient_nd, &
        evaluate_spline_plan_grid_with_gradient_nd
    implicit none
    character(len=16), parameter :: modes(8)=[character(len=16) :: &
        'mirror','constant','grid-constant','nearest', &
        'reflect','grid-mirror','wrap','grid-wrap']
    integer, parameter :: degrees(4)=[0,1,3,5]
    type(spline_interpolation_plan) :: plan
    real(dp) :: x(17,16),a(2,2),off(2),pos(2,6,5),r(6,5),rg(2,6,5),rdm(2,2,30),rdo(2,30)
    complex(dp) :: z(17,16),cz(6,5),cg(2,6,5),cdm(2,2,30),cdo(2,30)
    integer :: i,j,k,mi,di,kind,icase,code

    do j=1,16
        do i=1,17
            x(i,j)=sin(0.19_dp*i+0.23_dp*j)+0.003_dp*i*j
        end do
    end do
    z=cmplx(x,0.31_dp*x-0.17_dp,dp)
    a(1,:)=[0.48_dp,0.19_dp]
    a(2,:)=[-0.15_dp,0.31_dp]
    do mi=1,8
        do di=1,4
            do kind=0,1
                if (kind==0) then
                    call prepare_spline_plan_nd(x,plan,order=degrees(di), &
                        mode=trim(modes(mi)),cval=0.43_dp,status=code)
                else
                    call prepare_spline_plan_nd(z,plan,order=degrees(di), &
                        mode=trim(modes(mi)),cval=cmplx(0.43_dp,-0.29_dp,dp),status=code)
                end if
                if (code/=0) error stop 'prepare'
                do icase=1,2
                    if (icase==1) then
                        off=[1.13_dp,1.07_dp]
                    else
                        off=[-2.73_dp,13.68_dp]
                    end if
                    do j=1,5
                        do i=1,6
                            pos(:,i,j)=1.0_dp+off+matmul(a,real([i-1,j-1],dp))
                        end do
                    end do
                    if (kind==0) then
                        call evaluate_spline_plan_grid_with_gradient_nd(plan,pos,r,rg,code)
                        if (code/=0) error stop 'grid real'
                        call evaluate_spline_plan_affine_with_gradient_nd(plan,a,off,r,rdm,rdo,code, &
                            batch_size=4)
                        if (code/=0) error stop 'affine real'
                        do j=1,5
                            do i=1,6
                                k=i+6*(j-1)
                                write(*,'(6(i0,1x),*(es25.16,1x))') mi,degrees(di),kind,icase,i,j, &
                                    r(i,j),0.0_dp,rdm(1,1,k),0.0_dp,rdm(2,1,k),0.0_dp, &
                                    rdm(1,2,k),0.0_dp,rdm(2,2,k),0.0_dp, &
                                    rdo(1,k),0.0_dp,rdo(2,k),0.0_dp
                            end do
                        end do
                    else
                        call evaluate_spline_plan_grid_with_gradient_nd(plan,pos,cz,cg,code)
                        if (code/=0) error stop 'grid complex'
                        call evaluate_spline_plan_affine_with_gradient_nd(plan,a,off,cz,cdm,cdo,code, &
                            batch_size=4)
                        if (code/=0) error stop 'affine complex'
                        do j=1,5
                            do i=1,6
                                k=i+6*(j-1)
                                write(*,'(6(i0,1x),*(es25.16,1x))') mi,degrees(di),kind,icase,i,j, &
                                    real(cz(i,j)),aimag(cz(i,j)), &
                                    real(cdm(1,1,k)),aimag(cdm(1,1,k)), &
                                    real(cdm(2,1,k)),aimag(cdm(2,1,k)), &
                                    real(cdm(1,2,k)),aimag(cdm(1,2,k)), &
                                    real(cdm(2,2,k)),aimag(cdm(2,2,k)), &
                                    real(cdo(1,k)),aimag(cdo(1,k)),real(cdo(2,k)),aimag(cdo(2,k))
                            end do
                        end do
                    end if
                end do
            end do
        end do
    end do
end program spline_grid_affine_gradient_driver
