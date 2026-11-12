! SPDX-License-Identifier: MIT
! Test transport for randomized SciPy comparisons of Krogh and Floater-Hormann.
program global_interpolators_driver
    use scifort_kinds, only : dp
    use scifort_interpolate, only : krogh_interpolator, floater_hormann_interpolator, &
        krogh_interpolator_type, complex_krogh_interpolator_type, &
        floater_hormann_interpolator_type, complex_floater_hormann_interpolator_type
    implicit none
    integer :: ncases, case_id, algorithm, dtype, n, nch, nq, order, degree, code, nu
    real(dp), allocatable :: x(:), yr(:,:), yi(:,:), query(:), re(:), im(:)
    complex(dp), allocatable :: z(:,:), out(:,:), dall(:,:,:)
    type(krogh_interpolator_type) :: kp
    type(complex_krogh_interpolator_type) :: kz
    type(floater_hormann_interpolator_type) :: fp
    type(complex_floater_hormann_interpolator_type) :: fz
    read(*,*) ncases
    do case_id=1,ncases
        read(*,*) algorithm,dtype,n,nch,nq,order,degree
        allocate(x(n),yr(n,nch),yi(n,nch),query(nq),z(n,nch),re(nq*nch),im(nq*nch))
        read(*,*) x
        read(*,*) yr
        read(*,*) yi
        read(*,*) query
        z=cmplx(yr,yi,dp)
        code=1
        if (algorithm==1) then
            if (dtype==1) then
                kp=krogh_interpolator(x,yr,status=code)
            else
                kz=krogh_interpolator(x,z,status=code)
            end if
        else
            if (dtype==1) then
                fp=floater_hormann_interpolator(x,yr,d=degree,status=code)
            else
                fz=floater_hormann_interpolator(x,z,d=degree,status=code)
            end if
        end if
        if (code/=0) then
            print *, 'INVALID INPUT CASE',case_id
            error stop 1
        end if
        print '(a,1x,i0)', 'CASE',case_id
        do nu=0,order
            if (algorithm==1) then
                if (dtype==1) then
                    out=cmplx(kp%values(query,nu=nu),0.0_dp,dp)
                else
                    out=kz%values(query,nu=nu)
                end if
            else
                if (dtype==1) then
                    out=cmplx(fp%values(query),0.0_dp,dp)
                else
                    out=fz%values(query)
                end if
            end if
            print '(a,1x,i0)', 'ORDER',nu
            re=reshape(real(out,dp),[nq*nch])
            im=reshape(aimag(out),[nq*nch])
            print '(*(es26.17e3,1x))', re
            print '(*(es26.17e3,1x))', im
        end do
        if (algorithm==1) then
            block
                integer :: nder
                nder=min(n,4)
                if (dtype==1) then
                    dall=cmplx(kp%derivatives(query,der=nder),0.0_dp,dp)
                else
                    dall=kz%derivatives(query,der=nder)
                end if
                print '(a,1x,i0)', 'ALL',nder
                print '(*(es26.17e3,1x))', reshape(real(dall,dp),[size(dall)])
                print '(*(es26.17e3,1x))', reshape(aimag(dall),[size(dall)])
            end block
        end if
        deallocate(x,yr,yi,query,z,re,im)
    end do
end program global_interpolators_driver
