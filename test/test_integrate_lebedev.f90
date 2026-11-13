! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_integrate_lebedev
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_integrate, only : dp,lebedev_rule,integrate_status_success, &
        integrate_status_invalid_argument
    implicit none
    integer, parameter :: orders(*)=[3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,35, &
        41,47,53,59,65,71,77,83,89,95,101,107,113,119,125,131]
    integer, parameter :: sizes(*)=[6,14,26,38,50,74,86,110,146,170,194,230,266,302,350, &
        434,590,770,974,1202,1454,1730,2030,2354,2702,3074,3470,3890,4334,4802,5294,5810]
    real(dp), parameter :: pi=acos(-1.0_dp),eps=epsilon(1.0_dp)
    real(dp), allocatable :: x(:,:),w(:),powers(:,:,:),terms(:)
    real(dp) :: target,actual,bound,error_before,error_after,previous
    integer :: j,n,m,a,b,c,k,info,checks,moments,axis
    complex(dp) :: integral
    checks=0
    moments=0
    do j=1,size(orders)
        n=orders(j)
        m=sizes(j)
        call lebedev_rule(n,x,w,info)
        call check(info==integrate_status_success,'success status')
        call check(all(shape(x)==[3,m]) .and. size(w)==m,'array shapes')
        call check(all(ieee_is_finite(x)) .and. all(ieee_is_finite(w)),'finite data')
        call check(maxval(abs(sum(x*x,dim=1)-1.0_dp))<=4.0_dp*eps,'unit sphere')
        call check(abs(sum(w)-4.0_dp*pi)<5e-13_dp,'constant surface integral')
        select case(n)
        case(13,25,27)
            call check(any(w<0.0_dp),'documented signed weights')
        case default
            call check(all(w>0.0_dp),'positive weights on other rules')
        end select
        do axis=1,3
            call check(abs(sum(w*x(axis,:)))<32.0_dp*eps,'odd first moment')
            call check(abs(sum(w*x(axis,:)**n))<32.0_dp*eps,'highest odd axial moment')
        end do
        allocate(powers(0:n/2,3,m),terms(m))
        powers(0,:,:)=1.0_dp
        do k=1,n/2
            powers(k,:,:)=powers(k-1,:,:)*x*x
        end do
        ! All inequivalent even monomials; coordinate symmetry checked independently in Python.
        do a=0,n/2
            do b=0,a
                do c=0,b
                    if (a+b+c>n/2) cycle
                    target=moment(a,b,c)
                    terms=w*powers(a,1,:)*powers(b,2,:)*powers(c,3,:)
                    actual=sum(terms)
                    bound=3e-13_dp*target+16.0_dp*eps*sum(abs(terms))
                    call check(abs(actual-target)<=bound,'even polynomial moment')
                    moments=moments+1
                end do
            end do
        end do
        ! Nonsymmetric assignments exercise coordinate permutations at full degree.
        a=n/4
        b=n/2-a
        do axis=1,3
            k=mod(axis,3)+1
            actual=sum(w*powers(a,axis,:)*powers(b,k,:))
            target=moment(a,b,0)
            call check(abs(actual-target)<4e-13_dp*target,'permuted highest-degree moment')
        end do
        deallocate(powers,terms)
    end do
    call lebedev_rule(3,x,w)
    call check(all(abs(sum(x,dim=2))==0.0_dp),'six axial nodes sum to zero')
    call check(all(abs(w-2.0_dp*pi/3.0_dp)<4.0_dp*eps*(2.0_dp*pi/3.0_dp)), &
        'lowest-rule weights')
    call lebedev_rule(5,x,w,info)
    call check(abs(sum(w*x(1,:)**2*x(2,:)**2)-4.0_dp*pi/15.0_dp)<2e-15_dp, &
        'mixed fourth moment reference')
    previous=huge(1.0_dp)
    error_before=0.0_dp
    error_after=0.0_dp
    do n=3,17,2
        call lebedev_rule(n,x,w,info)
        actual=sum(w*exp(x(1,:)))
        target=4.0_dp*pi*sinh(1.0_dp)
        if (n==3) error_before=abs(actual-target)
        if (n==17) error_after=abs(actual-target)
        call check(abs(actual-target)<=previous+1e-14_dp,'analytic exponential convergence')
        previous=abs(actual-target)
    end do
    call check(error_after<1e-13_dp .and. error_before>1e-3_dp,'nonpolynomial accuracy')
    integral=sum(w*cmplx(x(1,:)**2,x(2,:)**2,dp))
    call check(abs(integral-cmplx(4.0_dp*pi/3.0_dp,4.0_dp*pi/3.0_dp,dp))<1e-13_dp, &
        'complex-valued integration')
    ! A second call replaces old allocated arrays, including on invalid input.
    do n=-2,132
        if (any(orders==n)) cycle
        call lebedev_rule(n,x,w,info)
        call check(info==integrate_status_invalid_argument,'unsupported order')
        call check(.not.allocated(x) .and. .not.allocated(w),'no output on invalid order')
    end do
    call lebedev_rule(huge(1),x,w,info)
    call check(info==integrate_status_invalid_argument,'large invalid integer')
    call lebedev_rule(131,x,w,info)
    call check(info==0 .and. size(w)==5810,'valid call after invalid input')
    call lebedev_rule(4,x,w)
    call check(.not.allocated(x) .and. .not.allocated(w),'invalid call without status')
    print '(a,i0,a,i0,a)','Lebedev: ',checks,' checks, ',moments,' polynomial moments passed.'
contains
    subroutine check(ok,label)
        implicit none
        logical, intent(in) :: ok !! assertion result
        character(*), intent(in) :: label !! assertion context
        if (.not.ok) then
            print *, 'FAIL: ',label
            error stop 1
        end if
        checks=checks+1
    end subroutine check
    pure function moment(a,b,c) result(value)
        implicit none
        integer, intent(in) :: a,b,c !! half of the nonnegative even Cartesian exponents
        real(dp) :: value
        integer :: i
        value=4.0_dp*pi
        do i=1,a
            value=value*real(2*i-1,dp)/real(2*i+1,dp)
        end do
        do i=1,b
            value=value*real(2*i-1,dp)/real(2*(a+i)+1,dp)
        end do
        do i=1,c
            value=value*real(2*i-1,dp)/real(2*(a+b+i)+1,dp)
        end do
    end function moment
end program test_integrate_lebedev
