! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_interpolate_aaa
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan,ieee_value, &
        ieee_quiet_nan,ieee_positive_inf
    use scifort_kinds, only : dp
    use scifort_interpolate, only : aaa,aaa_type,complex_aaa_type,aaa_success, &
        aaa_invalid_input,aaa_max_terms,aaa_undefined, &
        aaa_numerical_failure
    implicit none
    type(aaa_type) :: r,other
    type(complex_aaa_type) :: z
    real(dp) :: x(101),f(101),q(99),nan,inf,val
    real(dp), allocatable :: values(:),points(:),history(:),weights(:)
    complex(dp) :: cx(101),cf(101),cv,exact_pole
    complex(dp), allocatable :: poles(:),roots(:),res(:),cvalues(:)
    integer :: i,info,removed,checks
    checks=0
    nan=ieee_value(0.0_dp,ieee_quiet_nan)
    inf=ieee_value(0.0_dp,ieee_positive_inf)
    x=[(-2.0_dp+real(i-1,dp)/25.0_dp,i=1,101)]
    q=[(-1.97_dp+real(i-1,dp)*3.94_dp/98.0_dp,i=1,99)]
    f=x/(1.0_dp+x*x)
    r=aaa(x,f,status=info)
    call check(r%valid() .and. r%converged(),'rational fit is valid and converged')
    call check(info==aaa_success,'automatic cleanup successful')
    values=r%values(q,info)
    call check(info==0 .and. maxval(abs(values-q/(1+q*q)))<1e-12_dp,'rational evaluation')
    points=r%support_points()
    call check(size(points)==3,'three-term rational representation')
    call check(maxval(abs(r%values(points)-r%support_values()))==0.0_dp,'exact support values')
    history=r%errors()
    call check(size(history)==3 .and. history(3)<1e-12_dp,'error history')
    call check(r%max_error()<1e-12_dp .and. r%fit_status()==0,'current sample error')
    poles=r%poles(info)
    call check(info==0 .and. size(poles)==2,'two finite poles')
    call check(minval(abs(poles-cmplx(0.0_dp,1.0_dp,dp)))<1e-10_dp,'positive imaginary pole')
    call check(minval(abs(poles-cmplx(0.0_dp,-1.0_dp,dp)))<1e-10_dp,'negative imaginary pole')
    res=r%residues(info)
    call check(info==0 .and. maxval(abs(res-0.5_dp))<1e-10_dp,'analytic residues')
    roots=r%roots(info)
    call check(info==0 .and. size(roots)==1,'one finite zero')
    call check(abs(roots(1))<1e-10_dp,'analytic zero')
    cv=r%value(cmplx(0.2_dp,0.7_dp,dp),info)
    call check(info==0 .and. abs(cv-cmplx(0.2_dp,0.7_dp,dp)/ &
        (1+cmplx(0.2_dp,0.7_dp,dp)**2))<1e-12_dp,'complex query of real fit')
    cv=r%value(cmplx(1e-5_dp,1.0_dp,dp),info)
    call check(info==0 .and. abs(cv/(cmplx(1e-5_dp,1.0_dp,dp)/ &
        (1+cmplx(1e-5_dp,1.0_dp,dp)**2))-1)<1e-8_dp,'near-pole evaluation')
    do i=1,99
        call check(abs(values(i)-r%value(q(i)))<1e-15_dp,'scalar/vector agreement')
    end do
    val=r%value(inf,info)
    call check(info==0 .and. abs(val)<1e-12_dp,'finite limit at infinity')
    val=r%value(nan,info)
    call check(info==aaa_invalid_input .and. ieee_is_nan(val),'NaN query')
    val=other%value(0.0_dp,info)
    call check(info==aaa_invalid_input .and. ieee_is_nan(val),'uninitialized query')

    f=exp(x)
    r=aaa(x,f,status=info,rtol=1e-11_dp)
    call check(r%valid() .and. r%converged(),'exponential fit')
    points=r%support_points()
    call check(all(r%values(points)==r%support_values()),'exact stored support samples')
    values=r%values(q)
    call check(maxval(abs(values-exp(q)))<2e-10_dp,'exponential accuracy')
    r=aaa(x,f,status=info,max_terms=2,clean_up=.false.,rtol=1e-14_dp)
    call check(info==aaa_max_terms .and. r%valid() .and. .not.r%converged(),'iteration budget')

    exact_pole=cmplx(0.3_dp,0.4_dp,dp)
    cx=exp(cmplx(0.0_dp,acos(-1.0_dp)*(x+2.0_dp)/2.0_dp,dp))
    ! Avoid the repeated periodic endpoint, while retaining complex-coordinate coverage.
    cf=(cx+0.2_dp)/(cx-exact_pole)
    z=aaa(cx(:100),cf(:100),status=info)
    call check(z%valid() .and. z%converged(),'complex-domain fit')
    poles=z%poles(info)
    call check(info==0 .and. size(poles)==1,'complex finite pole')
    call check(abs(poles(1)-exact_pole)<1e-10_dp,'complex pole location')
    roots=z%roots(info)
    call check(info==0 .and. size(roots)==1 .and. abs(roots(1)+0.2_dp)<1e-10_dp,'complex zero')
    res=z%residues(info)
    call check(info==0 .and. abs(res(1)-exact_pole-0.2_dp)<1e-10_dp,'complex residue')
    cvalues=z%values(cmplx(q,0.5_dp,dp),info)
    call check(info==0 .and. maxval(abs(cvalues-(cmplx(q,0.5_dp,dp)+0.2_dp)/ &
        (cmplx(q,0.5_dp,dp)-exact_pole)))<1e-10_dp,'complex vector values')
    z=aaa(x,cmplx(f,0.5_dp,dp),status=info,clean_up=.false.)
    call check(z%valid() .and. z%converged(),'real coordinates, complex values')
    z=aaa(cx(:100),real(cf(:100),dp),status=info,max_terms=4,clean_up=.false.)
    call check(z%valid(),'complex coordinates, real values')

    r=aaa([0.0_dp],[7.0_dp],status=info)
    call check(info==0 .and. r%value(22.0_dp)==7.0_dp,'one sample is constant')
    poles=r%poles(info)
    call check(info==0 .and. size(poles)==0,'constant has no poles')
    roots=r%roots(info)
    call check(info==0 .and. size(roots)==0,'nonzero constant has no zeros')
    r=aaa(x,0.0_dp*x,status=info)
    call check(info==0 .and. r%value(1.0_dp)==0.0_dp,'zero function')
    roots=r%roots(info)
    call check(info==aaa_undefined .and. size(roots)==0,'identically zero roots undefined')
    r=aaa([1.0_dp,0.0_dp,1.0_dp],[3.0_dp,2.0_dp,99.0_dp],status=info)
    call check(r%converged() .and. r%value(1.0_dp)==3.0_dp,'first duplicate value retained')
    r=aaa([0.0_dp,1.0_dp,2.0_dp],[2.0_dp,nan,inf],status=info)
    call check(info==0 .and. r%value(0.2_dp)==2.0_dp,'nonfinite values discarded')
    r=aaa([0.0_dp,1.0_dp],[nan,inf],status=info)
    call check(info==aaa_invalid_input .and. .not.r%valid(),'all values nonfinite')
    r=aaa([0.0_dp,nan],[1.0_dp,inf],status=info)
    call check(info==aaa_invalid_input,'nonfinite coordinate rejected before filtering')
    r=aaa(x,f,status=info,rtol=-1.0_dp)
    call check(info==aaa_invalid_input,'negative tolerance')
    r=aaa(x,f,status=info,rtol=nan)
    call check(info==aaa_invalid_input,'NaN tolerance')
    r=aaa(x,f,status=info,max_terms=0)
    call check(info==aaa_invalid_input,'invalid support budget')
    r=aaa(x,f(:100),status=info)
    call check(info==aaa_invalid_input,'shape mismatch')
    r=aaa(x(:0),f(:0),status=info)
    call check(info==aaa_invalid_input,'empty input')
    r=aaa(x,f,status=info,clean_up_tol=-1.0_dp)
    call check(info==aaa_invalid_input,'invalid cleanup threshold')
    r=aaa([-1.0_dp,0.0_dp,1e-320_dp,1.0_dp],[-1.0_dp,0.0_dp,1.0_dp,2.0_dp], &
        status=info,rtol=0.0_dp,clean_up=.false.)
    call check(info==aaa_numerical_failure,'unrepresentable Loewner quotient rejected')
    call check(r%fit_status()==aaa_numerical_failure,'failed-fit diagnostic retained')
    r=aaa([-1.0_dp,0.0_dp,0.5_dp,1.0_dp],[1.0_dp,0.0_dp,1e-320_dp,2e-320_dp], &
        status=info,clean_up=.false.)
    call check(info==0 .and. r%valid() .and. r%converged(),'subnormal Loewner column scaling')
    val=r%value(0.5_dp,info)
    call check(info==0 .and. abs(val-1e-320_dp)<1e-12_dp,'finite equilibrated evaluation')

    r=aaa(1e100_dp*x,1e150_dp*x/(1+x*x),status=info)
    call check(r%valid() .and. r%converged(),'large coordinate and data scales')
    val=r%value(0.2e100_dp,info)
    call check(info==0 .and. abs(val/(1e150_dp*0.2_dp/1.04_dp)-1)<1e-10_dp, &
        'scaled physical evaluation')
    r=aaa(1e-100_dp*x,1e-150_dp*x/(1+x*x),status=info)
    call check(r%valid() .and. r%converged(),'small coordinate and data scales')
    val=r%value(0.2e-100_dp,info)
    call check(info==0 .and. abs(val/(1e-150_dp*0.2_dp/1.04_dp)-1)<1e-10_dp, &
        'small-scale evaluation')
    r=aaa(1e-200_dp*x,1e-200_dp*x/(1+x*x),status=info)
    call check(info==0 .and. r%converged(),'cleanup independent of physical residue underflow')
    f=(1.3e-320_dp*x)/1.3e-320_dp
    r=aaa(1.3e-320_dp*x,1e308_dp*(f/(1+f*f)),status=info)
    call check(info==0 .and. r%converged(),'subnormal coordinates with large data')
    res=r%residues(info)
    call check(info==0 .and. size(res)==2,'subnormal-scale residues available')
    call check(maxval(abs(res/(0.5_dp*(1.3e-320_dp*1e308_dp))-1))<1e-8_dp, &
        'combined residue scale avoids intermediate loss')
    f=(1e6_dp+1e-6_dp*x-1e6_dp)*1e6_dp
    r=aaa(1e6_dp+1e-6_dp*x,f/(1+f*f),status=info,clean_up=.false.)
    val=r%value(1e6_dp+0.2e-6_dp,info)
    cv=cmplx((1e6_dp+0.2e-6_dp-1e6_dp)*1e6_dp,0.0_dp,dp)
    call check(info==0 .and. abs(val-real(cv/(1+cv*cv),dp))<1e-10_dp, &
        'closely spaced coordinates with a large offset')

    f=sin(1000.0_dp*x)
    r=aaa(x(:41),f(:41),rtol=0.0_dp,max_terms=41,clean_up=.false.,status=info)
    call check(info==0 .and. r%converged(),'large retained representation')
    points=r%support_points()
    call check(size(points)>39,'large retained support set')
    poles=r%poles(info)
    call check(info==0 .and. size(poles)==40,'spectrum beyond former size limit')
    removed=r%clean_up(status=info)
    call check(info==0 .and. r%valid(),'large cleanup succeeds and preserves fit')
    call check(r%cleanup_status()==0,'stored successful cleanup diagnostic')
    other=aaa(x(:41),f(:41),rtol=0.0_dp,max_terms=41,status=info)
    call check(info==0 .and. other%valid() .and. other%converged(), &
        'automatic cleanup supports larger converged fits')

    f=1.0_dp/(x+3.0_dp)+1e-10_dp/(x-0.3_dp)
    other=aaa(x,f,status=info,rtol=1e-14_dp,clean_up_tol=1e-6_dp)
    call check(info==aaa_max_terms .and. other%valid() .and. .not.other%converged(), &
        'automatic cleanup reports final tolerance loss')
    r=aaa(x,f,status=info,rtol=1e-14_dp,clean_up=.false.)
    call check(r%valid() .and. r%converged(),'weak-pole fit')
    removed=r%clean_up(1e-6_dp,info)
    call check(info==0 .and. removed>=1,'weak-pole cleanup')
    poles=r%poles(info)
    call check(info==0 .and. size(poles)==1 .and. abs(poles(1)+3.0_dp)<1e-6_dp, &
        'dominant pole retained after cleanup')
    call check(r%max_error()<1e-7_dp .and. .not.r%converged(),'cleanup updates actual error')
    points=r%support_points()
    removed=r%clean_up(nan,info)
    call check(info==aaa_invalid_input .and. removed==0,'invalid manual cleanup')
    call check(all(r%support_points()==points),'failed cleanup preserves representation')

    print '(a,i0,a)','AAA: ',checks,' checks passed.'
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
end program test_interpolate_aaa
