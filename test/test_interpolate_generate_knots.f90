! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_interpolate_generate_knots
    use scifort_kinds, only: dp
    use scifort_interpolate, only: generate_knots,knot_generator_type,make_lsq_spline, &
        bspline_interpolator_type
    use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
    implicit none
    type(knot_generator_type) :: gen,copy,empty
    type(bspline_interpolator_type) :: model
    real(dp) :: x(41),y(41),z(41,2),w(41),nan,inf,previous,fp
    real(dp), allocatable :: t(:),other(:),last(:)
    integer :: i,k,code,count,oldn,checks,j
    logical :: done,finished
    checks=0
    nan=ieee_value(0.0_dp,ieee_quiet_nan)
    inf=ieee_value(0.0_dp,ieee_positive_inf)
    x=[(real(i-1,dp)/40.0_dp,i=1,41)]
    w=0.5_dp+x
    y=sin(25.0_dp*x)+0.2_dp*cos(71.0_dp*x)
    y(41)=y(1)
    z(:,1)=y
    z(:,2)=cos(16.0_dp*x)
    z(41,:)=z(1,:)
    do k=1,5
        gen=generate_knots(x,y,k=k,status=code)
        call require(code==0,'interpolation initialization')
        call gen%next(t,done,code)
        call require(code==0 .and. .not.done,'single interpolation vector')
        call require(size(t)==size(x)+k+1,'not-a-knot count')
        call require(all(t(:k+1)==x(1)) .and. all(t(42:)==x(41)),'clamped endpoints')
        model=make_lsq_spline(x,y,t,k=k,status=code)
        call require(code==0,'interpolation knot validity')
        call require(maxval(abs(model%evaluate(x)-y))<1e-10_dp,'interpolation recovery')
        call gen%next(t,done,code)
        call require(code==0 .and. done .and. size(t)==0,'interpolation exhausted')
        gen=generate_knots(x(:k+1),y(:k+1),k=k,s=1e-280_dp,status=code)
        call require(code==0,'minimal sample count initialized')
        call gen%next(t,done,code)
        call require(code==0 .and. .not.done,'minimal sample count first vector')
        call gen%next(t,done,code)
        call require(code==0 .and. done,'minimal sample count terminates')
        gen=generate_knots(x,y,k=k,bc_type='periodic',status=code)
        call gen%next(t,done,code)
        call require(code==0 .and. .not.done,'periodic interpolation vector')
        if (k==1) then
            call require(t(1)==t(2) .and. t(size(t))==t(size(t)-1), &
                'linear zero-target periodic clamped convention')
        else
            j=size(t)-2*k-1
            call require(maxval(abs(t(j+1:)-t(:size(t)-j)-1.0_dp))<2e-15_dp, &
                'periodic interpolation cyclic convention')
        end if

        gen=generate_knots(x,y,k=k,s=0.1_dp,w=w,status=code)
        call require(code==0,'weighted initialization')
        copy=gen
        oldn=0
        previous=huge(1.0_dp)
        count=0
        do
            call gen%next(t,done,code)
            call require(code==0,'weighted next status')
            call copy%next(other,finished,code)
            call require(code==0 .and. (done.eqv.finished),'independent state copy')
            if (done) exit
            call require(size(t)==size(other),'copy vector length')
            call require(all(t==other),'copy vector identity')
            call require(size(t)>oldn,'strictly increasing sequence sizes')
            oldn=size(t)
            count=count+1
            call require(count<=size(x),'finite sequence')
            model=make_lsq_spline(x,y,t,k=k,w=w,status=code)
            call require(code==0,'generated weighted LSQ knots')
            fp=sum((w*(model%evaluate(x)-y))**2)
            call require(fp<=previous+1e-10_dp,'residual monotonicity')
            previous=fp
        end do
        call require(fp<0.1001_dp,'final weighted residual')

        gen=generate_knots(x,z,k=k,s=0.001_dp,bc_type='periodic',status=code)
        call require(code==0,'shared periodic initialization')
        oldn=0
        do
            call gen%next(t,done,code)
            call require(code==0,'periodic next status')
            if (done) exit
            call require(size(t)>oldn,'periodic sizes increase')
            oldn=size(t)
            j=size(t)-2*k-1
            call require(maxval(abs(t(j+1:)-t(:size(t)-j)-1.0_dp))<2e-15_dp, &
                'canonical periodic knot cycles')
            last=t
        end do
        call require(size(last)<=size(x)+2*k,'periodic count bound')
        gen=generate_knots(x,y,k=k,s=1e8_dp,bc_type='periodic',status=code)
        call require(code==0,'constant-limit initialization')
        call gen%next(t,done,code)
        call require(code==0 .and. .not.done .and. size(t)==2*(k+1),'periodic constant limit')
        call gen%next(t,done,code)
        call require(code==0 .and. done,'constant limit exhausted')
        gen=generate_knots(x,y,k=k,s=1e-20_dp,nest=2*(k+1),status=code)
        call gen%next(t,done,code)
        call require(code==0 .and. size(t)==2*(k+1),'soft count initial vector')
        call gen%next(t,done,code)
        call require(code==0 .and. size(t)==2*(k+1)+1,'soft count may exceed initial capacity')
        call gen%next(t,done,code)
        call require(code==0 .and. done,'soft count exhausted')
    end do
    gen=generate_knots(x,2.0_dp+x,s=1e-10_dp,status=code)
    call gen%next(t,done,code)
    call require(code==0 .and. size(t)==8,'polynomial initial knots')
    call gen%next(t,done,code)
    call require(code==0 .and. done,'polynomial requires no extra knots')
    w=1.0_dp
    w(41)=1e300_dp
    gen=generate_knots(x,y,w=w,s=0.1_dp,bc_type='periodic',status=code)
    copy=generate_knots(x,y,s=0.1_dp,bc_type='periodic',status=j)
    call require(code==0 .and. j==0,'finite closing weight validated then ignored')
    do
        call gen%next(t,done,code)
        call copy%next(other,finished,j)
        call require(code==0 .and. j==0 .and. (done.eqv.finished),'ignored last weight sequence')
        if (done) exit
        call require(size(t)==size(other),'ignored last weight knot count')
        call require(all(t==other),'ignored last weight knots')
    end do
    call empty%next(t,done,code)
    call require(code==1 .and. done .and. size(t)==0,'default state rejected')
    gen=generate_knots([(real(i,dp),i=1,10001)],[(real(i,dp),i=1,10001)],status=code)
    call require(code==3,'sample size limit')
    gen=generate_knots(x,spread(y,2,257),status=code)
    call require(code==3,'channel size limit')
    gen=generate_knots(x,y*1e150_dp,s=1e-100_dp,status=code)
    call require(code==3,'unrepresentable normalized target')
    gen=generate_knots(x*1e-310_dp,y,status=code)
    call require(code==3,'subnormal interval rejected')
    w=1.0_dp
    w(2)=1e-100_dp
    gen=generate_knots(x,y,w=w,s=1.0_dp,status=code)
    call require(code==3,'excessive active weight range rejected')
    gen=generate_knots(1e15_dp+0.125_dp*[(real(i,dp),i=0,4)],y(:5),k=2,status=code)
    call require(code==0,'rounded midpoint case initializes')
    call gen%next(t,done,code)
    call require(code==3 .and. done .and. size(t)==0,'collapsed midpoint output rejected')
    gen=generate_knots(x,y,s=-1.0_dp,status=code)
    call require(code==1,'negative target')
    gen=generate_knots(x,y,s=nan,status=code)
    call require(code==1,'nan target')
    gen=generate_knots(x,y,s=inf,status=code)
    call require(code==1,'infinite target')
    gen=generate_knots(x,y,k=0,status=code)
    call require(code==1,'degree zero')
    gen=generate_knots(x,y,k=6,status=code)
    call require(code==1,'degree six')
    gen=generate_knots(x,y,w=w,status=code)
    call require(code==1,'interpolation weights rejected')
    gen=generate_knots(x,y,nest=9,status=code)
    call require(code==1,'interpolation capacity rejected')
    gen=generate_knots(x,y,s=1.0_dp,nest=7,status=code)
    call require(code==1,'small capacity rejected')
    gen=generate_knots(x,y,s=1.0_dp,bc_type='invalid',status=code)
    call require(code==1,'unknown boundary')
    gen=generate_knots(x,y,s=1.0_dp,xb=0.1_dp,status=code)
    call require(code==1,'nonenclosing bounds')
    gen=generate_knots(x,y,s=1.0_dp,xb=-1.0_dp,bc_type='periodic',status=code)
    call require(code==1,'extended periodic bounds')
    gen=generate_knots(x,y(:40),s=1.0_dp,status=code)
    call require(code==1,'shape mismatch')
    w(2)=0.0_dp
    gen=generate_knots(x,y,s=1.0_dp,w=w,status=code)
    call require(code==1,'zero weight rejected')
    w(2)=nan
    gen=generate_knots(x,y,s=1.0_dp,w=w,status=code)
    call require(code==1,'nan weight')
    w(2)=inf
    gen=generate_knots(x,y,s=1.0_dp,w=w,status=code)
    call require(code==1,'infinite weight')
    y(41)=2.0_dp
    gen=generate_knots(x,y,s=1.0_dp,bc_type='periodic',status=code)
    call require(code==1,'periodic closure required')
    y(2)=nan
    gen=generate_knots(x,y,status=code)
    call require(code==1,'nan observation')
    call gen%next(t,done,code)
    call require(code==1 .and. done .and. size(t)==0,'failed initialization leaves unusable state')
    y(2)=inf
    gen=generate_knots(x,y,status=code)
    call require(code==1,'infinite observation')
    x(2)=x(1)
    gen=generate_knots(x,z,status=code)
    call require(code==1,'duplicate sites')
    x(2)=nan
    gen=generate_knots(x,z,status=code)
    call require(code==1,'nan site')
    write(*,'(a,i0,a)') 'generate_knots: ',checks,' checks passed'
contains
    subroutine require(ok,label)
        implicit none
        logical, intent(in) :: ok !! Whether the tested contract holds.
        character(*), intent(in) :: label !! Contract description printed on failure.
        checks=checks+1
        if (.not.ok) then
            write(*,'(a)') 'FAIL: '//label
            stop 1
        end if
    end subroutine require
end program test_interpolate_generate_knots
