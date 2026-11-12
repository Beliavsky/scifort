! SPDX-License-Identifier: MIT
program test_interpolate_rbf
    use scifort_kinds, only:dp
    use scifort_interpolate, only: rbf_interpolator,rbf_interpolator_type, &
        complex_rbf_interpolator_type,rbf_status_success,rbf_status_invalid_input,rbf_status_singular
    implicit none
    type(rbf_interpolator_type) :: p
    type(complex_rbf_interpolator_type) :: pc
    real(dp) :: xy(9,2),f(9,3),z(2),q(4,2),expected(3),testmat(3,3)
    complex(dp) :: complex_f(9,3),expected_complex(3)
    real(dp), allocatable :: a(:,:)
    complex(dp), allocatable :: ca(:,:)
    integer :: i,code,j
    character(16) :: names(3)
    do i=1,9
        xy(i,1)=cos(real(i,dp))
        xy(i,2)=sin(real(i,dp))+0.05_dp*real(i,dp)
        f(i,1)=1.0_dp+2.0_dp*xy(i,1)-3.0_dp*xy(i,2)
        f(i,2)=5.0_dp-1.2_dp*xy(i,1)+4.0_dp*xy(i,2)
        f(i,3)=-0.3_dp+0.7_dp*xy(i,1)+0.8_dp*xy(i,2)
    end do
    q(1,:)=[0.1_dp,0.3_dp]
    q(2,:)=[-0.2_dp,-0.5_dp]
    q(3,:)=[0.85_dp,-0.4_dp]
    q(4,:)=[-1.2_dp,1.0_dp]
    p=rbf_interpolator(xy,f,kernel='thin_plate_spline',status=code)
    if (code/=rbf_status_success .or. .not.p%valid()) error stop 'RBF global fit'
    if (p%dimensions()/=2 .or. p%channels()/=3) error stop 'RBF metadata'
    a=p%evaluate(q,status=code)
    if (code/=0 .or. any(shape(a)/=[4,3])) error stop 'RBF shape'
    do i=1,4
        expected=[1.0_dp+2*q(i,1)-3*q(i,2), &
                  5.0_dp-1.2_dp*q(i,1)+4*q(i,2), &
                  -0.3_dp+0.7_dp*q(i,1)+0.8_dp*q(i,2)]
        if (maxval(abs(a(i,:)-expected))>2.0e-10_dp) error stop 'RBF polynomial reproduction'
        if (abs(p%value(q(i,:),channel=2)-expected(2))>2.0e-10_dp) error stop 'RBF scalar value'
    end do
    p=rbf_interpolator(xy,f,neighbors=7,degree=1,kernel='cubic',status=code)
    if (code/=0) error stop 'RBF local construction'
    a=p%evaluate(q,status=code)
    if (code/=0) error stop 'RBF local query'
    do i=1,4
        expected=[1.0_dp+2*q(i,1)-3*q(i,2), &
                  5.0_dp-1.2_dp*q(i,1)+4*q(i,2), &
                  -0.3_dp+0.7_dp*q(i,1)+0.8_dp*q(i,2)]
        if (maxval(abs(a(i,:)-expected))>2.0e-10_dp) error stop 'RBF local polynomial'
    end do
    complex_f=cmplx(f,2.0_dp*f,dp)
    pc=rbf_interpolator(xy,complex_f,kernel='thin_plate_spline',status=code)
    if (code/=0 .or. .not.pc%valid() .or. pc%channels()/=3) error stop 'RBF complex fit'
    ca=pc%evaluate(q,status=code)
    if (code/=0 .or. any(shape(ca)/=[4,3])) error stop 'RBF complex shape'
    do i=1,4
        expected=[1.0_dp+2*q(i,1)-3*q(i,2), &
                  5.0_dp-1.2_dp*q(i,1)+4*q(i,2), &
                  -0.3_dp+0.7_dp*q(i,1)+0.8_dp*q(i,2)]
        expected_complex=cmplx(expected,2.0_dp*expected,dp)
        if (maxval(abs(ca(i,:)-expected_complex))>2.0e-10_dp) error stop 'RBF complex polynomial'
    end do
    if (abs(pc%value(q(1,:),channel=3)-ca(1,3))>1.0e-12_dp) error stop 'RBF complex value'
    ! Separate SciFort-only Matern support, checked against the 2x2 closed form.
    names=[character(len=16) :: 'matern1_2','matern3_2','matern5_2']
    do j=1,3
        block
            real(dp) :: sites(2,1), d(2), w, k01, k0, k1, denom, truth
            sites(:,1)=[0.0_dp,1.0_dp]
            d=[2.0_dp,3.0_dp]
            select case(j)
            case(1)
                k01=exp(-1.0_dp)
                k0=exp(-0.4_dp)
                k1=exp(-0.6_dp)
            case(2)
                w=sqrt(3.0_dp)
                k01=(1.0_dp+w)*exp(-w)
                k0=(1.0_dp+0.4_dp*w)*exp(-0.4_dp*w)
                k1=(1.0_dp+0.6_dp*w)*exp(-0.6_dp*w)
            case(3)
                w=sqrt(5.0_dp)
                k01=(1.0_dp+w+5.0_dp/3.0_dp)*exp(-w)
                k0=(1.0_dp+0.4_dp*w+5.0_dp*0.4_dp**2/3.0_dp)*exp(-0.4_dp*w)
                k1=(1.0_dp+0.6_dp*w+5.0_dp*0.6_dp**2/3.0_dp)*exp(-0.6_dp*w)
            end select
            denom=1.1_dp**2-k01*k01
            truth=(k0*(1.1_dp*2.0_dp-k01*3.0_dp) &
                  +k1*(1.1_dp*3.0_dp-k01*2.0_dp))/denom
            p=rbf_interpolator(sites,d,kernel=trim(names(j)),epsilon=1.0_dp, &
                degree=-1,smoothing=0.1_dp,status=code)
            if (code/=0 .or. abs(p%value([0.4_dp])-truth)>1.0e-12_dp) &
                error stop 'RBF Matern analytically verified'
        end block
    end do
    p=rbf_interpolator(xy,f,kernel='gaussian',status=code)
    if (code/=rbf_status_invalid_input .or. p%valid()) error stop 'missing epsilon accepted'
    p=rbf_interpolator(xy,f,kernel='nonexistent',status=code)
    if (code/=rbf_status_invalid_input .or. p%valid()) error stop 'invalid kernel accepted'
    p=rbf_interpolator(xy,f,smoothing=-1.0_dp,status=code)
    if (code/=rbf_status_invalid_input .or. p%valid()) error stop 'negative smoothing accepted'
    p=rbf_interpolator(xy,f,degree=4,status=code)
    if (code/=rbf_status_invalid_input .or. p%valid()) error stop 'too few points accepted'
    block
        real(dp) :: bad_x(4,1)
        bad_x(:,1)=[0.0_dp,0.0_dp,1.0_dp,2.0_dp]
        p=rbf_interpolator(bad_x,[1.0_dp,1.0_dp,2.0_dp,3.0_dp], &
            kernel='cubic',status=code)
        if (code/=rbf_status_singular .or. p%valid()) error stop 'duplicate exact data accepted'
    end block
    print '(a)', 'RBFInterpolator global/local/complex/Matern/validation PASS'
end program test_interpolate_rbf
