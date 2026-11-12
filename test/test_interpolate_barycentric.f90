! SPDX-License-Identifier: MIT
program test_interpolate_barycentric
    use scifort_kinds, only : dp
    use scifort_interpolate, only : barycentric_interpolator, barycentric_interpolate, &
        barycentric_interpolator_type, complex_barycentric_interpolator_type
    implicit none
    type(barycentric_interpolator_type) :: p
    type(complex_barycentric_interpolator_type) :: z
    real(dp) :: x(4), y(4,2), q, s(4)
    real(dp), allocatable :: yvals(:,:), weights(:), nodes(:)
    complex(dp) :: yc(4)
    integer :: i,code
    x=[-1.0_dp,0.2_dp,1.2_dp,2.0_dp]
    y(:,1)=x**3-2.0_dp*x**2+1.0_dp
    y(:,2)=2.0_dp*y(:,1)+x
    p=barycentric_interpolator(x,y,status=code)
    if (.not.p%valid() .or. code/=0) error stop 'barycentric construction'
    do i=0,30
        q=-1.2_dp+real(i,dp)*0.12_dp
        if (abs(p%value(q)-(q**3-2*q**2+1))>5.0e-12_dp) error stop 'barycentric value'
        if (abs(p%derivative(q)-(3*q**2-4*q))>1.0e-11_dp) error stop 'barycentric first derivative'
        if (abs(p%derivative(q,nu=2)-(6*q-4))>1.0e-10_dp) error stop 'barycentric second derivative'
        if (abs(p%derivative(q,nu=3)-6.0_dp)>1.0e-10_dp) error stop 'barycentric third derivative'
        if (abs(p%value(q,channel=2)-(2*(q**3-2*q**2+1)+q))>5.0e-12_dp) &
            error stop 'barycentric channel 2'
    end do
    yvals=p%values([0.25_dp,0.75_dp],nu=2)
    if (any(abs(yvals(:,2)-2.0_dp*(6.0_dp*[0.25_dp,0.75_dp]-4.0_dp))>1.0e-10_dp)) &
        error stop 'matrix derivatives'
    weights=p%weights()
    if (size(weights)/=4 .or. abs(maxval(abs(weights))-1.0_dp)>1.0e-14_dp) &
        error stop 'normalized weights'
    call p%add_nodes([0.25_dp],reshape([0.25_dp**3-2*0.25_dp**2+1, &
        2*(0.25_dp**3-2*0.25_dp**2+1)+0.25_dp],[1,2]),status=code)
    if (code/=0 .or. abs(p%value(0.7_dp)-(0.7_dp**3-2*0.7_dp**2+1))>1.0e-11_dp) &
        error stop 'add nodes'
    call p%add_nodes([0.25_dp],reshape([0.0_dp,0.0_dp],[1,2]),status=code)
    if (code==0) error stop 'duplicate node accepted'
    nodes=p%nodes()
    if (size(nodes)/=5) error stop 'transactional failed add'
    yc=cmplx(y(:,1),3.0_dp*y(:,1),dp)
    z=barycentric_interpolator(x,yc,status=code)
    if (.not.z%valid() .or. code/=0) error stop 'complex constructor'
    if (abs(z%value(0.5_dp)-cmplx(0.5_dp**3-2*0.5_dp**2+1, &
        3*(0.5_dp**3-2*0.5_dp**2+1),dp))>1.0e-12_dp) error stop 'complex value'
    if (abs(z%derivative(0.7_dp,nu=3)-cmplx(6.0_dp,18.0_dp,dp))>1.0e-10_dp) &
        error stop 'complex derivative'
    if (size(barycentric_interpolate(x,y(:,1),[0.7_dp]))/=1) error stop 'convenience interface'
    p=barycentric_interpolator([0.0_dp,0.0_dp],[1.0_dp,2.0_dp],status=code)
    if (code==0 .or. p%valid()) error stop 'duplicate input not rejected'
    p=barycentric_interpolator([1.0_dp],[2.0_dp],status=code)
    if (code/=0 .or. p%derivative(0.1_dp)/=0.0_dp) error stop 'constant polynomial'
    print '(a)', 'barycentric interpolation native regression PASS'
end program test_interpolate_barycentric
