! SPDX-License-Identifier: MIT
! Regression for accelerated 3D simplex location at interior points and faces.
program test_interpolate_tetra_bvh
    use scifort_kinds, only: dp
    use scifort_interpolate, only: linear_nd_interpolator, &
        linear_nd_interpolator_type, complex_linear_nd_interpolator_type
    implicit none

    integer, parameter :: nsite=90, nchannels=3
    real(dp) :: x(nsite,3),y(nsite,nchannels),q(3),weights(4)
    real(dp), allocatable :: v(:), expected(:), batch(:,:),queries(:,:)
    complex(dp) :: cy(nsite,nchannels)
    complex(dp), allocatable :: cv(:)
    integer, allocatable :: simplices(:,:),seed(:)
    integer :: i,j,nt,nseed,code
    type(linear_nd_interpolator_type) :: p
    type(complex_linear_nd_interpolator_type) :: cp

    call random_seed(size=nseed)
    allocate(seed(nseed))
    seed=20261031
    call random_seed(put=seed)
    call random_number(x)
    x(:,1)=10.0_dp+2.0_dp*x(:,1)
    x(:,2)=-20.0_dp+3.0_dp*x(:,2)
    x(:,3)=50.0_dp+4.0_dp*x(:,3)
    y(:,1)=sin(x(:,1))*cos(2.0_dp*x(:,3))
    y(:,2)=cos(20.0_dp*x(:,2))+0.1_dp*x(:,1)
    y(:,3)=sin(x(:,1)+x(:,3))
    cy=cmplx(y,-0.5_dp*y,dp)

    p=linear_nd_interpolator(x,y,rescale=.true.,fill_value=-42.0_dp,status=code)
    if (code/=0 .or. .not.p%valid()) error stop 'tetra BVH real constructor'
    simplices=p%simplices()
    nt=size(simplices,2)
    if (size(simplices,1)/=4 .or. nt<100) error stop 'tetra BVH simplices'
    allocate(queries(nt,3))
    weights=[0.1_dp,0.2_dp,0.3_dp,0.4_dp]

    ! Interior points of every simplex: they must interpolate its four
    ! vertices with exactly the same barycentric weights as the original scan.
    do j=1,nt
        q=matmul(weights,x(simplices(:,j),:))
        expected=matmul(weights,y(simplices(:,j),:))
        v=p%evaluate(q,status=code)
        if (code/=0 .or. maxval(abs(v-expected))>2.0e-8_dp) &
            error stop 'tetra BVH interior interpolation'
        queries(j,:)=q
    end do
    batch=p%evaluate(queries,status=code)
    if (code/=0) error stop 'tetra BVH batch status'
    do j=1,nt
        v=p%evaluate(queries(j,:))
        if (maxval(abs(v-batch(j,:)))>1.0e-12_dp) &
            error stop 'tetra BVH batch/single consistency'
    end do

    ! Points extremely close to shared faces and edges are important. Move
    ! slightly inward to avoid conditioning-sensitive exact-boundary rounding
    ! present even in the predecessor's full-scan implementation.
    do j=1,nt
        q=(1.0_dp-1.0e-7_dp)*sum(x(simplices(1:3,j),:),dim=1)/3.0_dp + &
          1.0e-7_dp*x(simplices(4,j),:)
        expected=(1.0_dp-1.0e-7_dp)*sum(y(simplices(1:3,j),:),dim=1)/3.0_dp + &
          1.0e-7_dp*y(simplices(4,j),:)
        v=p%evaluate(q)
        if (maxval(abs(v-expected))>2.0e-8_dp) then
            error stop 'tetra BVH shared-face interpolation'
        end if
        q=(0.5_dp-1.0e-7_dp)*(x(simplices(1,j),:)+x(simplices(4,j),:)) + &
          1.0e-7_dp*(x(simplices(2,j),:)+x(simplices(3,j),:))
        expected=(0.5_dp-1.0e-7_dp)*(y(simplices(1,j),:)+y(simplices(4,j),:)) + &
          1.0e-7_dp*(y(simplices(2,j),:)+y(simplices(3,j),:))
        v=p%evaluate(q)
        if (maxval(abs(v-expected))>2.0e-8_dp) &
            error stop 'tetra BVH shared-edge interpolation'
    end do
    do i=1,nsite
        v=p%evaluate(x(i,:))
        if (maxval(abs(v-y(i,:)))>2.0e-8_dp) error stop 'tetra BVH vertex interpolation'
    end do
    q=[-30.0_dp,-30.0_dp,40.0_dp]
    v=p%evaluate(q,status=code)
    if (code/=0 .or. any(v/=-42.0_dp)) error stop 'tetra BVH outside fill'

    cp=linear_nd_interpolator(x,cy,rescale=.true.,fill_value=cmplx(-7.0_dp,2.0_dp,dp),status=code)
    if (code/=0 .or. .not.cp%valid()) error stop 'tetra BVH complex constructor'
    do j=1,nt,7
        q=matmul(weights,x(simplices(:,j),:))
        cv=cp%evaluate(q)
        if (maxval(abs(cv-cmplx(matmul(weights,y(simplices(:,j),:)), &
                    -0.5_dp*matmul(weights,y(simplices(:,j),:)),dp)))>2.0e-8_dp) &
            error stop 'tetra BVH complex interpolation'
    end do
    cv=cp%evaluate([-30.0_dp,-30.0_dp,40.0_dp])
    if (any(abs(cv-cmplx(-7.0_dp,2.0_dp,dp))>1.0e-12_dp)) &
        error stop 'tetra BVH complex outside fill'

    print '(a,i0,a)', '3D indexed tetrahedron location: ',nt,' simplex, face, edge and complex regressions PASS'
end program test_interpolate_tetra_bvh
