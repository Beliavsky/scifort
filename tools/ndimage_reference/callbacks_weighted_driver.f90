! SPDX-License-Identifier: MIT
! SciPy comparison driver: axis-weighted distances and callback filtering.
program callbacks_weighted_driver
    use scifort_kinds, only : dp
    use scifort_ndimage, only : correlate1d_nd, convolve1d_nd, &
        generic_filter_nd, labeled_comprehension_nd, distance_transform_cdt_nd, &
        distance_transform_weighted_cdt_nd
    implicit none
    character(len=14), parameter :: modes(8)=[character(len=14) :: &
        'reflect','grid-mirror','mirror','nearest','wrap','grid-wrap','constant','grid-constant']
    real(dp) :: a(5,6),a3(4,3,2),out(5,6),out3(4,3,2)
    complex(dp) :: z(5,6),zout(5,6)
    real(dp) :: weights(4),answers(6)
    complex(dp) :: cweights(3)
    integer :: maskdist(5,6),features(2,30),maskdist3(4,3,2),features3(3,24)
    logical :: mask(5,6),mask3(4,3,2),foot(2,3),foot3(3,2,2)
    integer :: labels(5,6),labels3(4,3,2),ids(6),st,iu,i,j,k,m,ax,w,org,what,n
    integer :: cost2(2),cost3(3),stencil2(3,3),stencil3(3,3,3)

    open(newunit=iu,file='callbacks_weighted_results.bin',form='unformatted', &
         access='stream',status='replace')
    do j=1,6
        do i=1,5
            a(i,j)=sin(0.3_dp*i+0.7_dp*j)+i-0.4_dp*j
            z(i,j)=cmplx(a(i,j),cos(0.7_dp*i-0.2_dp*j),dp)
            labels(i,j)=mod(i+3*j,7)-3
            mask(i,j)=mod(i*i+2*j,7)/=0
        end do
    end do
    do k=1,2
        do j=1,3
            do i=1,4
                a3(i,j,k)=cos(0.6_dp*i-j+2*k)+0.7_dp*i*j-k
                labels3(i,j,k)=mod(i+2*j+3*k,6)-2
                mask3(i,j,k)=mod(i+2*j+4*k,5)/=0
            end do
        end do
    end do
    weights=[1.3_dp,-0.4_dp,0.2_dp,0.7_dp]
    cweights=[cmplx(0.4_dp,0.7_dp,dp),cmplx(1.1_dp,-0.3_dp,dp),cmplx(-0.8_dp,0.9_dp,dp)]
    do m=1,8
        do ax=1,2
            do w=2,4
                do org=-(w/2),(w-1)/2
                    call correlate1d_nd(a,weights(:w),ax,out,st,trim(modes(m)), -1.25_dp,org)
                    if (st/=0) error stop 'real correlate1d'
                    write(iu) out
                    call convolve1d_nd(a,weights(:w),ax,out,st,trim(modes(m)), -1.25_dp,org)
                    if (st/=0) error stop 'real convolve1d'
                    write(iu) out
                end do
            end do
        end do
        do ax=1,3
            do w=2,3
                do org=-(w/2),(w-1)/2
                    call correlate1d_nd(a3,weights(:w),ax,out3,st,trim(modes(m)), 1.75_dp,org)
                    if (st/=0) error stop '3d correlate1d'
                    write(iu) out3
                    call convolve1d_nd(a3,weights(:w),ax,out3,st,trim(modes(m)), 1.75_dp,org)
                    if (st/=0) error stop '3d convolve1d'
                    write(iu) out3
                end do
            end do
        end do
        do ax=1,2
            do w=2,3
                do org=-(w/2),(w-1)/2
                    call correlate1d_nd(z,cweights(:w),ax,zout,st,trim(modes(m)), &
                        cmplx(-1.25_dp,0.3_dp,dp),org)
                    if (st/=0) error stop 'complex correlate1d'
                    write(iu) zout
                    call convolve1d_nd(z,cweights(:w),ax,zout,st,trim(modes(m)), &
                        cmplx(-1.25_dp,0.3_dp,dp),org)
                    if (st/=0) error stop 'complex convolve1d'
                    write(iu) zout
                end do
            end do
        end do
    end do
    foot=reshape([.true.,.false.,.true.,.true.,.false.,.true.],[2,3])
    foot3=reshape([.true.,.false.,.true.,.true.,.true.,.false.,.false., &
         .true.,.true.,.false.,.true.,.false.],[3,2,2])
    do m=1,8
        do what=1,3
            select case(what)
            case(1)
                call generic_filter_nd(a,weighted_sum,out,st,sizes=[2,3], &
                    origin=[-1,1],mode=trim(modes(m)),cval=-2.5_dp)
            case(2)
                call generic_filter_nd(a,weighted_sum,out,st,footprint=foot, &
                    origin=[0,-1],mode=trim(modes(m)),cval=-2.5_dp)
            case(3)
                call generic_filter_nd(a,weighted_sum,out,st,sizes=[3,3], &
                    origin=[1,-1],mode=trim(modes(m)),cval=-2.5_dp)
            end select
            if (st/=0) error stop 'generic 2D'
            write(iu) out
        end do
        call generic_filter_nd(a3,weighted_sum,out3,st,footprint=foot3, &
             origin=[-1,0,-1],mode=trim(modes(m)),cval=1.125_dp)
        if (st/=0) error stop 'generic 3D'
        write(iu) out3
    end do
    ids=[-2,0,3,11,-2,2]
    do what=1,3
        select case(what)
        case(1)
            call labeled_comprehension_nd(a,labels,ids,reduce_max,-11.0_dp,answers,st)
        case(2)
            call labeled_comprehension_nd(a,labels,ids,reduce_var,-11.0_dp,answers,st)
        case(3)
            call labeled_comprehension_nd(a,labels,ids,weighted_sum,-11.0_dp,answers,st)
        end select
        if (st/=0) error stop 'labeled 2D'
        write(iu) answers
    end do
    do what=1,3
        select case(what)
        case(1)
            call labeled_comprehension_nd(a3,labels3,ids,reduce_max,-11.0_dp,answers,st)
        case(2)
            call labeled_comprehension_nd(a3,labels3,ids,reduce_var,-11.0_dp,answers,st)
        case(3)
            call labeled_comprehension_nd(a3,labels3,ids,weighted_sum,-11.0_dp,answers,st)
        end select
        if (st/=0) error stop 'labeled 3D'
        write(iu) answers
    end do
    do w=1,5
        select case(w)
        case(1)
            cost2=[1,1]
            cost3=[1,1,1]
        case(2)
            cost2=[2,3]
            cost3=[1,3,5]
        case(3)
            cost2=[3,2]
            cost3=[4,2,1]
        case(4)
            cost2=[1,6]
            cost3=[6,1,2]
        case(5)
            cost2=[4,5]
            cost3=[2,4,3]
        end select
        call distance_transform_cdt_nd(mask,maskdist,st,axis_weights=cost2,feature_indices=features)
        if (st/=0) error stop 'weighted chamfer 2D'
        write(iu) maskdist,features
        call distance_transform_cdt_nd(mask3,maskdist3,st,axis_weights=cost3,feature_indices=features3)
        if (st/=0) error stop 'weighted chamfer 3D'
        write(iu) maskdist3,features3
    end do
    do w=1,5
        stencil2=0
        select case(w)
        case(1)
            stencil2(1,2)=1
            stencil2(3,2)=1
            stencil2(2,1)=1
            stencil2(2,3)=1
        case(2)
            stencil2=1
            stencil2(2,2)=0
        case(3)
            stencil2=4
            stencil2(2,2)=0
            stencil2(1,2)=3
            stencil2(3,2)=3
            stencil2(2,1)=3
            stencil2(2,3)=3
        case(4)
            stencil2(1,2)=3
            stencil2(3,2)=2
            stencil2(2,1)=4
            stencil2(2,3)=1
        case(5)
            stencil2(1,1)=1
            stencil2(1,3)=2
            stencil2(3,1)=3
            stencil2(3,3)=4
        end select
        call distance_transform_weighted_cdt_nd(mask,stencil2,maskdist,st,features)
        if (st/=0) error stop 'stencil chamfer 2D'
        write(iu) maskdist,features
    end do
    stencil3=0
    do k=1,3
        do j=1,3
            do i=1,3
                n=abs(i-2)+abs(j-2)+abs(k-2)
                if (n>0) stencil3(i,j,k)=n
            end do
        end do
    end do
    call distance_transform_weighted_cdt_nd(mask3,stencil3,maskdist3,st,features3)
    if (st/=0) error stop 'stencil chamfer 3D'
    write(iu) maskdist3,features3
    close(iu)
contains
    real(dp) function weighted_sum(values) result(y)
        real(dp), intent(in) :: values(:)
        integer :: l
        y=0.0_dp
        do l=1,size(values)
            y=y+real(l,dp)*values(l)
        end do
    end function weighted_sum

    real(dp) function reduce_max(values) result(y)
        real(dp), intent(in) :: values(:)
        y=maxval(values)
    end function reduce_max

    real(dp) function reduce_var(values) result(y)
        real(dp), intent(in) :: values(:)
        real(dp) :: ave
        ave=sum(values)/real(size(values),dp)
        y=sum((values-ave)**2)/real(size(values),dp)
    end function reduce_var
end program callbacks_weighted_driver
