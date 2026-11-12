! SPDX-License-Identifier: MIT
! Selected-axis grayscale morphology and generic real/complex neighborhoods.
program test_ndimage_gray_generic_axes
    use scifort_kinds, only : dp
    use scifort_ndimage, only : grey_dilation_nd, grey_erosion_nd, &
        grey_opening_nd, generic_filter_nd, ndimage_status_success, &
        ndimage_status_invalid_input
    implicit none
    real(dp) :: x(4,5,3),a(4,5,3),b(4,5,3),large(8,5,3)
    real(dp) :: x15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,3)
    real(dp) :: y15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,3)
    complex(dp) :: cx(4,5,3),cy(4,5,3)
    integer :: i,j,k,status
    logical :: fp(2,3)

    do k=1,3
        do j=1,5
            do i=1,4
                x(i,j,k)=real(i*i-3*j+2*k,dp)
            end do
        end do
    end do
    fp=reshape([.true.,.false.,.true.,.true.,.false.,.true.],[2,3])
    cx=cmplx(x,.2_dp*x,dp)

    call grey_dilation_nd(x,a,status,sizes=[2,3],axes=[3,1])
    call require(status==ndimage_status_success,'selected grayscale dilation')
    call grey_dilation_nd(x,large(1:8:2,:,:),status,sizes=[2,3],axes=[3,1])
    call require(status==ndimage_status_success,'noncontiguous grayscale output')
    call require(all(a==large(1:8:2,:,:)),'noncontiguous dilation equality')
    call grey_erosion_nd(x,b,status,sizes=[integer::],axes=[integer::])
    call require(status==ndimage_status_success .and. all(x==b),'empty axes grayscale identity')

    call grey_opening_nd(x,a,status,footprint=fp,axes=[3,1], &
        mode='constant',cval=-3.0_dp)
    call require(status==ndimage_status_success,'nonseparable selected opening')

    call generic_filter_nd(x,largest,b,status,footprint=fp,axes=[3,1])
    call require(status==ndimage_status_success,'selected neighborhood callback')
    call generic_filter_nd(cx,largest_complex,cy,status,sizes=[2,3],axes=[3,1])
    call require(status==ndimage_status_success,'selected complex callback')

    x15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,:)=[1.0_dp,4.0_dp,2.0_dp]
    call grey_dilation_nd(x15,y15,status,sizes=[2],axes=[15],mode='nearest')
    call require(status==ndimage_status_success,'rank-15 grayscale dilation')
    call require(all(y15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,:)== &
        [4.0_dp,4.0_dp,2.0_dp]),'rank-15 dilation values')
    call generic_filter_nd(x15,largest,y15,status,sizes=[3],axes=[15],mode='nearest')
    call require(status==ndimage_status_success,'rank-15 generic callback')
    call require(all(y15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,:)== &
        [4.0_dp,4.0_dp,4.0_dp]),'rank-15 callback values')

    call grey_dilation_nd(x,a,status,sizes=[2,3],axes=[1,1])
    call require(status==ndimage_status_invalid_input,'duplicate grayscale axes')
    call grey_erosion_nd(x,a,status,sizes=[2,3],axes=[0,2])
    call require(status==ndimage_status_invalid_input,'invalid grayscale axis')
    call grey_opening_nd(x,a,status,sizes=[3],axes=[2,3])
    call require(status==ndimage_status_invalid_input,'wrong grayscale sizes length')
    call grey_erosion_nd(x,a,status,sizes=[2,3],axes=[1,3],origin=[2,0])
    call require(status==ndimage_status_invalid_input,'invalid grayscale origin')
    call grey_dilation_nd(x,a,status,sizes=[2,3],axes=[1,3], &
        modes=[character(len=16):: 'wrap'])
    call require(status==ndimage_status_invalid_input,'invalid gray modes length')
    call grey_dilation_nd(x,a,status,sizes=[2,3],axes=[1,3], &
        mode='reflect',modes=[character(len=16):: 'wrap','nearest'])
    call require(status==ndimage_status_invalid_input,'conflicting grayscale mode arguments')
    call generic_filter_nd(x,largest,b,status,sizes=[2,3],axes=[1,1])
    call require(status==ndimage_status_invalid_input,'duplicate generic axes')
    call generic_filter_nd(x,largest,b,status,sizes=[2],axes=[1,3])
    call require(status==ndimage_status_invalid_input,'generic footprint length')
    call generic_filter_nd(x,largest,b,status,footprint=fp,axes=[1])
    call require(status==ndimage_status_invalid_input,'generic footprint rank')
    call generic_filter_nd(x,largest,b,status,footprint=fp,axes=[1,3],origin=[1,0])
    call require(status==ndimage_status_invalid_input,'generic origin invalid')
    call generic_filter_nd(x,largest,b,status,sizes=[2,3],axes=[1,3], &
        modes=[character(len=16):: 'nearest'])
    call require(status==ndimage_status_invalid_input,'generic modes invalid')

    print '(a)', 'PASS: selected-axis grayscale and generic filters'
contains
    subroutine require(ok,message)
        logical, intent(in) :: ok !! test condition
        character(len=*), intent(in) :: message !! assertion identifier
        if (.not.ok) then
            print '(a)', trim(message)
            error stop 'selected-axis gray/generic regression failed'
        end if
    end subroutine require
    function largest(v) result(r)
        real(dp), intent(in) :: v(:) !! active real neighborhood
        real(dp) :: r
        r=maxval(v)
    end function largest
    function largest_complex(v) result(r)
        complex(dp), intent(in) :: v(:) !! active complex neighborhood
        complex(dp) :: r
        r=sum(v)
    end function largest_complex
end program test_ndimage_gray_generic_axes
