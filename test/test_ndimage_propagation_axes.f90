! SPDX-License-Identifier: MIT
! Axis-selective binary propagation and hole filling regression.
program test_ndimage_propagation_axes
    use scifort_ndimage, only : binary_propagation_nd, binary_fill_holes_nd
    implicit none
    logical :: seeds(7,4,7),mask(7,4,7),res(7,4,7),expected(7,4,7)
    logical :: image(7,4,7),foot(3,3),out_big(14,4,7)
    logical :: x15(3,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    logical :: y15(3,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    integer :: i,j,k,status

    do k=1,7
        do j=1,4
            do i=1,7
                seeds(i,j,k)=modulo(7*i+11*j+3*k,23)==0
                mask(i,j,k)=modulo(i+3*j+5*k,13)/=0
            end do
        end do
    end do
    foot=reshape([.false.,.true.,.false.,.true.,.true.,.true., &
        .false.,.true.,.false.],shape(foot))
    call binary_propagation_nd(seeds,mask,res,status,axes=[3,1],structure=foot)
    call check(status==0,'axis-selective propagation status')
    do j=1,4
        call binary_propagation_nd(seeds(:,j,:),mask(:,j,:),expected(:,j,:),status,structure=foot)
        call check(status==0,'slice propagation status')
    end do
    call check(all(res.eqv.expected),'axis propagation matches sliced reference')
    call binary_propagation_nd(seeds,mask,res,status,axes=[integer::])
    call check(status==0.and.all(res.eqv.seeds),'empty axes propagation identity')
    call binary_propagation_nd(seeds,mask,res,status,axes=[1,1])
    call check(status/=0,'reject duplicate propagation axes')
    call binary_propagation_nd(seeds,mask,res,status,axes=[1,3],structure=spread(.true.,1,3))
    call check(status/=0,'reject invalid compact propagation footprint')
    call binary_propagation_nd(seeds,mask,res,status,axes=[4])
    call check(status/=0,'reject invalid propagation axis')

    image=.true.
    do j=1,4
        image(2:6,j,2:6)=.false.
        image(4,j,4)=.true.
        image(1,j,3)=.false.
    end do
    call binary_fill_holes_nd(image,res,status,axes=[1,3],structure=foot)
    call check(status==0,'selected-axis fill holes status')
    do j=1,4
        call binary_fill_holes_nd(image(:,j,:),expected(:,j,:),status,structure=foot)
        call check(status==0,'slice fill holes status')
    end do
    call check(all(res.eqv.expected),'selected-axis holes match sliced reference')
    out_big=.false.
    call binary_fill_holes_nd(image,out_big(1:14:2,:,:),status,axes=[1,3])
    call check(status==0,'noncontiguous fill result')
    call binary_fill_holes_nd(image,expected,status,axes=[1,3])
    call check(all(out_big(1:14:2,:,:).eqv.expected),'noncontiguous fill values')
    call binary_fill_holes_nd(image,res,status,axes=[integer::])
    call check(status==0.and.all(res.eqv.image),'empty axes fill holes identity')
    call binary_fill_holes_nd(image,res,status,axes=[2,2])
    call check(status/=0,'reject duplicate fill axes')
    call binary_fill_holes_nd(image,res,status,axes=[2,3],structure=spread(.true.,1,3))
    call check(status/=0,'reject invalid fill footprint')

    x15=.true.
    x15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=.false.
    call binary_fill_holes_nd(x15,y15,status,axes=[1])
    call check(status==0.and.all(y15),'rank15 selected-axis hole filling')
    call binary_propagation_nd(x15,x15,y15,status,axes=[1])
    call check(status==0.and.all(y15.eqv.x15),'rank15 selected propagation')
    print '(a)', 'PASS: selected-axis binary propagation and hole filling'
contains
    subroutine check(condition,message)
        logical, intent(in) :: condition !! assertion
        character(len=*), intent(in) :: message !! assertion meaning
        if (.not.condition) then
            print '(2a)', 'FAIL: ',message
            error stop 1
        end if
    end subroutine check
end program test_ndimage_propagation_axes
