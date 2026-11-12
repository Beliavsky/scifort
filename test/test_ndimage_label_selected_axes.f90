! SPDX-License-Identifier: MIT
! Connectivity restricted to selected axes (a SciFort extension).
program test_ndimage_label_selected_axes
    use scifort_kinds, only : dp
    use scifort_ndimage, only : label_nd
    implicit none
    logical :: mask(6,4,7),mask15(3,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    integer :: res(6,4,7),compare(6,4,7),plain(6,7),out15(3,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    integer :: numeric(6,4,7),c,j,k,i,u,v,total_features,features,status
    real(dp) :: rmask(6,4,7)
    logical :: struct(3,3)
    integer :: many

    do k=1,7
        do j=1,4
            do i=1,6
                mask(i,j,k)=(modulo(2*i+3*j+7*k,9)<4)
            end do
        end do
    end do
    struct=reshape([.true.,.false.,.true.,.false.,.true.,.false., &
        .true.,.false.,.true.],[3,3])
    call label_nd(mask,res,features,status,axes=[3,1],structure=struct)
    call check(status==0,'rank-3 selected labeling')
    total_features=0
    do j=1,4
        call label_nd(mask(:,j,:),plain,c,status,structure=struct)
        call check(status==0,'sliced labeling')
        total_features=total_features+c
        do k=1,7
            do i=1,6
                if (.not.mask(i,j,k)) then
                    call check(res(i,j,k)==0,'background remains zero')
                    cycle
                end if
                do v=1,7
                    do u=1,6
                        if (.not.mask(u,j,v)) cycle
                        call check((res(i,j,k)==res(u,j,v)) .eqv. &
                            (plain(i,k)==plain(u,v)),'same connectivity within slice')
                    end do
                end do
            end do
        end do
    end do
    call check(features==total_features,'component count across independent slices')
    do j=1,3
        do k=j+1,4
            do u=1,6
                do v=1,7
                    if (.not.mask(u,j,v)) cycle
                    call check(.not.any(res(:,k,:)==res(u,j,v)), &
                        'global label IDs distinct across independent slices')
                end do
            end do
        end do
    end do
    numeric=merge(1,0,mask)
    rmask=real(numeric,dp)
    call label_nd(numeric,compare,many,status,axes=[1,3],structure=struct)
    call check(status==0.and.many==features,'numeric axis labeling')
    call check(all(compare==res),'numeric and logical partitions')
    call label_nd(rmask,compare,many,status,axes=[1,3],structure=struct)
    call check(status==0.and.many==features,'real axis labeling')
    call check(all(compare==res),'real and logical partitions')

    call label_nd(mask,compare,many,status,axes=[integer::])
    call check(status==0.and.many==count(mask),'zero selected axes disconnect pixels')
    call label_nd(mask,compare,many,status,axes=[1,1])
    call check(status/=0,'reject repeated axes')
    call label_nd(mask,compare,many,status,axes=[4])
    call check(status/=0,'reject axes exceeding rank')
    call label_nd(mask,compare,many,status,axes=[1,3],structure=spread(.true.,1,3))
    call check(status/=0,'reject footprint with wrong compact rank')

    mask15=.false.
    mask15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=.true.
    mask15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=.true.
    call label_nd(mask15,out15,many,status,axes=[1])
    call check(status==0.and.many==1,'rank-15 single-axis connected component')
    print '(a)', 'PASS: selected-axis connected component labeling'
contains
    subroutine check(ok,message)
        logical, intent(in) :: ok !! assertion result
        character(len=*), intent(in) :: message !! assertion description
        if (.not.ok) then
            print '(2a)', 'FAIL: ',message
            error stop 1
        end if
    end subroutine check
end program test_ndimage_label_selected_axes
