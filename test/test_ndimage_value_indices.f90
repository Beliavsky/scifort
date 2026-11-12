! SPDX-License-Identifier: MIT
! Compressed index lookup for integer region maps.
program test_ndimage_value_indices
    use scifort_ndimage, only : value_indices_nd, ndimage_status_invalid_input
    implicit none
    integer :: a(3,2),status
    integer :: high(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    integer, allocatable :: keys(:),starts(:),coords(:,:)
    integer :: negative(2),empty(0)

    a=reshape([2,0,3,2,0,2],shape(a))
    call value_indices_nd(a,keys,starts,coords,status)
    if (status/=0) error stop 'value indices status'
    if (any(keys/=[0,2,3])) error stop 'value indices keys'
    if (any(starts/=[1,3,6,7])) error stop 'value indices starts'
    if (any(coords(:,starts(2):starts(3)-1)/= &
            reshape([1,1,1,2,3,2],[2,3]))) then
        ! Label 2 has C-order coordinates (1,1), (1,2), (3,2).
        error stop 'value indices grouped coordinates'
    end if
    call value_indices_nd(a,keys,starts,coords,status,ignore_value=0)
    if (status/=0 .or. any(keys/=[2,3])) error stop 'ignore background keys'
    if (any(starts/=[1,4,5])) error stop 'ignore background starts'
    call value_indices_nd(a,keys,starts,coords,status,ignore_value=11)
    if (status/=0 .or. size(keys)/=3) error stop 'ignore absent key'
    negative=[-huge(1),huge(1)]
    call value_indices_nd(negative,keys,starts,coords,status)
    if (status/=0 .or. any(keys/=negative)) error stop 'extreme signed keys'
    high=7
    call value_indices_nd(high,keys,starts,coords,status,ignore_value=7)
    if (status/=0) error stop 'all ignored status'
    if (size(keys)/=0 .or. any(starts/=[1])) error stop 'all ignored starts'
    if (size(coords,1)/=15 .or. size(coords,2)/=0) error stop 'all ignored coords'
    high(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=8
    call value_indices_nd(high,keys,starts,coords,status,ignore_value=7)
    if (status/=0 .or. any(keys/=[8]) .or. any(starts/=[1,2])) &
        error stop 'rank15 unique key'
    if (any(coords(:,1)/=1)) error stop 'rank15 coordinates'
    call value_indices_nd(empty,keys,starts,coords,status)
    if (status/=ndimage_status_invalid_input) error stop 'empty input invalid'
    print '(a)', 'PASS: value indices grouped coordinates and rank15 cases'
end program test_ndimage_value_indices
