! SPDX-License-Identifier: MIT
! Reproducible binary I/O for randomized directed chamfer reference tests.
program weighted_stencil_random_driver
    use scifort_ndimage, only : distance_transform_weighted_cdt_nd
    implicit none
    integer :: iu,ou,case_id,status
    integer :: raw2(4,5),w2(3,3),d2(4,5),features2(2,20)
    integer :: raw3(3,3,2),w3(3,3,3),d3(3,3,2),features3(3,18)

    open(newunit=iu,file='weighted_stencil_random_inputs.bin',access='stream', &
        form='unformatted',status='old')
    open(newunit=ou,file='weighted_stencil_random_results.bin',access='stream', &
        form='unformatted',status='replace')
    do case_id=1,96
        read(iu) raw2,w2
        call distance_transform_weighted_cdt_nd(raw2/=0,w2,d2,status,features2)
        if (status/=0) error stop 'weighted random 2D status'
        write(ou) d2,features2
        read(iu) raw3,w3
        call distance_transform_weighted_cdt_nd(raw3/=0,w3,d3,status,features3)
        if (status/=0) error stop 'weighted random 3D status'
        write(ou) d3,features3
    end do
    close(iu)
    close(ou)
    print '(a)', 'PASS: randomized 2D and 3D weighted chamfer Fortran driver'
end program weighted_stencil_random_driver
