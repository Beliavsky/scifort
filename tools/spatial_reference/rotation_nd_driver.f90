! SPDX-License-Identifier: MIT
! Stream reference cases for arbitrary-rank rotation and SE(3) broadcasting.
program rotation_nd_driver
    use scifort_kinds, only : dp
    use scifort_spatial_transform_nd, only : rotation_nd_type, rigid_nd_type, &
        rotation_nd_from_quat, rigid_nd_from_components
    implicit none
    integer :: na(4), nb(4), nc(4), i, trial, ncases, code
    integer, allocatable :: shape_out(:)
    real(dp), allocatable :: qa(:,:,:,:,:), qb(:,:,:,:,:)
    real(dp), allocatable :: ta(:,:,:,:,:), tb(:,:,:,:,:), points(:,:,:,:,:)
    real(dp), allocatable :: rmat(:,:,:), tmat(:,:,:), rpoints(:,:), tpoints(:,:)
    type(rotation_nd_type) :: ra, rb, rc
    type(rigid_nd_type) :: aa, bb, cc
    read(*,*) ncases
    do trial=1,ncases
        read(*,*) na
        read(*,*) nb
        allocate(qa(na(1),na(2),na(3),na(4),4))
        allocate(qb(nb(1),nb(2),nb(3),nb(4),4))
        allocate(ta(na(1),na(2),na(3),na(4),3))
        allocate(tb(nb(1),nb(2),nb(3),nb(4),3))
        allocate(points(nb(1),nb(2),nb(3),nb(4),3))
        read(*,*) qa
        read(*,*) qb
        read(*,*) ta
        read(*,*) tb
        read(*,*) points
        ra=rotation_nd_from_quat(qa,status=code)
        if (code/=0) error stop 'A rotation constructor'
        rb=rotation_nd_from_quat(qb,status=code)
        if (code/=0) error stop 'B rotation constructor'
        rc=ra%compose(rb,status=code)
        if (code/=0) error stop 'rotation broadcasting'
        aa=rigid_nd_from_components(ta,ra,status=code)
        if (code/=0) error stop 'A rigid constructor'
        bb=rigid_nd_from_components(tb,rb,status=code)
        if (code/=0) error stop 'B rigid constructor'
        cc=aa%compose(bb,status=code)
        if (code/=0) error stop 'rigid broadcasting'
        rmat=rc%as_matrix()
        tmat=cc%as_matrix()
        rpoints=rc%apply(points,shape_out,code)
        if (code/=0) error stop 'rotation application'
        tpoints=cc%apply(points,status=code)
        if (code/=0) error stop 'rigid application'
        nc=rc%shape()
        if (any(shape_out/=nc)) error stop 'shape disagreement'
        write(*,'(A,1X,4(I0,1X))') 'CASE',nc
        do i=1,size(rmat,1)
            write(*,'(*(ES26.17E3,1X))') reshape(rmat(i,:,:),[9]), &
                reshape(tmat(i,:,:),[16]), rpoints(i,:), tpoints(i,:)
        end do
        deallocate(qa,qb,ta,tb,points)
    end do
end program rotation_nd_driver
