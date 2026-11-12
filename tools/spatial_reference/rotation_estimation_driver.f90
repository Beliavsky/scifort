! SPDX-License-Identifier: MIT
! Standard input/output driver for independent scipy.spatial.transform comparisons.
program rotation_estimation_driver
    use scifort_kinds, only : dp
    use scifort_spatial_transform, only : rotation_type, rotation_from_quat, &
        rotation_align_vectors, rotation_mean
    implicit none
    type(rotation_type) :: fit, mean_rotation
    type(rotation_type), allocatable :: rotations(:)
    real(dp), allocatable :: a(:,:), b(:,:), w(:), wm(:)
    real(dp) :: q(4), rssd, sensitivity(3,3)
    integer :: total, n, m, i, j, test, want_sens, status_fit, status_mean

    read(*,*) total
    do test=1,total
        read(*,*) n, m, want_sens
        allocate(a(n,3), b(n,3), w(n), rotations(m), wm(m))
        do i=1,n
            read(*,*) a(i,1:3), b(i,1:3), w(i)
        end do
        do i=1,m
            read(*,*) q, wm(i)
            rotations(i) = rotation_from_quat(q)
        end do
        if (want_sens == 1) then
            fit = rotation_align_vectors(a,b,weights=w,rssd=rssd, &
                                         sensitivity=sensitivity,status=status_fit)
        else
            fit = rotation_align_vectors(a,b,weights=w,rssd=rssd,status=status_fit)
            sensitivity = 0.0_dp
        end if
        mean_rotation = rotation_mean(rotations,weights=wm,status=status_mean)
        write(*,'(2(i0,1x))') status_fit, status_mean
        write(*,'(18(es25.16e3,1x))') fit%as_quat(),rssd,mean_rotation%as_quat(), &
                                       ((sensitivity(i,j),j=1,3),i=1,3)
        deallocate(a,b,w,rotations,wm)
    end do
end program rotation_estimation_driver
