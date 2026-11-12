! SPDX-License-Identifier: MIT
program linalg_subspaces_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_subspaces, only : orth, null_space, subspace_angles
    implicit none
    real(dp) :: a(4,5), b(4,1)
    real(dp), allocatable :: range_basis(:, :), kernel(:, :), angles(:)
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer :: status

    ! A has rank two, with three independent null directions.
    a = 0.0_dp
    a(1,1) = 1.0_dp
    a(2,2) = 1.0_dp
    a(1,3) = 1.0_dp
    call orth(a,range_basis,status)
    if (status /= linalg_status_success) error stop 'orth failed'
    call null_space(a,kernel,status)
    if (status /= linalg_status_success) error stop 'null_space failed'
    print '(A,I0)', 'dimension of column space: ',size(range_basis,2)
    print '(A,I0)', 'dimension of null space: ',size(kernel,2)
    print '(A,ES12.4)', 'maximum null-space residual: ',maxval(abs(matmul(a,kernel)))

    ! Angle between span{e1,e2} and a vector tilted from e1 by 30 degrees.
    b = 0.0_dp
    b(1,1) = cos(pi/6.0_dp)
    b(3,1) = sin(pi/6.0_dp)
    call subspace_angles(a,b,angles,status)
    if (status /= linalg_status_success) error stop 'subspace_angles failed'
    print '(A,F12.6)', 'principal angle (degrees): ',angles(1)*180.0_dp/pi
end program linalg_subspaces_example
