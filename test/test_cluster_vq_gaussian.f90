! SPDX-License-Identifier: MIT
program test_cluster_vq_gaussian
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_random, only : rng_state, rng_seed, rng_uniform
    use scifort_cluster_vq, only : kmeans2, cluster_vq_success, &
        cluster_vq_empty_cluster, cluster_vq_invalid_input
    implicit none

    real(dp) :: x(9,3), wide(3,6), constant(5,3), one(1,3), &
        initial(3,3), mu(3), normal(9), two_pi, u1, u2
    real(dp), allocatable :: centers(:,:), centers2(:,:), vector_centers(:)
    integer, allocatable :: labels(:), labels2(:)
    type(rng_state) :: state, state2
    integer :: i, j, status, status2

    do i = 1,9
        x(i,1) = real(mod(7*i,11)-5,dp)
        x(i,2) = real(i*i-33,dp)*0.11_dp
        x(i,3) = x(i,1)+2.0_dp*x(i,2)  ! Exactly rank deficient.
    end do
    call rng_seed(state,91245)
    call rng_seed(state2,91245)
    two_pi = 2.0_dp*acos(-1.0_dp)
    mu = sum(x,dim=1)/9.0_dp
    do i = 1,3
        do j = 1,9
            u1 = rng_uniform(state2)
            u2 = rng_uniform(state2)
            normal(j) = sqrt(-2.0_dp*log(u1))*cos(two_pi*u2)
        end do
        ! Independently compute N(mu,sample covariance) in the original
        ! input scale and use that initial codebook for deterministic Lloyd.
        initial(i,:) = mu + matmul(normal,x-spread(mu,1,9))/sqrt(8.0_dp)
    end do
    call kmeans2(x,3,centers,labels,status,state=state,init='random',iterations=4)
    call kmeans2(x,initial,centers2,labels2,status2,iterations=4)
    if (status /= status2 .or. &
        (status /= cluster_vq_success .and. status /= cluster_vq_empty_cluster)) &
        error stop 'Gaussian initialization status'
    if (any(shape(centers) /= shape(centers2))) &
        error stop 'Gaussian initialization center dimensions'
    if (maxval(abs(centers-centers2)) > 3.0e-13_dp .or. &
        any(labels /= labels2)) error stop 'Gaussian init vs independent sample'

    call rng_seed(state,91245)
    call kmeans2(x,3,centers2,labels2,status2,state=state,init='random',iterations=4)
    if (status2 /= status .or. any(centers2 /= centers) .or. any(labels2 /= labels)) &
        error stop 'Gaussian initialization repeatable'

    ! More features than observations: no singular-covariance breakdown.
    wide(1,:) = [1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp]
    wide(2,:) = [2.0_dp,4.0_dp,6.0_dp,8.0_dp,10.0_dp,12.0_dp]
    wide(3,:) = [3.0_dp,6.0_dp,9.0_dp,12.0_dp,15.0_dp,18.0_dp]
    call rng_seed(state,123456)
    call kmeans2(wide,2,centers,labels,status,state=state,init='random',iterations=4)
    if ((status /= cluster_vq_success .and. status /= cluster_vq_empty_cluster) .or. &
        any(shape(centers) /= [2,6]) .or. .not. all(ieee_is_finite(centers))) &
        error stop 'Gaussian rank-deficient wide observations'

    constant = 23.0_dp
    call rng_seed(state,7)
    call kmeans2(constant,1,centers,labels,status,state=state,init='random')
    if (status /= cluster_vq_success .or. any(centers /= 23.0_dp) .or. &
        any(labels /= 1)) error stop 'Gaussian constant observations'

    call rng_seed(state,91245)
    call kmeans2(x(:,1),2,vector_centers,labels,status, &
        state=state,init='random',iterations=3)
    if ((status /= cluster_vq_success .and. status /= cluster_vq_empty_cluster) .or. &
        size(vector_centers) /= 2) error stop 'Gaussian 1d observations'

    one = 1.0_dp
    call kmeans2(one,1,centers,labels,status,state=state,init='random')
    if (status /= cluster_vq_invalid_input .or. size(centers) /= 0) &
        error stop 'Gaussian n=1 needs sample variance'
    call kmeans2(x,3,centers,labels,status,init='random')
    if (status /= cluster_vq_invalid_input) error stop 'Gaussian explicit RNG required'
    call kmeans2(x,3,centers,labels,status,state=state,init='unknown')
    if (status /= cluster_vq_invalid_input) error stop 'Gaussian unknown init option'

    print '(a)', 'cluster vq Gaussian initialization regression PASS'
end program test_cluster_vq_gaussian
