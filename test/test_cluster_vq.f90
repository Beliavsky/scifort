! SPDX-License-Identifier: MIT
program test_cluster_vq
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_random, only : rng_state, rng_seed
    use scifort_cluster_vq, only : whiten, vq, kmeans, kmeans2, &
        cluster_vq_success, cluster_vq_invalid_input, cluster_vq_empty_cluster
    implicit none

    real(dp) :: x(6,2), initial(2,2), bad(3,2), flat(3,3), expected(2,2)
    real(dp) :: zeros(3,2), duplicate_centers(2,2)
    real(dp), allocatable :: normalized(:,:), vector(:), vector_centers(:), &
        centers(:,:), centers2(:,:), distances(:)
    integer, allocatable :: labels(:), labels2(:)
    type(rng_state) :: state, state2
    real(dp) :: distortion
    integer :: status

    x = transpose(reshape([0.0_dp,0.0_dp, 0.0_dp,1.0_dp, 1.0_dp,0.0_dp, &
        10.0_dp,10.0_dp, 10.0_dp,11.0_dp, 11.0_dp,10.0_dp],[2,6]))
    initial = transpose(reshape([0.0_dp,0.0_dp, 10.0_dp,10.0_dp],[2,2]))
    expected = transpose(reshape([1.0_dp/3.0_dp,1.0_dp/3.0_dp, &
        31.0_dp/3.0_dp,31.0_dp/3.0_dp],[2,2]))

    flat = transpose(reshape([1.0_dp,10.0_dp,9.0_dp, 3.0_dp,20.0_dp,9.0_dp, &
        5.0_dp,30.0_dp,9.0_dp],[3,3]))
    normalized = whiten(flat,status)
    if (status/=cluster_vq_success) error stop 'whiten matrix status'
    if (any(shape(normalized)/=[3,3])) error stop 'whiten matrix shape'
    if (maxval(abs(normalized(:,1)-flat(:,1)/sqrt(8.0_dp/3.0_dp)))>1.0e-13_dp) &
        error stop 'whiten first column'
    if (maxval(abs(normalized(:,2)-flat(:,2)/sqrt(200.0_dp/3.0_dp)))>1.0e-13_dp) &
        error stop 'whiten second column'
    if (any(normalized(:,3)/=9.0_dp)) error stop 'whiten constant column'
    vector = whiten([1.0_dp,3.0_dp,5.0_dp],status)
    if (status/=cluster_vq_success .or. size(vector)/=3) error stop 'whiten vector'
    if (abs(vector(3)-5.0_dp/sqrt(8.0_dp/3.0_dp))>1.0e-13_dp) &
        error stop 'whiten vector value'
    vector = whiten([7.0_dp,7.0_dp],status)
    if (any(vector/=7.0_dp)) error stop 'whiten vector constant'
    normalized = whiten(reshape([1.0e200_dp,2.0e200_dp,3.0e200_dp],[3,1]),status)
    if (status/=cluster_vq_success) error stop 'large whiten status'
    if (abs(normalized(3,1)-3.0_dp/sqrt(2.0_dp/3.0_dp))>1.0e-12_dp) &
        error stop 'large whiten values'

    call vq(x,initial,labels,distances,status)
    if (status/=cluster_vq_success) error stop 'vq status'
    if (any(labels/=[1,1,1,2,2,2])) error stop 'vq indices'
    if (any(abs(distances-[0.0_dp,1.0_dp,1.0_dp,0.0_dp,1.0_dp,1.0_dp]) &
        > 1.0e-14_dp)) error stop 'vq distances'
    call vq([0.0_dp,1.0_dp,6.0_dp,9.0_dp],[0.0_dp,10.0_dp], &
        labels,distances,status)
    if (status/=cluster_vq_success .or. any(labels/=[1,1,2,2])) &
        error stop 'one-dimensional vq indices'
    if (any(abs(distances-[0.0_dp,1.0_dp,4.0_dp,1.0_dp])>1.0e-13_dp)) error stop '1d vq distances'
    call vq([2.0_dp],[2.0_dp,2.0_dp],labels,distances,status)
    if (labels(1)/=1 .or. distances(1)/=0.0_dp) error stop 'vq tie breaking'
    call vq([1.0e200_dp],[0.0_dp,1.0e200_dp],labels,distances,status)
    if (status/=cluster_vq_success .or. labels(1)/=2 .or. distances(1)/=0.0_dp) &
        error stop 'large vq distance'

    call kmeans2(x,initial,centers,labels,status,iterations=10)
    if (status/=cluster_vq_success) error stop 'kmeans2 matrix status'
    if (maxval(abs(centers-expected))>1.0e-13_dp) error stop 'kmeans2 centers'
    if (any(labels/=[1,1,1,2,2,2])) error stop 'kmeans2 labels'
    call kmeans(x,initial,centers,distortion,status)
    if (status/=cluster_vq_success) error stop 'kmeans matrix status'
    if (maxval(abs(centers-expected))>1.0e-13_dp) error stop 'kmeans centers'
    if (abs(distortion-0.6540388352636305_dp)>1.0e-13_dp) &
        error stop 'mean non-squared distortion'

    call kmeans2([0.0_dp,1.0_dp,2.0_dp,10.0_dp,11.0_dp,12.0_dp], &
        [0.0_dp,10.0_dp],vector_centers,labels,status,iterations=10)
    if (status/=cluster_vq_success) error stop 'kmeans2 1d status'
    if (maxval(abs(vector_centers-[1.0_dp,11.0_dp]))>1.0e-13_dp) &
        error stop 'kmeans2 1d centroids'
    if (any(labels/=[1,1,1,2,2,2])) error stop 'kmeans2 1d labels'
    call kmeans([0.0_dp,1.0_dp,2.0_dp,10.0_dp,11.0_dp,12.0_dp], &
        [0.0_dp,10.0_dp],vector_centers,distortion,status)
    if (status/=cluster_vq_success) error stop 'kmeans 1d status'
    if (maxval(abs(vector_centers-[1.0_dp,11.0_dp]))>1.0e-13_dp) &
        error stop 'kmeans 1d centroids'
    if (abs(distortion-2.0_dp/3.0_dp)>1.0e-13_dp) &
        error stop 'kmeans 1d mean distortion'

    call rng_seed(state,20261103)
    call kmeans2([0.0_dp,1.0_dp,2.0_dp,10.0_dp,11.0_dp,12.0_dp], &
        2,vector_centers,labels,status,state=state,iterations=10,init='++')
    if (status/=cluster_vq_success) error stop 'kmeans2 1d plus plus'
    call rng_seed(state,20261103)
    call kmeans([0.0_dp,1.0_dp,2.0_dp,10.0_dp,11.0_dp,12.0_dp], &
        2,vector_centers,distortion,status,state=state,restarts=8)
    if (status/=cluster_vq_success) error stop 'kmeans 1d restarts'
    if (abs(distortion-2.0_dp/3.0_dp)>1.0e-13_dp) &
        error stop 'kmeans 1d restart distortion'

    call rng_seed(state,20261103)
    call rng_seed(state2,20261103)
    call kmeans2(x,2,centers,labels,status,state=state,init='++',iterations=8)
    if (status/=cluster_vq_success) error stop 'kmeans2 plus plus status'
    call kmeans2(x,2,centers2,labels2,status,state=state2,init='++',iterations=8)
    if (status/=cluster_vq_success) error stop 'kmeans2 repeat status'
    if (any(centers/=centers2) .or. any(labels/=labels2)) &
        error stop 'kmeans2 reproducibility'
    if (any(abs(centers(:,1)-[1.0_dp/3.0_dp,31.0_dp/3.0_dp])>1.0e-13_dp) .and. &
        any(abs(centers(:,1)-[31.0_dp/3.0_dp,1.0_dp/3.0_dp])>1.0e-13_dp)) &
        error stop 'kmeans2 plus plus partition'
    call kmeans2(x,2,centers,labels,status,state=state,init='points',iterations=8)
    if (status/=cluster_vq_success) error stop 'kmeans2 points status'

    call rng_seed(state,9431)
    call kmeans(x,2,centers,distortion,status,state=state,restarts=12)
    if (status/=cluster_vq_success .or. size(centers,1)/=2) &
        error stop 'kmeans restarts status'
    if (abs(distortion-0.6540388352636305_dp)>1.0e-12_dp) &
        error stop 'kmeans restarts distortion'

    zeros = 0.0_dp
    duplicate_centers = 0.0_dp
    call kmeans2(zeros,duplicate_centers,centers,labels,status,iterations=4)
    if (status/=cluster_vq_empty_cluster) error stop 'kmeans2 empty warning'
    if (any(shape(centers)/=[2,2]) .or. any(labels/=1)) &
        error stop 'kmeans2 empty warning results'
    call kmeans2(zeros,duplicate_centers,centers,labels,status,missing='raise')
    if (status/=cluster_vq_empty_cluster) error stop 'kmeans2 empty raise status'
    if (size(centers)/=0 .or. size(labels)/=0) error stop 'kmeans2 empty raise outputs'
    call kmeans(zeros,duplicate_centers,centers,distortion,status)
    if (status/=cluster_vq_success .or. size(centers,1)/=1) &
        error stop 'kmeans removes empty clusters'
    if (distortion/=0.0_dp) error stop 'kmeans empty distortion'

    bad = 0.0_dp
    bad(2,2) = ieee_value(0.0_dp,ieee_quiet_nan)
    normalized = whiten(bad,status)
    if (status/=cluster_vq_invalid_input .or. size(normalized)/=0) &
        error stop 'whiten invalid NaN'
    call vq(bad,initial,labels,distances,status)
    if (status/=cluster_vq_invalid_input .or. size(labels)/=0) error stop 'vq invalid NaN'
    call kmeans2(bad,initial,centers,labels,status)
    if (status/=cluster_vq_invalid_input .or. size(centers)/=0) &
        error stop 'kmeans2 invalid NaN'
    call kmeans2(x,0,centers,labels,status,state=state)
    if (status/=cluster_vq_invalid_input) error stop 'kmeans2 invalid k'
    call kmeans2(x,2,centers,labels,status,init='++')
    if (status/=cluster_vq_invalid_input) error stop 'kmeans2 requires rng'
    call kmeans2(x,2,centers,labels,status,state=state,iterations=0)
    if (status/=cluster_vq_invalid_input) error stop 'kmeans2 iterations=0'
    call kmeans2(x,initial,centers,labels,status,missing='bogus')
    if (status/=cluster_vq_invalid_input) error stop 'kmeans2 missing option'
    call kmeans(x,initial,centers,distortion,status,threshold=-1.0_dp)
    if (status/=cluster_vq_invalid_input) error stop 'kmeans invalid threshold'
    call kmeans(x,initial,centers,distortion,status,max_iterations=0)
    if (status/=cluster_vq_invalid_input) error stop 'kmeans invalid max_iter'

    print '(a)', 'cluster vq regression PASS'
end program test_cluster_vq
