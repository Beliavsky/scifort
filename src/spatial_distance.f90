! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Pairwise distances and condensed/square distance matrix conversion.
! Independent implementations of the standard metric definitions and
! the condensed indexing convention documented by scipy.spatial.distance.
module scifort_spatial_distance
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_quiet_nan, &
        ieee_positive_inf, ieee_value
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: spatial_distance_success = 0
    integer, parameter, public :: spatial_distance_invalid_input = 1

    public :: point_distance, pdist, cdist, squareform

    interface squareform
        module procedure squareform_to_matrix
        module procedure squareform_to_vector
    end interface squareform

contains

    pure function invalid_value() result(v)
        real(dp) :: v
        v = ieee_value(0.0_dp, ieee_quiet_nan)
    end function invalid_value

    pure logical function valid_metric(metric) result(ok)
        character(len=*), intent(in) :: metric !! SciPy metric name
        select case (metric)
        case ('euclidean', 'sqeuclidean', 'cityblock', 'chebyshev', 'minkowski', &
            'cosine', 'correlation', 'hamming', 'jaccard', 'canberra', 'braycurtis', &
            'seuclidean', 'mahalanobis', 'jensenshannon', &
            'dice', 'rogerstanimoto', 'russellrao', 'sokalsneath', 'yule')
            ok = .true.
        case default
            ok = .false.
        end select
    end function valid_metric

    pure logical function valid_options(n, metric, p, w, variance, inv_cov, base) result(ok)
        integer, intent(in) :: n !! positive number of components
        character(len=*), intent(in) :: metric !! supported metric
        real(dp), intent(in), optional :: p !! Minkowski order, >0
        real(dp), intent(in), optional :: w(:) !! optional coordinate weights
        real(dp), intent(in), optional :: variance(:) !! strictly positive component variances
        real(dp), intent(in), optional :: inv_cov(:,:) !! inverse covariance matrix
        real(dp), intent(in), optional :: base !! Jensen-Shannon logarithm base >1
        ok = .false.
        if (n <= 0 .or. .not. valid_metric(metric)) return
        if (present(p)) then
            if (metric /= 'minkowski') return
            if (.not. ieee_is_finite(p)) return
            if (p <= 0.0_dp) return
        end if
        if (present(w)) then
            if (metric=='seuclidean' .or. metric=='mahalanobis' .or. metric=='jensenshannon') return
            if (size(w) /= n) return
            if (.not. all(ieee_is_finite(w))) return
            if (any(w < 0.0_dp)) return
        end if
        if (present(variance)) then
            if (metric /= 'seuclidean') return
            if (size(variance) /= n) return
            if (.not. all(ieee_is_finite(variance))) return
            if (any(variance <= 0.0_dp)) return
        end if
        if (present(inv_cov)) then
            if (metric /= 'mahalanobis') return
            if (size(inv_cov,1) /= n .or. size(inv_cov,2) /= n) return
            if (.not. all(ieee_is_finite(inv_cov))) return
        end if
        if (present(base)) then
            if (metric /= 'jensenshannon') return
            if (.not. ieee_is_finite(base)) return
            if (base <= 1.0_dp) return
        end if
        ok = .true.
    end function valid_options

    ! Sample variance uses ddof=1, matching SciPy's defaults for pdist/cdist.
    ! y (if present) supplies a second batch of observations for cdist.
    pure subroutine estimate_variance(x, variances, ok, y)
        real(dp), intent(in) :: x(:,:) !! first observations (rows, features)
        real(dp), intent(out) :: variances(:) !! estimated component variances
        logical, intent(out) :: ok !! true if every variance is positive and finite
        real(dp), intent(in), optional :: y(:,:) !! second batch of observations
        real(dp) :: averages(size(x,2)), total
        integer :: n, j
        ok = .false.
        n = size(x,1)
        if (present(y)) n = n + size(y,1)
        if (n < 2) return
        total = real(n,dp)
        averages = sum(x,dim=1)/total
        if (present(y)) averages = averages + sum(y,dim=1)/total
        if (.not. all(ieee_is_finite(averages))) return
        do j=1,size(x,2)
            variances(j) = sum((x(:,j)-averages(j))**2)
            if (present(y)) variances(j) = variances(j) + sum((y(:,j)-averages(j))**2)
        end do
        variances = variances/real(n-1,dp)
        ok = all(ieee_is_finite(variances)) .and. all(variances>0.0_dp)
    end subroutine estimate_variance

    ! Estimate the sample covariance and invert its Cholesky factor without
    ! LAPACK. Rank deficient data yield invalid_input instead of NaNs/infs.
    pure subroutine estimate_inverse_covariance(x, inverse, ok, y)
        real(dp), intent(in) :: x(:,:) !! first observations (rows, features)
        real(dp), intent(out) :: inverse(:,:) !! inverse sample covariance
        logical, intent(out) :: ok !! false for singular covariance/insufficient data
        real(dp), intent(in), optional :: y(:,:) !! additional observations
        real(dp) :: averages(size(x,2)), centered(size(x,2))
        real(dp) :: cov(size(x,2),size(x,2)), lower(size(x,2),size(x,2))
        real(dp) :: lower_inv(size(x,2),size(x,2)), stddev(size(x,2)), s
        integer :: n, ndim, i, j, k, batch, obs
        ok = .false.
        ndim = size(x,2)
        n = size(x,1)
        if (present(y)) n = n + size(y,1)
        if (n <= ndim) return
        averages = sum(x,dim=1)/real(n,dp)
        if (present(y)) averages = averages + sum(y,dim=1)/real(n,dp)
        if (.not. all(ieee_is_finite(averages))) return
        cov = 0.0_dp
        do batch=1,2
            if (batch==2 .and. .not. present(y)) exit
            obs = size(x,1)
            if (batch==2) obs = size(y,1)
            do k=1,obs
                if (batch==1) then
                    centered = x(k,:)-averages
                else
                    centered = y(k,:)-averages
                end if
                do j=1,ndim
                    do i=1,ndim
                        cov(i,j) = cov(i,j) + centered(i)*centered(j)
                    end do
                end do
            end do
        end do
        cov = cov/real(n-1,dp)
        if (.not. all(ieee_is_finite(cov))) return
        ! Normalize to a correlation matrix before factorization. Otherwise
        ! unequal feature units can trigger a false singularity (for example
        ! variances of 1e+16 and 1e-16 in the same observation matrix).
        do i=1,ndim
            if (cov(i,i) <= 0.0_dp) return
            stddev(i) = sqrt(cov(i,i))
        end do
        do j=1,ndim
            do i=1,ndim
                cov(i,j) = (cov(i,j)/stddev(i))/stddev(j)
            end do
        end do
        if (.not. all(ieee_is_finite(cov))) return
        lower = 0.0_dp
        do i=1,ndim
            do j=1,i
                s = cov(i,j)
                do k=1,j-1
                    s = s - lower(i,k)*lower(j,k)
                end do
                if (i==j) then
                    if (.not. ieee_is_finite(s)) return
                    if (s <= 64.0_dp*epsilon(s)) return
                    lower(i,j) = sqrt(s)
                else
                    lower(i,j) = s/lower(j,j)
                end if
            end do
        end do
        lower_inv = 0.0_dp
        do j=1,ndim
            do i=j,ndim
                s = 0.0_dp
                if (i==j) s = 1.0_dp
                do k=j,i-1
                    s = s - lower(i,k)*lower_inv(k,j)
                end do
                lower_inv(i,j) = s/lower(i,i)
            end do
        end do
        inverse = matmul(transpose(lower_inv),lower_inv)
        do j=1,ndim
            do i=1,ndim
                inverse(i,j) = (inverse(i,j)/stddev(i))/stddev(j)
            end do
        end do
        ok = all(ieee_is_finite(inverse))
    end subroutine estimate_inverse_covariance

    pure logical function valid_probabilities(u) result(ok)
        real(dp), intent(in) :: u(:) !! probability weights, not necessarily normalized
        ok = .false.
        if (any(u<0.0_dp)) return
        if (.not. any(u>0.0_dp)) return
        ok = .true.
    end function valid_probabilities

    ! The weights are prepared once per matrix, not once per pair.
    pure function prepare_weights(n,w) result(weights)
        integer, intent(in) :: n !! coordinate count
        real(dp), intent(in), optional :: w(:) !! validated optional weights
        real(dp) :: weights(n)
        weights = 1.0_dp
        if (present(w)) weights = w
    end function prepare_weights

    ! The pointwise kernel accepts already validated equal-length inputs.
    ! Normalize norms before summation to avoid unnecessary overflow/underflow.
    pure function distance_kernel(u, v, metric, order, w, variance, inv_cov, log_base) result(dist)
        real(dp), intent(in) :: u(:) !! first observation
        real(dp), intent(in) :: v(:) !! second observation
        character(len=*), intent(in) :: metric !! validated metric name
        real(dp), intent(in) :: order !! validated Minkowski order (default 2)
        real(dp), intent(in) :: w(:) !! validated nonnegative weights
        real(dp), intent(in), optional :: variance(:) !! validated component variance
        real(dp), intent(in), optional :: inv_cov(:,:) !! validated inverse covariance
        real(dp), intent(in) :: log_base !! log(JS base), 1.0 for natural logarithms
        real(dp) :: dist
        real(dp) :: delta(size(u)), a(size(u)), b(size(u)), denom(size(u))
        real(dp) :: peak, su, sv, numerator, denominator, mean_u, mean_v, m, divergence
        real(dp) :: tt, ff, tf, ft, mismatch, contribution
        integer :: i
        select case (metric)
        case ('euclidean', 'sqeuclidean', 'minkowski', 'seuclidean', 'mahalanobis')
            delta = abs(u-v)
            peak = 0.0_dp
            do i=1,size(u)
                if (w(i) > 0.0_dp) peak = max(peak,delta(i))
            end do
            if (peak == 0.0_dp) then
                dist = 0.0_dp
                return
            end if
            a = delta/peak
            if (metric == 'minkowski') then
                dist = peak * sum(w*a**order)**(1.0_dp/order)
            else if (metric == 'seuclidean') then
                dist = peak*sqrt(sum(a*a/variance))
            else if (metric == 'mahalanobis') then
                ! Preserve the sign of differences, which matters for covariance.
                a = (u-v)/peak
                numerator = dot_product(a,matmul(inv_cov,a))
                if (numerator < 0.0_dp) then
                    dist = invalid_value()
                else
                    dist = peak*sqrt(numerator)
                end if
            else
                dist = peak*sqrt(sum(w*a*a))
                if (metric == 'sqeuclidean') dist = dist*dist
            end if
        case ('cityblock')
            dist = sum(w*abs(u-v))
        case ('chebyshev')
            dist = 0.0_dp
            do i=1,size(u)
                if (w(i) > 0.0_dp) dist = max(dist,abs(u(i)-v(i)))
            end do
        case ('cosine','correlation')
            a = u
            b = v
            if (metric == 'correlation') then
                ! One positively weighted coordinate (or a constant vector)
                ! has zero weighted variance, regardless of rounding in the mean.
                if (count(w > 0.0_dp) < 2) then
                    dist = invalid_value()
                    return
                end if
                if (maxval(u,mask=w>0.0_dp) == minval(u,mask=w>0.0_dp) .or. &
                    maxval(v,mask=w>0.0_dp) == minval(v,mask=w>0.0_dp)) then
                    dist = invalid_value()
                    return
                end if
                denominator = sum(w)
                if (denominator == 0.0_dp) then
                    dist = invalid_value()
                    return
                end if
                mean_u = sum(w*u)/denominator
                mean_v = sum(w*v)/denominator
                a = u-mean_u
                b = v-mean_v
            end if
            su = maxval(abs(a))
            sv = maxval(abs(b))
            if (su == 0.0_dp .or. sv == 0.0_dp) then
                dist = invalid_value()
                return
            end if
            a = a/su
            b = b/sv
            denominator = sqrt(sum(w*a*a))*sqrt(sum(w*b*b))
            if (denominator == 0.0_dp) then
                dist = invalid_value()
                return
            end if
            dist = max(0.0_dp,min(2.0_dp,1.0_dp-sum(w*a*b)/denominator))
        case ('hamming')
            denominator = sum(w)
            if (denominator == 0.0_dp) then
                dist = invalid_value()
            else
                dist = sum(w,mask=u/=v)/denominator
            end if
        case ('jaccard')
            denominator = sum(w,mask=(u/=0.0_dp .or. v/=0.0_dp))
            if (denominator == 0.0_dp) then
                dist = 0.0_dp
            else
                dist = sum(w,mask=((u/=0.0_dp) .neqv. (v/=0.0_dp)))/denominator
            end if
        case ('canberra')
            denom = abs(u)+abs(v)
            dist = 0.0_dp
            do i=1,size(u)
                if (denom(i) > 0.0_dp) dist = dist + w(i)*abs(u(i)-v(i))/denom(i)
            end do
        case ('braycurtis')
            numerator = sum(w*abs(u-v))
            denominator = sum(w*abs(u+v))
            if (denominator == 0.0_dp) then
                dist = invalid_value()
                if (numerator > 0.0_dp) dist = ieee_value(0.0_dp, ieee_positive_inf)
            else
                dist = numerator/denominator
            end if
        case ('dice', 'rogerstanimoto', 'russellrao', 'sokalsneath', 'yule')
            ! Weighted Boolean contingency table: TT, FF, TF, FT.
            ! The native real-valued API uses the bulk SciPy convention of
            ! converting every nonzero value to true (not multiplying values).
            ! Rescale positive weights first, avoiding overflow when totals of
            ! otherwise finite weights exceed the largest representable real.
            peak = maxval(w)
            if (peak == 0.0_dp) then
                dist = invalid_value()
                if (metric == 'yule') dist = 0.0_dp
                return
            end if
            tt = 0.0_dp
            ff = 0.0_dp
            tf = 0.0_dp
            ft = 0.0_dp
            do i=1,size(u)
                if (w(i) == 0.0_dp) cycle
                contribution = w(i)/peak
                if (u(i) /= 0.0_dp) then
                    if (v(i) /= 0.0_dp) then
                        tt = tt + contribution
                    else
                        tf = tf + contribution
                    end if
                else
                    if (v(i) /= 0.0_dp) then
                        ft = ft + contribution
                    else
                        ff = ff + contribution
                    end if
                end if
            end do
            denominator = (tt+ff)+(tf+ft)
            tt = tt/denominator
            ff = ff/denominator
            tf = tf/denominator
            ft = ft/denominator
            mismatch = tf+ft
            select case (metric)
            case ('dice')
                denominator = 2.0_dp*tt+mismatch
                dist = invalid_value()
                if (denominator > 0.0_dp) dist = mismatch/denominator
            case ('rogerstanimoto')
                dist = (2.0_dp*mismatch)/(tt+ff+2.0_dp*mismatch)
            case ('russellrao')
                dist = 1.0_dp-tt
            case ('sokalsneath')
                denominator = tt+2.0_dp*mismatch
                dist = invalid_value()
                if (denominator > 0.0_dp) dist = 2.0_dp*mismatch/denominator
            case ('yule')
                ! SciPy returns zero when neither off-diagonal product nor
                ! the diagonal product is positive (no defined odds ratio).
                numerator = tf*ft
                denominator = tt*ff+numerator
                dist = 0.0_dp
                if (denominator > 0.0_dp) dist = 2.0_dp*numerator/denominator
            end select
        case ('jensenshannon')
            ! Scale before normalizing: large finite weights need not overflow.
            a = (u/maxval(u))
            b = (v/maxval(v))
            a = a/sum(a)
            b = b/sum(b)
            divergence = 0.0_dp
            do i=1,size(u)
                m = 0.5_dp*(a(i)+b(i))
                if (a(i)>0.0_dp) divergence = divergence + a(i)*log(a(i)/m)
                if (b(i)>0.0_dp) divergence = divergence + b(i)*log(b(i)/m)
            end do
            dist = sqrt(max(0.0_dp,0.5_dp*divergence/log_base))
        case default
            dist = invalid_value()
        end select
    end function distance_kernel

    function point_distance(u, v, metric, p, w, status, variance, inv_cov, base) result(dist)
        real(dp), intent(in) :: u(:) !! first coordinate vector
        real(dp), intent(in) :: v(:) !! second coordinate vector, same length as u
        character(len=*), intent(in), optional :: metric !! supported metric, default euclidean
        real(dp), intent(in), optional :: p !! Minkowski exponent > 0, default 2
        real(dp), intent(in), optional :: w(:) !! nonnegative coordinate weights for supported metrics
        integer, intent(out), optional :: status !! zero on success; one for invalid arguments
        real(dp), intent(in), optional :: variance(:) !! required for pointwise seuclidean
        real(dp), intent(in), optional :: inv_cov(:,:) !! required for pointwise mahalanobis
        real(dp), intent(in), optional :: base !! Jensen-Shannon logarithm base > 1 (default e)
        real(dp) :: dist
        real(dp) :: weights(size(u)), order, logarithm_base
        character(len=:), allocatable :: name
        integer :: code
        code = spatial_distance_invalid_input
        dist = invalid_value()
        name = 'euclidean'
        if (present(metric)) name = trim(adjustl(metric))
        if (size(v) /= size(u)) goto 900
        if (.not. all(ieee_is_finite(u)) .or. .not. all(ieee_is_finite(v))) goto 900
        if (.not. valid_options(size(u),name,p,w,variance,inv_cov,base)) goto 900
        if (name=='seuclidean' .and. .not. present(variance)) goto 900
        if (name=='mahalanobis' .and. .not. present(inv_cov)) goto 900
        if (name=='jensenshannon') then
            if (.not. valid_probabilities(u) .or. .not. valid_probabilities(v)) goto 900
        end if
        order = 2.0_dp
        if (present(p)) order = p
        weights = prepare_weights(size(u),w)
        logarithm_base = 1.0_dp
        if (present(base)) logarithm_base = log(base)
        dist = distance_kernel(u,v,name,order,weights,variance,inv_cov,logarithm_base)
        code = spatial_distance_success
900     if (present(status)) status = code
    end function point_distance

    function cdist(xa, xb, metric, p, w, status, variance, inv_cov, base) result(dist)
        real(dp), intent(in) :: xa(:,:) !! first points (observations, dimensions)
        real(dp), intent(in) :: xb(:,:) !! second points (observations, dimensions)
        character(len=*), intent(in), optional :: metric !! supported metric, default euclidean
        real(dp), intent(in), optional :: p !! Minkowski exponent > 0, default 2
        real(dp), intent(in), optional :: w(:) !! nonnegative coordinate weights for supported metrics
        integer, intent(out), optional :: status !! zero on success; one for invalid arguments
        real(dp), intent(in), optional :: variance(:) !! seuclidean component variance; estimated if absent
        real(dp), intent(in), optional :: inv_cov(:,:) !! Mahalanobis inverse covariance; estimated if absent
        real(dp), intent(in), optional :: base !! Jensen-Shannon logarithm base >1 (default e)
        real(dp), allocatable :: dist(:,:)
        real(dp), allocatable :: variances(:), inverse(:,:)
        real(dp) :: weights(size(xa,2)), order, logarithm_base
        character(len=:), allocatable :: name
        integer :: i,j,code,ndim
        logical :: ok
        code = spatial_distance_invalid_input
        name = 'euclidean'
        if (present(metric)) name = trim(adjustl(metric))
        ndim = size(xa,2)
        if (ndim /= size(xb,2)) goto 900
        if (.not. valid_options(ndim,name,p,w,variance,inv_cov,base)) goto 900
        if (.not. all(ieee_is_finite(xa)) .or. .not. all(ieee_is_finite(xb))) goto 900
        if (name=='jensenshannon') then
            do i=1,size(xa,1)
                if (.not. valid_probabilities(xa(i,:))) goto 900
            end do
            do i=1,size(xb,1)
                if (.not. valid_probabilities(xb(i,:))) goto 900
            end do
        end if
        if (name=='seuclidean') then
            allocate(variances(ndim))
            if (present(variance)) then
                variances = variance
            else
                call estimate_variance(xa, variances, ok, xb)
                if (.not. ok) goto 900
            end if
        else if (name=='mahalanobis') then
            allocate(inverse(ndim,ndim))
            if (present(inv_cov)) then
                inverse = inv_cov
            else
                call estimate_inverse_covariance(xa, inverse, ok, xb)
                if (.not. ok) goto 900
            end if
        end if
        weights = prepare_weights(ndim,w)
        order = 2.0_dp
        if (present(p)) order = p
        logarithm_base = 1.0_dp
        if (present(base)) logarithm_base = log(base)
        allocate(dist(size(xa,1),size(xb,1)))
        do j=1,size(xb,1)
            do i=1,size(xa,1)
                dist(i,j) = distance_kernel(xa(i,:),xb(j,:),name,order,weights, &
                    variances,inverse,logarithm_base)
            end do
        end do
        code = spatial_distance_success
900     if (code /= spatial_distance_success) allocate(dist(0,0))
        if (present(status)) status = code
    end function cdist

    function pdist(x, metric, p, w, status, variance, inv_cov, base) result(dist)
        real(dp), intent(in) :: x(:,:) !! points (observations, dimensions)
        character(len=*), intent(in), optional :: metric !! supported metric, default euclidean
        real(dp), intent(in), optional :: p !! Minkowski exponent > 0, default 2
        real(dp), intent(in), optional :: w(:) !! nonnegative coordinate weights for supported metrics
        integer, intent(out), optional :: status !! zero on success; one for invalid arguments
        real(dp), intent(in), optional :: variance(:) !! seuclidean component variance; estimated if absent
        real(dp), intent(in), optional :: inv_cov(:,:) !! Mahalanobis inverse covariance; estimated if absent
        real(dp), intent(in), optional :: base !! Jensen-Shannon logarithm base >1 (default e)
        real(dp), allocatable :: dist(:)
        real(dp), allocatable :: variances(:), inverse(:,:)
        real(dp) :: weights(size(x,2)), order, logarithm_base
        character(len=:), allocatable :: name
        integer :: i,j,k,n,code,ndim
        logical :: ok
        code = spatial_distance_invalid_input
        name = 'euclidean'
        if (present(metric)) name = trim(adjustl(metric))
        n = size(x,1)
        ndim = size(x,2)
        if (.not. valid_options(ndim,name,p,w,variance,inv_cov,base)) goto 900
        if (.not. all(ieee_is_finite(x))) goto 900
        if ((int(n,int64)*(int(n,int64)-1_int64))/2_int64 > int(huge(n),int64)) goto 900
        if (name=='jensenshannon') then
            do i=1,n
                if (.not. valid_probabilities(x(i,:))) goto 900
            end do
        end if
        if (name=='seuclidean') then
            allocate(variances(ndim))
            if (present(variance)) then
                variances = variance
            else
                call estimate_variance(x, variances, ok)
                if (.not. ok) goto 900
            end if
        else if (name=='mahalanobis') then
            allocate(inverse(ndim,ndim))
            if (present(inv_cov)) then
                inverse = inv_cov
            else
                call estimate_inverse_covariance(x, inverse, ok)
                if (.not. ok) goto 900
            end if
        end if
        weights = prepare_weights(ndim,w)
        order = 2.0_dp
        if (present(p)) order = p
        logarithm_base = 1.0_dp
        if (present(base)) logarithm_base = log(base)
        allocate(dist(n*(n-1)/2))
        k = 0
        do i=1,n-1
            do j=i+1,n
                k = k+1
                dist(k) = distance_kernel(x(i,:),x(j,:),name,order,weights, &
                    variances,inverse,logarithm_base)
            end do
        end do
        code = spatial_distance_success
900     if (code /= spatial_distance_success) allocate(dist(0))
        if (present(status)) status = code
    end function pdist

    ! SciPy's condensed layout visits pairs in lexicographic (i,j) order:
    ! (1,2),(1,3),...,(1,n),(2,3),...,(n-1,n).
    function squareform_to_matrix(v,status) result(dist)
        real(dp), intent(in) :: v(:) !! condensed upper triangle, length n*(n-1)/2
        integer, intent(out), optional :: status !! zero if length is triangular
        real(dp), allocatable :: dist(:,:)
        integer(int64) :: n, nvalue
        integer :: i,j,k,code
        code = spatial_distance_invalid_input
        nvalue = int(size(v),int64)
        n = int((1.0_dp+sqrt(1.0_dp+8.0_dp*real(nvalue,dp)))/2.0_dp,int64)
        ! Correct rounding when a large condensed vector is near a triangular length.
        do while (n*(n-1_int64)/2_int64 < nvalue)
            n = n+1_int64
        end do
        do while (n>1_int64)
            if (n*(n-1_int64)/2_int64 <= nvalue) exit
            n = n-1_int64
        end do
        if (n*(n-1_int64)/2_int64 /= nvalue) goto 900
        if (n>int(huge(i),int64)) goto 900
        allocate(dist(int(n),int(n)))
        dist = 0.0_dp
        k = 0
        do i=1,int(n)-1
            do j=i+1,int(n)
                k = k+1
                dist(i,j) = v(k)
                dist(j,i) = v(k)
            end do
        end do
        code = spatial_distance_success
900     if (code /= spatial_distance_success) allocate(dist(0,0))
        if (present(status)) status = code
    end function squareform_to_matrix

    function squareform_to_vector(m,checks,status) result(dist)
        real(dp), intent(in) :: m(:,:) !! symmetric square distance matrix
        logical, intent(in), optional :: checks !! verify zero diagonal and symmetry, default true
        integer, intent(out), optional :: status !! zero if matrix accepted
        real(dp), allocatable :: dist(:)
        integer :: n,i,j,k,code
        logical :: verify
        code = spatial_distance_invalid_input
        verify = .true.
        if (present(checks)) verify = checks
        n = size(m,1)
        if (size(m,2) /= n) goto 900
        if (verify) then
            if (.not. all(ieee_is_finite(m))) goto 900
            do i=1,n
                if (m(i,i) /= 0.0_dp) goto 900
                do j=i+1,n
                    if (m(i,j) /= m(j,i)) goto 900
                end do
            end do
        end if
        if ((int(n,int64)*(int(n,int64)-1_int64))/2_int64 > int(huge(n),int64)) goto 900
        allocate(dist(n*(n-1)/2))
        k = 0
        do i=1,n-1
            do j=i+1,n
                k = k+1
                dist(k) = m(i,j)
            end do
        end do
        code = spatial_distance_success
900     if (code /= spatial_distance_success) allocate(dist(0))
        if (present(status)) status = code
    end function squareform_to_vector

end module scifort_spatial_distance
