! SPDX-License-Identifier: MIT
! Complex generic derivatives and indexed label reductions.
program test_ndimage_complex_derivatives_label_lookup
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_ndimage, only : generic_gradient_magnitude_nd, generic_laplace_nd, &
        correlate1d_nd, sum_labels_nd, mean_labels_nd, variance_labels_nd, &
        standard_deviation_labels_nd, center_of_mass_nd, count_labels_nd, &
        histogram_labels_nd, extrema_labels_nd, bounding_boxes_nd
    implicit none
    complex(dp) :: x(5,4,3),g(5,4,3),lap(5,4,3),work(5,4,3),expected(5,4,3)
    complex(dp) :: z15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    complex(dp) :: out15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    complex(dp) :: sliced(10,8,3)
    real(dp) :: data(5,4),totals(8),means(8),vars(8),stds(8)
    real(dp) :: minval(8),maxima(8),centers(2,8),sumcheck,meancheck,varcheck
    integer :: labs(5,4),requests(8),counts(8),hist(7,8)
    integer :: labels15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1),tiny_counts(3)
    real(dp) :: real15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1),tiny_totals(3)
    integer :: minpos(2,8),maxpos(2,8),lower(2,8),upper(2,8)
    integer :: i,j,k,t,s,n,bin,coords(2),small(2),large(2)
    character(len=16) :: boundary(2)
    real(dp), parameter :: tol=1.0e-11_dp

    do k=1,3
        do j=1,4
            do i=1,5
                x(i,j,k)=cmplx(sin(0.31_dp*i+0.2_dp*j+0.47_dp*k), &
                    cos(0.11_dp*i-0.44_dp*j+0.13_dp*k),dp)
            end do
        end do
    end do
    boundary=['constant        ','nearest         ']
    call generic_gradient_magnitude_nd(x,complex_first,g,s, &
        modes=boundary,axes=[3,1],cval=cmplx(0.37_dp,-0.6_dp,dp))
    call check(s==0,'complex selected gradient')
    expected=cmplx(0.0_dp,0.0_dp,dp)
    call complex_first(x,3,work,trim(boundary(1)),cmplx(0.37_dp,-0.6_dp,dp),s)
    expected=expected+work**2
    call complex_first(x,1,work,trim(boundary(2)),cmplx(0.37_dp,-0.6_dp,dp),s)
    expected=sqrt(expected+work**2)
    call check(maxval(abs(expected-g))<tol,'complex gradient sum of squares')
    call generic_laplace_nd(x,complex_second,lap,s, &
        modes=boundary,axes=[3,1],cval=cmplx(0.37_dp,-0.6_dp,dp))
    call check(s==0,'complex selected Laplace')
    call complex_second(x,3,expected,trim(boundary(1)),cmplx(0.37_dp,-0.6_dp,dp),s)
    call complex_second(x,1,work,trim(boundary(2)),cmplx(0.37_dp,-0.6_dp,dp),s)
    expected=expected+work
    call check(maxval(abs(expected-lap))<tol,'complex Laplace sum')
    call generic_laplace_nd(x,complex_second,lap,s,axes=[1,1])
    call check(s/=0,'reject repeated axes')
    call generic_laplace_nd(x,complex_second,lap,s,axes=[4])
    call check(s/=0,'reject invalid axis')
    call generic_gradient_magnitude_nd(x,complex_first,g,s,modes=['bad_mode'])
    call check(s/=0,'reject invalid modes')
    call generic_laplace_nd(x,complex_fail,lap,s)
    call check(s==1,'propagate callback failure')
    call generic_gradient_magnitude_nd(x,complex_first,g,s,mode='wrap')
    call check(s==0,'complex gradient all axes')
    call generic_laplace_nd(x,complex_second,sliced(1:10:2,1:8:2,:),s,mode='reflect')
    call check(s==0,'noncontiguous complex output')
    z15=cmplx(2.0_dp,-1.0_dp,dp)
    z15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=cmplx(3.0_dp,4.0_dp,dp)
    call generic_gradient_magnitude_nd(z15,complex_first,out15,s,axes=[1])
    call check(s==0,'rank-15 complex gradient')
    call generic_laplace_nd(z15,complex_second,out15,s,axes=[1])
    call check(s==0,'rank-15 complex Laplace')

    do j=1,4
        do i=1,5
            data(i,j)=real(3*i-j,dp)/5.0_dp+cos(real(i*j,dp))
            select case(modulo(2*i+3*j,6))
            case(0)
                labs(i,j)=-7
            case(1)
                labs(i,j)=0
            case(2)
                labs(i,j)=500000001
            case default
                labs(i,j)=3
            end select
        end do
    end do
    requests=[3,-7,0,3,500000001,-999,-7,huge(1)]
    call sum_labels_nd(data,labs,requests,totals,s)
    call check(s==0,'label sums')
    call mean_labels_nd(data,labs,requests,means,s)
    call check(s==0,'label means')
    call variance_labels_nd(data,labs,requests,vars,s)
    call check(s==0,'label variances')
    call standard_deviation_labels_nd(data,labs,requests,stds,s)
    call check(s==0,'label standard deviations')
    call center_of_mass_nd(data,labs,requests,centers,s)
    call check(s==0,'label centers')
    call count_labels_nd(labs,requests,counts,s)
    call check(s==0,'label counts')
    call histogram_labels_nd(data,labs,requests,-5.0_dp,5.0_dp,hist,s)
    call check(s==0,'label histograms')
    call extrema_labels_nd(data,labs,requests,minval,maxima,minpos,maxpos,s)
    call check(s==0,'label extrema')
    call bounding_boxes_nd(labs,requests,lower,upper,s)
    call check(s==0,'label boxes')
    do t=1,size(requests)
        n=0
        sumcheck=0.0_dp
        meancheck=0.0_dp
        varcheck=0.0_dp
        small=huge(1)
        large=0
        do j=1,4
            do i=1,5
                if (labs(i,j)/=requests(t)) cycle
                n=n+1
                sumcheck=sumcheck+data(i,j)
                small=min(small,[i,j])
                large=max(large,[i,j])
            end do
        end do
        call check(counts(t)==n,'count matches simple scan')
        call check(abs(totals(t)-sumcheck)<tol,'sum matches simple scan')
        call check(sum(hist(:,t))==n,'histogram mass matches count')
        if (n==0) then
            call check(all(lower(:,t)==0).and.all(upper(:,t)==0),'missing box')
            call check(all(minpos(:,t)==0).and.all(maxpos(:,t)==0),'missing extrema')
            call check(ieee_is_nan(means(t)).and.ieee_is_nan(vars(t)),'missing statistics')
        else
            meancheck=sumcheck/real(n,dp)
            do j=1,4
                do i=1,5
                    if (labs(i,j)/=requests(t)) cycle
                    varcheck=varcheck+(data(i,j)-meancheck)**2
                end do
            end do
            varcheck=varcheck/real(n,dp)
            call check(abs(means(t)-meancheck)<tol,'mean matches direct calculation')
            call check(abs(vars(t)-varcheck)<tol,'variance matches direct calculation')
            call check(abs(stds(t)-sqrt(varcheck))<tol,'standard deviation matches direct calculation')
            call check(all(lower(:,t)==small).and.all(upper(:,t)==large),'box matches direct calculation')
            if (abs(sumcheck)>1.0e-9_dp) then
                do k=1,2
                    meancheck=0.0_dp
                    do j=1,4
                        do i=1,5
                            if (labs(i,j)/=requests(t)) cycle
                            coords=[i,j]
                            meancheck=meancheck+data(i,j)*real(coords(k),dp)
                        end do
                    end do
                    call check(abs(centers(k,t)-meancheck/sumcheck)<tol,'center of mass')
                end do
            end if
        end if
    end do
    call check(all(counts([1,2])==counts([4,7])),'duplicate-label counts')
    call check(maxval(abs(totals([1,2])-totals([4,7])))<tol,'duplicate-label sums')
    call check(all(hist(:,1)==hist(:,4)),'duplicate-label histogram')
    call check(all(minpos(:,2)==minpos(:,7)),'duplicate-label extrema')
    call check(all(lower(:,2)==lower(:,7)),'duplicate-label bounds')
    labels15=huge(1)
    labels15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=-huge(1)
    real15=3.0_dp
    call count_labels_nd(labels15,[-huge(1),huge(1),-huge(1)],tiny_counts,s)
    call check(s==0.and.all(tiny_counts==1),'rank-15 indexed label counts')
    call sum_labels_nd(real15,labels15,[-huge(1),huge(1),-huge(1)],tiny_totals,s)
    call check(s==0.and.all(abs(tiny_totals-3.0_dp)<tol),'rank-15 indexed label sums')
    call count_labels_nd(labels15,[integer::],tiny_counts(:0),s)
    call check(s==0,'zero label requests')
    call sum_labels_nd(real15,labels15,[integer::],tiny_totals(:0),s)
    call check(s==0,'zero-label sum requests')
    print '(a)', 'PASS: complex derivatives and indexed label statistics'
contains
    subroutine check(test,message)
        logical, intent(in) :: test !! assertion
        character(len=*), intent(in) :: message !! assertion description
        if (.not.test) then
            print '(2a)', 'FAIL: ',message
            error stop 1
        end if
    end subroutine check

    subroutine complex_first(input,axis,output,mode,cval,status)
        complex(dp), intent(in) :: input(..) !! input field
        integer, intent(in) :: axis !! differentiation dimension
        complex(dp), intent(out) :: output(..) !! first derivative
        character(len=*), intent(in) :: mode !! boundary mode
        complex(dp), intent(in) :: cval !! boundary value
        integer, intent(out) :: status !! success/error
        call correlate1d_nd(input,cmplx([-0.5_dp,0.0_dp,0.5_dp],0.0_dp,dp), &
            axis,output,status=status,mode=mode,cval=cval)
    end subroutine complex_first

    subroutine complex_second(input,axis,output,mode,cval,status)
        complex(dp), intent(in) :: input(..) !! input field
        integer, intent(in) :: axis !! differentiation dimension
        complex(dp), intent(out) :: output(..) !! second derivative
        character(len=*), intent(in) :: mode !! boundary mode
        complex(dp), intent(in) :: cval !! boundary value
        integer, intent(out) :: status !! success/error
        call correlate1d_nd(input,cmplx([1.0_dp,-2.0_dp,1.0_dp],0.0_dp,dp), &
            axis,output,status=status,mode=mode,cval=cval)
    end subroutine complex_second

    subroutine complex_fail(input,axis,output,mode,cval,status)
        complex(dp), intent(in) :: input(..) !! input field
        integer, intent(in) :: axis !! ignored axis
        complex(dp), intent(out) :: output(..) !! ignored output
        character(len=*), intent(in) :: mode !! ignored mode
        complex(dp), intent(in) :: cval !! ignored value
        integer, intent(out) :: status !! intentional failure
        status=1
    end subroutine complex_fail
end program test_ndimage_complex_derivatives_label_lookup
