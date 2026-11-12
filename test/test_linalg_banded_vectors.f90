! SPDX-License-Identifier: MIT
program test_linalg_banded_vectors
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape
    use scifort_linalg_banded_eigen, only : eig_banded
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none
    real(dp) :: a(7,7), ab(3,7), w(3), v(7,3), short(1)
    real(dp) :: full(7), scale
    complex(dp) :: c(7,7), cb(3,7), z(7,3)
    real(dp) :: wd(7), zero(3,7)
    real(dp), allocatable :: large(:, :), largez(:, :)
    real(dp) :: largew(3)
    integer :: status, nfound, i, j, k, n
    a=0.0_dp
    do i=1,7
        a(i,i)=3.0_dp+0.07_dp*i
        if (i<7) then
            a(i+1,i)=-0.8_dp
            a(i,i+1)=-0.8_dp
        end if
        if (i<6) then
            a(i+2,i)=0.2_dp
            a(i,i+2)=0.2_dp
        end if
    end do
    ab=0.0_dp
    do j=1,7
        do k=0,min(2,7-j)
            ab(k+1,j)=a(j+k,j)
        end do
    end do
    call eig_banded(ab,w,status,v,lower=.true.,subset_by_index=[2,4], &
        num_selected=nfound,method='inverse')
    call expect(status==linalg_status_success .and. nfound==3,'selected real status')
    call eig_banded(ab,full,status,lower=.true.)
    call expect(status==linalg_status_success,'full real status')
    call expect(maxval(abs(w-full(2:4)))<1.0e-9_dp,'selected real eigenvalues')
    do j=1,3
        call expect(sqrt(sum((matmul(a,v(:,j))-w(j)*v(:,j))**2))<1.0e-9_dp, &
            'selected real residual')
    end do
    call expect(maxval(abs(matmul(transpose(v),v)-eye(3)))<2.0e-9_dp, &
        'selected real orthogonality')

    ! The real-LU and complex-LU paths should agree for a real symmetric
    ! input while retaining separate arithmetic and storage conventions.
    cb=cmplx(ab,0.0_dp,dp)
    call eig_banded(cb,wd(1:3),status,z,lower=.true., &
        subset_by_index=[2,4],num_selected=nfound,method='inverse')
    call expect(status==linalg_status_success .and. nfound==3, &
        'real-vs-complex LU status')
    call expect(maxval(abs(wd(1:3)-w))<1.0e-10_dp, &
        'real-vs-complex LU eigenvalues')
    do j=1,3
        call expect(abs(dot_product(cmplx(v(:,j),0.0_dp,dp),z(:,j))) &
            >1.0_dp-5.0e-8_dp,'real-vs-complex LU eigenvectors')
    end do

    ! Upper storage and lower/upper padding equivalence.
    do j=1,7
        do k=0,min(2,j-1)
            ab(3-k,j)=a(j-k,j)
        end do
    end do
    ab(1,1)=ieee_value(0.0_dp,ieee_quiet_nan)
    call eig_banded(ab,w,status,v,subset_by_index=[5,7],num_selected=nfound,method='inverse')
    call expect(status==linalg_status_success .and. nfound==3,'real upper status')
    call expect(maxval(abs(w-full(5:7)))<1.0e-9_dp,'real upper values')
    call expect(maxval(abs(matmul(a,v)-v*spread(w,1,7)))<1.0e-9_dp,'real upper residual')

    c=cmplx(a,0.0_dp,dp)
    do j=1,7
        do i=j+1,min(7,j+2)
            c(i,j)=c(i,j)+cmplx(0.0_dp,0.25_dp*real(i-j,dp),dp)
            c(j,i)=conjg(c(i,j))
        end do
    end do
    cb=cmplx(0.0_dp,0.0_dp,dp)
    do j=1,7
        do k=0,min(2,7-j)
            cb(k+1,j)=c(j+k,j)
        end do
    end do
    call eig_banded(cb,w,status,z,lower=.true.,subset_by_index=[3,5], &
        num_selected=nfound,method='inverse')
    call expect(status==linalg_status_success .and. nfound==3,'complex selected status')
    call eig_banded(cb,wd,status,lower=.true.)
    call expect(status==linalg_status_success,'complex full status')
    call expect(maxval(abs(w-wd(3:5)))<1.0e-9_dp,'complex selected values')
    do j=1,3
        call expect(sqrt(sum(abs(matmul(c,z(:,j))-w(j)*z(:,j))**2))<2.0e-9_dp, &
            'complex selected residual')
    end do
    call expect(maxval(abs(matmul(conjg(transpose(z)),z)-cmplx(eye(3),0.0_dp,dp))) &
        <2.0e-9_dp,'complex orthogonality')

    call eig_banded(cb,w,status,z,lower=.true., &
        subset_by_value=[0.5_dp*(wd(1)+wd(2)),0.5_dp*(wd(4)+wd(5))], &
        num_selected=nfound,method='inverse')
    call expect(status==linalg_status_success .and. nfound==3, &
        'complex inverse value selection')
    call expect(maxval(abs(w-wd(2:4)))<1.0e-8_dp, &
        'complex inverse value eigenvalues')
    call eig_banded(cb,w,status,z,lower=.true.,subset_by_value=[wd(2),wd(4)], &
        num_selected=nfound,method='dense')
    call expect(status==linalg_status_success .and. nfound>=1,'dense value selection')
    call eig_banded(cb,w,status,z,lower=.true., &
        subset_by_value=[100.0_dp,101.0_dp], &
        num_selected=nfound,method='inverse')
    call expect(status==linalg_status_success .and. nfound==0, 'empty inverse selection')
    call expect(all(w==0.0_dp) .and. all(z==cmplx(0.0_dp,0.0_dp,dp)), &
        'empty inverse output initialized')
    call eig_banded(cb,w,status,z(:,1:1),lower=.true.,subset_by_index=[2,4], &
        num_selected=nfound,method='inverse')
    call expect(status==linalg_status_invalid_shape .and. nfound==3, &
        'too few selected eigenvector columns')
    call eig_banded(cb,short,status,lower=.true.,subset_by_index=[2,4], &
        num_selected=nfound,method='inverse')
    call expect(status==linalg_status_invalid_shape .and. nfound==3, &
        'too small selected value output')
    call eig_banded(cb,w,status,z,lower=.true.,subset_by_index=[2,4], &
        subset_by_value=[0.0_dp,10.0_dp])
    call expect(status==linalg_status_invalid_shape,'mutually exclusive subsets')
    call eig_banded(cb,w,status,z,method='unknown')
    call expect(status==linalg_status_invalid_shape,'unknown method')
    call eig_banded(cb,w,status,z,method='inverse')
    call expect(status==linalg_status_invalid_shape,'inverse requires selection')
    call eig_banded(cb,w,status,z,lower=.true.,subset_by_index=[1,3],tol=-1.0_dp)
    call expect(status==linalg_status_invalid_shape,'invalid bisection tolerance')

    ! Exact repeated eigenvalues are resolved without LU breakdowns.
    zero=0.0_dp
    zero(1,:)=[2.0_dp,2.0_dp,2.0_dp,-1.0_dp,3.0_dp,3.0_dp,3.0_dp]
    call eig_banded(zero,w,status,v,lower=.true.,subset_by_index=[2,4], &
        num_selected=nfound,method='inverse')
    call expect(status==linalg_status_success .and. nfound==3,'diagonal multiplicity')
    call expect(maxval(abs(matmul(transpose(v),v)-eye(3)))<1.0e-13_dp, &
        'diagonal orthogonality')
    zero(1,:)=0.0_dp
    call eig_banded(zero,w,status,v,lower=.true.,subset_by_index=[2,4], &
        num_selected=nfound,method='inverse')
    call expect(status==linalg_status_success,'zero spectrum')

    ! Demonstrate O(n*kd+n*k) workspace on order 6000, bandwidth 1.
    n=6000
    allocate(large(2,n),largez(n,3))
    large=0.0_dp
    large(1,:)=2.0_dp
    large(2,1:n-1)=-1.0_dp
    call eig_banded(large,largew,status,largez,lower=.true., &
        subset_by_index=[2,4],num_selected=nfound,method='inverse')
    call expect(status==linalg_status_success .and. nfound==3,'large narrow band')
    do j=1,3
        scale=0.0_dp
        do i=1,n
            if (i==1) then
                scale=max(scale,abs((2.0_dp-largew(j))*largez(i,j)-largez(i+1,j)))
            else if (i==n) then
                scale=max(scale,abs((2.0_dp-largew(j))*largez(i,j)-largez(i-1,j)))
            else
                scale=max(scale,abs((2.0_dp-largew(j))*largez(i,j) &
                    -largez(i-1,j)-largez(i+1,j)))
            end if
        end do
        call expect(scale<3.0e-8_dp,'large eigenvector residual')
    end do
    write(*,'(a)') 'selected banded eigenvector native tests passed'
contains
    subroutine expect(ok, description)
        logical, intent(in) :: ok !! required condition
        character(len=*), intent(in) :: description !! assertion label
        if (.not. ok) then
            print '(a)', 'FAILED: '//description
            error stop 1
        end if
    end subroutine expect
    function eye(n) result(out)
        integer, intent(in) :: n !! square matrix order
        real(dp) :: out(n,n)
        integer :: index
        out=0.0_dp
        do index=1,n
            out(index,index)=1.0_dp
        end do
    end function eye
end program test_linalg_banded_vectors
