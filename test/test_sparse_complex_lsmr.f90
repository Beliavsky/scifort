program test_sparse_complex_lsmr
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csc_complex_matrix, &
        csr_complex_from_coo, csr_complex_to_csc, csr_complex_matvec
    use scifort_sparse_linalg_complex_lsmr, only : lsmr_complex, complex_lsmr_invalid_input
    implicit none
    type(csr_complex_matrix) :: a, tall, wide, large
    type(csc_complex_matrix) :: ac
    complex(dp), allocatable :: b(:), x(:), xtrue(:), centers(:), vals(:), rhs(:)
    complex(dp) :: diag(3), expected(3), noisy(4), wideb(2)
    integer, allocatable :: rows(:), cols(:)
    integer :: info, status, iterations, j, i
    real(dp) :: normr, normar, norma, conda, normx, damping

    rows = [1,2,3,4,1,2,3]
    cols = [1,2,3,3,1,2,3]
    vals = [cmplx(1.0_dp,0.5_dp,dp),cmplx(2.0_dp,-0.2_dp,dp), &
            cmplx(0.5_dp,1.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp), &
            cmplx(0.1_dp,-0.1_dp,dp),cmplx(0.2_dp,0.2_dp,dp), &
            cmplx(0.3_dp,-0.1_dp,dp)]
    a = csr_complex_from_coo(4,3,rows,cols,vals,info)
    call require(info == 0,'CSR construction')
    ac = csr_complex_to_csc(a,info)
    call require(info == 0,'CSC conversion')
    xtrue = [cmplx(0.5_dp,-0.7_dp,dp), cmplx(-1.0_dp,0.2_dp,dp), cmplx(0.4_dp,1.3_dp,dp)]
    b = csr_complex_matvec(a,xtrue,info)
    call require(info == 0,'CSR matvec')
    call lsmr_complex(a,b,x,status,itn=iterations,normr=normr,normar=normar, &
        norma=norma,conda=conda,normx=normx,atol=1.0e-12_dp,btol=1.0e-12_dp, &
        maxiter=60)
    call require(status == 1 .or. status == 2,'CSR consistent solve status')
    call require(maxval(abs(x-xtrue)) < 1.0e-9_dp,'CSR consistent solution')
    call require(normr < 1.0e-9_dp .and. normar < 1.0e-9_dp,'CSR exact residuals')
    call require(iterations > 0 .and. iterations <= 60,'iterations returned')
    call require(norma > 0.0_dp .and. conda >= 1.0_dp .and. normx > 0.0_dp,'diagnostics')
    call lsmr_complex(ac,b,x,status,atol=1.0e-12_dp,btol=1.0e-12_dp,maxiter=60)
    call require(status == 1 .or. status == 2,'CSC consistent solve status')
    call require(maxval(abs(x-xtrue)) < 1.0e-9_dp,'CSC solution')

    ! Check damping and nonzero x0 against the analytic diagonal ridge solution.
    diag = [cmplx(1.0_dp,0.5_dp,dp),cmplx(2.0_dp,-0.2_dp,dp),cmplx(0.5_dp,1.0_dp,dp)]
    rows = [1,2,3]
    cols = [1,2,3]
    tall = csr_complex_from_coo(4,3,rows,cols,diag,info)
    call require(info == 0,'diagonal rectangular construction')
    noisy = [cmplx(1.2_dp,-0.3_dp,dp),cmplx(-0.9_dp,2.0_dp,dp), &
        cmplx(0.7_dp,0.9_dp,dp),cmplx(1.0_dp,-2.0_dp,dp)]
    centers = [cmplx(0.3_dp,-0.2_dp,dp),cmplx(-0.1_dp,-0.7_dp,dp),cmplx(1.0_dp,1.1_dp,dp)]
    damping = 0.7_dp
    do j = 1,3
        expected(j) = centers(j)+conjg(diag(j))*(noisy(j)-diag(j)*centers(j)) &
            /(abs(diag(j))**2+damping**2)
    end do
    call lsmr_complex(tall,noisy,x,status,damp=damping,x0=centers, &
        atol=1.0e-12_dp,btol=1.0e-12_dp,maxiter=80,normar=normar)
    call require(status == 1 .or. status == 2,'damped centered solve status')
    call require(maxval(abs(x-expected)) < 2.0e-9_dp,'centered Tikhonov solution')
    call require(normar < 1.0e-8_dp,'centered Tikhonov gradient')
    call lsmr_complex(tall,noisy,x,status,atol=1.0e-12_dp,btol=1.0e-12_dp,maxiter=80)
    call require(status == 1 .or. status == 2,'inconsistent solve status')
    call require(maxval(abs(x-noisy(:3)/diag)) < 2.0e-9_dp,'inconsistent min norm')

    ! Underdetermined system: x0 controls the unconstrained components.
    rows = [1,2]
    cols = [1,2]
    vals = [cmplx(1.0_dp,0.5_dp,dp),cmplx(0.8_dp,-0.3_dp,dp)]
    wide = csr_complex_from_coo(2,4,rows,cols,vals,info)
    call require(info == 0,'wide matrix')
    wideb = [cmplx(0.7_dp,-1.1_dp,dp), cmplx(1.0_dp,1.5_dp,dp)]
    centers = [cmplx(0.1_dp,0.2_dp,dp),cmplx(0.2_dp,0.3_dp,dp), &
        cmplx(0.4_dp,-0.5_dp,dp),cmplx(-0.6_dp,0.7_dp,dp)]
    call lsmr_complex(wide,wideb,x,status,x0=centers,atol=1.0e-12_dp, &
        btol=1.0e-12_dp,maxiter=80)
    call require(status == 1 .or. status == 2,'underdetermined status')
    call require(abs(x(1)-wideb(1)/vals(1)) < 2.0e-10_dp,'wide first')
    call require(abs(x(2)-wideb(2)/vals(2)) < 2.0e-10_dp,'wide second')
    call require(maxval(abs(x(3:4)-centers(3:4))) < 2.0e-12_dp,'wide nullspace preserves x0')

    call lsmr_complex(a,b,x,status,maxiter=0,itn=iterations)
    call require(status == 7 .and. iterations == 0,'zero iteration limit')
    call lsmr_complex(a,b,x,status,x0=xtrue,itn=iterations,normr=normr)
    call require(status == 0 .and. iterations == 0,'exact starting solution')
    call require(normr < 1.0e-13_dp,'exact initial residual')
    call lsmr_complex(a,b,x,status,damp=-1.0_dp)
    call require(status == complex_lsmr_invalid_input .and. size(x) == 0,'negative damping')
    call lsmr_complex(a,b,x,status,atol=-1.0_dp)
    call require(status == complex_lsmr_invalid_input,'negative tolerance')
    call lsmr_complex(a,b,x,status,maxiter=-1)
    call require(status == complex_lsmr_invalid_input,'negative iterations')
    call lsmr_complex(a,b(:3),x,status)
    call require(status == complex_lsmr_invalid_input,'invalid RHS length')
    call lsmr_complex(a,b,x,status,x0=xtrue(:2))
    call require(status == complex_lsmr_invalid_input,'invalid x0')
    b(2) = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
    call lsmr_complex(a,b,x,status)
    call require(status == complex_lsmr_invalid_input,'nonfinite RHS')
    ! Large sparse tall system with no dense allocation (6,000 by 120).
    deallocate(rows,cols,vals)
    allocate(rows(3*6000),cols(3*6000),vals(3*6000))
    do i = 1,6000
        rows(3*i-2:3*i) = i
        cols(3*i-2) = mod(i-1,120)+1
        cols(3*i-1) = mod(7*i+5,120)+1
        cols(3*i) = mod(13*i+17,120)+1
        vals(3*i-2) = cmplx(2.0_dp,0.5_dp,dp)
        vals(3*i-1) = cmplx(-0.4_dp,0.25_dp,dp)
        vals(3*i) = cmplx(0.2_dp,-0.45_dp,dp)
    end do
    xtrue = [(cmplx(sin(real(j,dp)),cos(real(j,dp))/2.0_dp,dp), j=1,120)]
    large = csr_complex_from_coo(6000,120,rows,cols,vals,info)
    call require(info == 0,'large COO to CSR')
    rhs = csr_complex_matvec(large,xtrue,info)
    call require(info == 0,'large sparse matvec')
    call lsmr_complex(large,rhs,x,status,atol=1.0e-11_dp,btol=1.0e-11_dp, &
        maxiter=720,normr=normr,itn=iterations)
    call require(status == 1 .or. status == 2,'large sparse convergence')
    call require(iterations > 0 .and. iterations < 720,'large sparse iteration bound')
    call require(maxval(abs(x-xtrue)) < 5.0e-8_dp,'large sparse solution')
    call require(normr < 1.0e-7_dp,'large sparse true residual')

    print *, 'PASS complex CSR/CSC LSMR tests'
contains
    subroutine require(condition,message)
        logical, intent(in) :: condition !! Whether expected condition holds.
        character(len=*), intent(in) :: message !! Failure diagnostic.
        if (.not. condition) then
            print *, 'FAIL: ',message
            error stop 1
        end if
    end subroutine require
end program test_sparse_complex_lsmr
