program example

use scifort_kinds, only : dp
use scifort_sparse_complex, only : csr_complex_matrix, csr_complex_from_dense
use scifort_sparse_linalg_complex_lsmr, only : lsmr_complex
implicit none
type(csr_complex_matrix) :: a
complex(dp), allocatable :: x(:)
complex(dp) :: b(3)
integer :: istop, itn
real(dp) :: normr, normar

a = csr_complex_from_dense(reshape([ &
    cmplx(1.0_dp,0.0_dp,dp), cmplx(0.0_dp,0.0_dp,dp), cmplx(2.0_dp,1.0_dp,dp), &
    cmplx(0.0_dp,0.0_dp,dp), cmplx(3.0_dp,0.0_dp,dp), cmplx(0.0_dp,0.0_dp,dp) &
], [3,2]))
b = [cmplx(1.0_dp,1.0_dp,dp), cmplx(2.0_dp,0.0_dp,dp), cmplx(1.0_dp,-1.0_dp,dp)]
call lsmr_complex(a,b,x,istop,itn=itn,normr=normr,normar=normar, &
    atol=1.0e-12_dp,btol=1.0e-12_dp,maxiter=100)

if (istop < 0) error stop 1
print *, "PASS LSMR docs", istop
end program example
