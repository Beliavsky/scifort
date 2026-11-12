! Generated from SciPy 1.17 _numdiff.approx_derivative by
! tools/generate_optimize_trust_constr_constraint_rel_step_reference.py
module test_optimize_trust_constr_constraint_rel_step_reference
    use scifort_kinds, only : dp
    implicit none
    real(dp), parameter :: dense_two_jac(2) = [3.79999999999999982e+00_dp, -6.60000000000000053e+00_dp]
    integer, parameter :: dense_two_nfev = 3
    real(dp), parameter :: sparse_two_jac(2) = [3.79999999999999982e+00_dp, -6.60000000000000053e+00_dp]
    integer, parameter :: sparse_two_nfev = 2
    real(dp), parameter :: dense_three_jac(2) = [4.00000000000000000e+00_dp, -6.00000000000000000e+00_dp]
    integer, parameter :: dense_three_nfev = 5
    real(dp), parameter :: sparse_three_jac(2) = [4.00000000000000000e+00_dp, -6.00000000000000000e+00_dp]
    integer, parameter :: sparse_three_nfev = 3
    real(dp), parameter :: bounded_two_jac(2) = [4.20000000000000018e+00_dp, -6.60000000000000053e+00_dp]
    integer, parameter :: bounded_two_nfev = 3
    real(dp), parameter :: bounded_sparse_two_jac(2) = [4.20000000000000018e+00_dp, -6.60000000000000053e+00_dp]
    integer, parameter :: bounded_sparse_two_nfev = 2
    real(dp), parameter :: bounded_sparse_three_jac(2) = [4.00000000000000888e+00_dp, -6.00000000000000178e+00_dp]
    integer, parameter :: bounded_sparse_three_nfev = 3
    real(dp), parameter :: zero_step_two_jac(2) = [1.49011611938476562e-08_dp, -6.60000000000000053e+00_dp]
    integer, parameter :: zero_step_two_nfev = 3
    real(dp), parameter :: small_two_jac(2) = [1.10000000000000009e+00_dp, -1.10000000000000009e+00_dp]
    integer, parameter :: small_two_nfev = 3
    real(dp), parameter :: reverse_two_jac(2) = [4.20000000000000018e+00_dp, -5.39999999999999947e+00_dp]
    integer, parameter :: reverse_two_nfev = 3
end module test_optimize_trust_constr_constraint_rel_step_reference
