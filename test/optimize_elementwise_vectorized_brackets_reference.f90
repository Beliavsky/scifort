! SPDX-License-Identifier: MIT
! Regenerate: python tools/generate_optimize_elementwise_vectorized_brackets_reference.py
! SciPy 1.17.0 scipy.optimize._elementwise bracket_root/bracket_minimum
module optimize_elementwise_vectorized_brackets_reference
    implicit none
    private
    public :: vector_bracket_root_status, vector_bracket_root_nit, vector_bracket_root_nfev
    public :: vector_bracket_min_status, vector_bracket_min_nit, vector_bracket_min_nfev
    integer, parameter :: vector_bracket_root_status(7) = [0, 0, 0, 0, 0, -2, 0]
    integer, parameter :: vector_bracket_root_nit(7) = [1, 1, 1, 1, 1, 30, 0]
    integer, parameter :: vector_bracket_root_nfev(7) = [4, 4, 4, 4, 4, 62, 2]
    integer, parameter :: vector_bracket_min_status(7) = [0, 0, 0, 0, 0, 0, -2]
    integer, parameter :: vector_bracket_min_nit(7) = [1, 0, 1, 2, 2, 2, 25]
    integer, parameter :: vector_bracket_min_nfev(7) = [4, 3, 4, 5, 5, 5, 28]
end module optimize_elementwise_vectorized_brackets_reference
