! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
module scifort_rng_c_bridge
    use, intrinsic :: iso_c_binding, only : c_int, c_int64_t, c_double
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_random, only : rng_state, rng_seed, rng_set_state, rng_get_state, rng_fill_uniform
    implicit none
    private
    public :: rng_load_c_state, rng_store_c_state
    public :: scifort_rng_seed_i64, scifort_rng_uniform_f64
contains
    subroutine rng_load_c_state(values, state, status)
        !! Decode and validate the six-component MRG32k3a state; no stream is advanced.
        implicit none
        integer(c_int64_t), intent(in) :: values(6) !! serialized recurrence components, in native order
        type(rng_state), intent(out) :: state !! decoded native state; unusable when status is nonzero
        integer, intent(out) :: status !! 0 valid; 1 invalid serialized components
        call rng_set_state(state,int(values,int64),status)
    end subroutine rng_load_c_state

    subroutine rng_store_c_state(state, values)
        !! Export all six native recurrence components without advancing the stream.
        implicit none
        type(rng_state), intent(in) :: state !! valid native state to serialize
        integer(c_int64_t), intent(out) :: values(6) !! recurrence components in native order
        integer(int64) :: native_values(6)
        call rng_get_state(state,native_values)
        values=int(native_values,c_int64_t)
    end subroutine rng_store_c_state

    subroutine scifort_rng_seed_i64(seed, values, status) bind(c, name='scifort_rng_seed_i64')
        !! Initialize caller-owned serialized state from any signed 64-bit seed.
        implicit none
        integer(c_int64_t), value, intent(in) :: seed !! deterministic seed; full signed int64 range
        integer(c_int64_t), intent(out) :: values(6) !! initialized MRG32k3a recurrence components
        integer(c_int), intent(out) :: status !! zero on success
        type(rng_state) :: state
        integer, parameter :: kind_check=1/merge(1,0,int64==c_int64_t)
        call rng_seed(state,int(seed,int64))
        call rng_store_c_state(state,values)
        status=0_c_int
    end subroutine scifort_rng_seed_i64

    subroutine scifort_rng_uniform_f64(n, values, output, status) &
            bind(c, name='scifort_rng_uniform_f64')
        !! Fill caller memory with native uniforms and advance explicit serialized state.
        !! n >= 0; empty output validates state but does not advance it.
        implicit none
        integer(c_int), value, intent(in) :: n !! requested number of draws, >= 0
        integer(c_int64_t), intent(inout) :: values(6) !! caller-owned recurrence state; valid on entry
        real(c_double), intent(inout) :: output(n) !! shape (n); uniforms in the open interval (0,1)
        integer(c_int), intent(out) :: status !! 0 success; 1 invalid extent; 3 invalid state
        type(rng_state) :: state
        integer :: info
        if (n<0_c_int) then
            status=1_c_int
            return
        end if
        call rng_load_c_state(values,state,info)
        if (info/=0) then
            status=3_c_int
            return
        end if
        call rng_fill_uniform(state,output)
        call rng_store_c_state(state,values)
        status=0_c_int
    end subroutine scifort_rng_uniform_f64
end module scifort_rng_c_bridge
