! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_special_reductions
    use scifort_special_reductions_low, only : log_softmax, logsumexp, &
        logsumexp_keepdims, logsumexp_signed, logsumexp_signed_keepdims, softmax
    use scifort_special_reductions_high_rank, only : log_softmax, logsumexp, &
        logsumexp_keepdims, logsumexp_signed, logsumexp_signed_keepdims, softmax
    use scifort_special_reductions_high_rank_weighted_6_10, only : logsumexp, &
        logsumexp_keepdims, logsumexp_signed, logsumexp_signed_keepdims
    use scifort_special_reductions_high_rank_weighted_11_15, only : logsumexp, &
        logsumexp_keepdims, logsumexp_signed, logsumexp_signed_keepdims
    implicit none
    private

    public :: log_softmax
    public :: logsumexp
    public :: logsumexp_keepdims
    public :: logsumexp_signed
    public :: logsumexp_signed_keepdims
    public :: softmax
end module scifort_special_reductions
