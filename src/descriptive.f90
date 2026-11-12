! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive
    use scifort_descriptive_lowrank
    use scifort_descriptive_highrank_means
    use scifort_descriptive_highrank_weights_4_9
    use scifort_descriptive_highrank_weights_10_15
    use scifort_descriptive_highrank_stats_5_9
    use scifort_descriptive_highrank_stats_10_15
    use scifort_descriptive_named_scale_lowrank
    use scifort_descriptive_named_scale_highrank
    use scifort_descriptive_standardize_1_9
    use scifort_descriptive_standardize_10_15
    use scifort_descriptive_trimmed_1_8
    use scifort_descriptive_trimmed_9_15
    use scifort_descriptive_proportion_trim_1_8
    use scifort_descriptive_proportion_trim_9_15
    use scifort_descriptive_kstat_1_8
    use scifort_descriptive_kstat_9_15
    use scifort_descriptive_circular_1_8
    use scifort_descriptive_circular_9_15
    use scifort_descriptive_frequency_engine, only : cumfreq, cumfreq_result, relfreq, relfreq_result
    use scifort_descriptive_frequency_1_8
    use scifort_descriptive_frequency_9_15
    use scifort_descriptive_rank_1_8
    use scifort_descriptive_rank_9_15
    use scifort_descriptive_entropy_1_8
    use scifort_descriptive_entropy_9_15
    use scifort_descriptive_summary_engine, only : describe, describe_all, mode, mode_all, mode_keepdims
    use scifort_descriptive_expectile_sigmaclip, only : expectile, sigmaclip
    use scifort_descriptive_lmoment_1_8
    use scifort_descriptive_lmoment_9_15
    implicit none
end module scifort_descriptive
