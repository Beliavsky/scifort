! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Public scipy.spatial subset: pairwise distances and KD-tree queries.
module scifort_spatial
    use scifort_spatial_distance, only : point_distance, pdist, cdist, squareform, &
        spatial_distance_success, spatial_distance_invalid_input
    use scifort_spatial_hausdorff, only : directed_hausdorff, hausdorff_distance
    use scifort_spatial_slerp, only : geometric_slerp
    use scifort_spatial_voronoi, only : spherical_voronoi, spherical_voronoi_type
    use scifort_spatial_convex_hull, only : convex_hull, convex_hull_type
    use scifort_spatial_halfspace, only : halfspace_intersection, halfspace_intersection_type
    use scifort_spatial_procrustes, only : procrustes, orthogonal_procrustes, &
        spatial_procrustes_no_convergence
    use scifort_interpolate_kdtree, only : kd_tree_type, kd_tree, kdtree_success, kdtree_invalid_input
    implicit none
    private
    public :: point_distance, pdist, cdist, squareform
    public :: directed_hausdorff, hausdorff_distance, geometric_slerp
    public :: procrustes, orthogonal_procrustes, spatial_procrustes_no_convergence
    public :: spherical_voronoi, spherical_voronoi_type
    public :: convex_hull, convex_hull_type
    public :: halfspace_intersection, halfspace_intersection_type
    public :: spatial_distance_success, spatial_distance_invalid_input
    public :: kd_tree_type, kd_tree, kdtree_success, kdtree_invalid_input
end module scifort_spatial
