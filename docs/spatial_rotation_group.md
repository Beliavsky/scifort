# SciFort finite SO(3) symmetry groups and rotation reduction

`scifort_spatial_rotation_group` adds native Fortran implementations of
[`scipy.spatial.transform.Rotation.create_group`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.create_group.html)
and [`Rotation.reduce`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.reduce.html).
No additional libraries, BLAS or LAPACK are required. It builds on the
existing scalar `rotation_type` and one-dimensional `rotation_stack_type`.

## Construct finite rotational symmetry groups

```fortran
use scifort_spatial_rotation_stack, only : rotation_stack_type
use scifort_spatial_rotation_group, only : rotation_create_group

type(rotation_stack_type) :: cubic, axial
integer :: status

cubic = rotation_create_group('O', status=status)     ! 24 cubic rotations
axial = rotation_create_group('D7', axis='X')         ! 14 rotations, X axis
if (status /= 0) error stop 'Invalid rotation group'
```

| `group` | Order | Mathematical group |
| --- | ---: | --- |
| `'Cn'` | `n` | Cyclic group of rotations about the specified axis |
| `'Dn'` | `2n` | Proper dihedral rotation group (no reflections) |
| `'T'` | 12 | Tetrahedral rotation group |
| `'O'` | 24 | Octahedral/cubic rotation group |
| `'I'` | 60 | Icosahedral/dodecahedral rotation group |

Here `n` is a positive decimal integer, currently capped at 10,000 to bound
resource use. `axis` is `'X'`, `'Y'` or `'Z'` (case insensitive), defaulting
to `'Z'`; as in SciPy, the axis is ignored for `'T'`, `'O'` and `'I'`.
Group names themselves are case sensitive. The algorithms return a
`rotation_stack_type` whose `%at(i)` and `%as_quat()` expose the members.

`'T'` and `'O'` use the standard unit quaternion descriptions of their
binary rotation-group covers. `'I'` uses the unit quaternions corresponding
to the 600-cell vertex configuration, identifying `q` with `-q`. No code is
copied from SciPy or an external numerical package.

## Symmetry reduction

Given input rotation `P`, `rotation_reduce(P, left, right)` finds the
minimum-principal-angle representative among

\[
Q = L_i P R_j,
\]

where `L_i` and `R_j` are members of the optional left and right groups.
It returns a `rotation_type` for scalar `P` or `rotation_stack_type` for a
stack, respectively. This mirrors SciPy's active rotation multiplication
convention. The search maximizes the absolute scalar quaternion component,
which minimizes the corresponding principal SO(3) angle.

```fortran
use scifort_kinds, only : dp
use scifort_spatial_transform, only : rotation_type, rotation_from_rotvec
use scifort_spatial_rotation_stack, only : rotation_stack_type
use scifort_spatial_rotation_group, only : rotation_create_group, rotation_reduce

type(rotation_stack_type) :: left, right
type(rotation_type) :: p, reduced
integer :: i_left, i_right, status

left = rotation_create_group('O')
right = rotation_create_group('D5')
p = rotation_from_rotvec([0.2_dp,-0.6_dp,1.4_dp])
reduced = rotation_reduce(p, left, right, &
                          index_left=i_left, index_right=i_right, status=status)
if (status /= 0) error stop 'Invalid symmetry reduction'
```

The optional `index_left` and `index_right` are **one-based indices** into
the *returned SciFort group ordering*, not SciPy's ordering. If the
corresponding symmetry group is omitted, the index is `0` and the identity
rotation is used. For a rotation stack, the optional index outputs must be
allocatable `integer(:)` arrays and are one-dimensional (one per input
rotation). The original rotation is returned unchanged if both groups are
omitted, up to quaternion sign equivalence.

Multiple group-element pairs can achieve the same minimum angle. Different
witness indices or equally optimal quaternions are therefore possible, even
when the minimum angle agrees with SciPy. The group member enumeration is
not guaranteed to coincide with SciPy's ordering.

## Limits and error behavior

- Scalar inputs and **one-dimensional** rotation stacks are supported;
  arbitrary-rank rotation-grid reduction is not yet implemented.
- Invalid group names, nonpositive or excessive orders, invalid axis names,
  undefined rotations, empty/invalid symmetry stacks, or an empty input
  rotation stack return status `spatial_rotation_invalid_input` (`1`).
- Invalid scalar results have `%valid() == .false.`; invalid stack results
  are similarly undefined. Successful calls return status `0`.
- The exhaustive reduction costs `O(N * L * R)` quaternion compositions
  for `N` input rotations, `L` left group elements, and `R` right group
  elements. This is practical for finite groups of order at most 60, but
  may be expensive for large user-requested cyclic groups.

## Regression and independent comparisons

```text
python tools/spatial_reference/run_native_spatial.py --mode both
python tools/spatial_reference/compare_rotation_group_scipy.py --mode both
```

Native regression tests cover group order, distinctness, inverses,
composition closure, all three dihedral axes, invalid inputs, scalar and
stack reduction, exact witness reconstruction and minimum-angle optimality.
The independent SciPy 1.17.0 comparator uses two randomized seeds, all five
group families, one- and two-sided reduction, near-identity/near-pi inputs,
and both compiler configurations. It checks group membership independent of
ordering, symmetry-angle minima, and reconstructed witnesses for each batch.
