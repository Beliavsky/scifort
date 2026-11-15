! SPDX-License-Identifier: MIT AND BSD-3-Clause
! Copyright (c) 2026 SciFort contributors (module adapters and error hook).
! Adapted from Reference-LAPACK v3.12.1. See CODE_PROVENANCE.md and
! LICENSES/lapack-BSD-3-Clause.txt for full terms and source hashes.
! Copyright (c) 1992-2023 The University of Tennessee and The University
! of Tennessee Research Foundation. All rights reserved.
! Copyright (c) 2000-2023 The University of California Berkeley. All rights reserved.
! Copyright (c) 2006-2023 The University of Colorado Denver. All rights reserved.
! LAPACK authors' and NAG attribution is retained below.
! Mechanical adaptation: lowercase free form, project dp, explicit interfaces,
! documented dummy intents and no printing/STOP in the private error hook.
module scifort_linalg_qz_lapack
    use scifort_kinds, only: dp
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    implicit none
    private
    public :: dhgeqz,dgghrd,dtgexc,dlag2
    real(dp), parameter :: zero=0.0_dp, half=0.5_dp, one=1.0_dp
    real(dp), parameter :: safmin=real(radix(0.0_dp),dp)** &
        max(minexponent(0.0_dp)-1,1-maxexponent(0.0_dp))
    real(dp), parameter :: safmax=one/safmin
    real(dp), parameter :: tsml=real(radix(0.0_dp),dp)** &
        ceiling((minexponent(0.0_dp)-1)*0.5_dp)
    real(dp), parameter :: tbig=real(radix(0.0_dp),dp)** &
        floor((maxexponent(0.0_dp)-digits(0.0_dp)+1)*0.5_dp)
    real(dp), parameter :: ssml=real(radix(0.0_dp),dp)** &
        (-floor((minexponent(0.0_dp)-digits(0.0_dp))*0.5_dp))
    real(dp), parameter :: sbig=real(radix(0.0_dp),dp)** &
        (-ceiling((maxexponent(0.0_dp)+digits(0.0_dp)-1)*0.5_dp))
contains

! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dhgeqz.f
!> \brief \b dhgeqz
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download dhgeqz + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dhgeqz.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dhgeqz.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dhgeqz.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       subroutine dhgeqz( job, compq, compz, n, ilo, ihi, h, ldh, t, ldt,
!                          alphar, alphai, beta, q, ldq, z, ldz, work,
!                          lwork, info )
!
!       .. scalar arguments ..
!       character          compq, compz, job
!       integer            ihi, ilo, info, ldh, ldq, ldt, ldz, lwork, n
!       ..
!       .. array arguments ..
!       double precision   alphai( * ), alphar( * ), beta( * ),
!      $                   h( ldh, * ), q( ldq, * ), t( ldt, * ),
!      $                   work( * ), z( ldz, * )
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dhgeqz computes the eigenvalues of a real matrix pair (h,t),
!> where h is an upper hessenberg matrix and t is upper triangular,
!> using the double-shift qz method.
!> matrix pairs of this type are produced by the reduction to
!> generalized upper hessenberg form of a real matrix pair (a,b):
!>
!>    a = q1*h*z1**t,  b = q1*t*z1**t,
!>
!> as computed by dgghrd.
!>
!> if job='s', then the hessenberg-triangular pair (h,t) is
!> also reduced to generalized schur form,
!>
!>    h = q*s*z**t,  t = q*p*z**t,
!>
!> where q and z are orthogonal matrices, p is an upper triangular
!> matrix, and s is a quasi-triangular matrix with 1-by-1 and 2-by-2
!> diagonal blocks.
!>
!> the 1-by-1 blocks correspond to real eigenvalues of the matrix pair
!> (h,t) and the 2-by-2 blocks correspond to complex conjugate pairs of
!> eigenvalues.
!>
!> additionally, the 2-by-2 upper triangular diagonal blocks of p
!> corresponding to 2-by-2 blocks of s are reduced to positive diagonal
!> form, i.e., if s(j+1,j) is non-zero, then p(j+1,j) = p(j,j+1) = 0,
!> p(j,j) > 0, and p(j+1,j+1) > 0.
!>
!> optionally, the orthogonal matrix q from the generalized schur
!> factorization may be postmultiplied into an input matrix q1, and the
!> orthogonal matrix z may be postmultiplied into an input matrix z1.
!> if q1 and z1 are the orthogonal matrices from dgghrd that reduced
!> the matrix pair (a,b) to generalized upper hessenberg form, then the
!> output matrices q1*q and z1*z are the orthogonal factors from the
!> generalized schur factorization of (a,b):
!>
!>    a = (q1*q)*s*(z1*z)**t,  b = (q1*q)*p*(z1*z)**t.
!>
!> to avoid overflow, eigenvalues of the matrix pair (h,t) (equivalently,
!> of (a,b)) are computed as a pair of values (alpha,beta), where alpha is
!> complex and beta real.
!> if beta is nonzero, lambda = alpha / beta is an eigenvalue of the
!> generalized nonsymmetric eigenvalue problem (gnep)
!>    a*x = lambda*b*x
!> and if alpha is nonzero, mu = beta / alpha is an eigenvalue of the
!> alternate form of the gnep
!>    mu*a*y = b*y.
!> real eigenvalues can be read directly from the generalized schur
!> form:
!>   alpha = s(i,i), beta = p(i,i).
!>
!> ref: c.b. moler & g.w. stewart, "an algorithm for generalized matrix
!>      eigenvalue problems", siam j. numer. anal., 10(1973),
!>      pp. 241--256.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] job
!> \verbatim
!>          job is character*1
!>          = 'e': compute eigenvalues only;
!>          = 's': compute eigenvalues and the schur form.
!> \endverbatim
!>
!> \param[in] compq
!> \verbatim
!>          compq is character*1
!>          = 'n': left schur vectors (q) are not computed;
!>          = 'i': q is initialized to the unit matrix and the matrix q
!>                 of left schur vectors of (h,t) is returned;
!>          = 'v': q must contain an orthogonal matrix q1 on entry and
!>                 the product q1*q is returned.
!> \endverbatim
!>
!> \param[in] compz
!> \verbatim
!>          compz is character*1
!>          = 'n': right schur vectors (z) are not computed;
!>          = 'i': z is initialized to the unit matrix and the matrix z
!>                 of right schur vectors of (h,t) is returned;
!>          = 'v': z must contain an orthogonal matrix z1 on entry and
!>                 the product z1*z is returned.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is integer
!>          the order of the matrices h, t, q, and z.  n >= 0.
!> \endverbatim
!>
!> \param[in] ilo
!> \verbatim
!>          ilo is integer
!> \endverbatim
!>
!> \param[in] ihi
!> \verbatim
!>          ihi is integer
!>          ilo and ihi mark the rows and columns of h which are in
!>          hessenberg form.  it is assumed that a is already upper
!>          triangular in rows and columns 1:ilo-1 and ihi+1:n.
!>          if n > 0, 1 <= ilo <= ihi <= n; if n = 0, ilo=1 and ihi=0.
!> \endverbatim
!>
!> \param[in,out] h
!> \verbatim
!>          h is double precision array, dimension (ldh, n)
!>          on entry, the n-by-n upper hessenberg matrix h.
!>          on exit, if job = 's', h contains the upper quasi-triangular
!>          matrix s from the generalized schur factorization.
!>          if job = 'e', the diagonal blocks of h match those of s, but
!>          the rest of h is unspecified.
!> \endverbatim
!>
!> \param[in] ldh
!> \verbatim
!>          ldh is integer
!>          the leading dimension of the array h.  ldh >= max( 1, n ).
!> \endverbatim
!>
!> \param[in,out] t
!> \verbatim
!>          t is double precision array, dimension (ldt, n)
!>          on entry, the n-by-n upper triangular matrix t.
!>          on exit, if job = 's', t contains the upper triangular
!>          matrix p from the generalized schur factorization;
!>          2-by-2 diagonal blocks of p corresponding to 2-by-2 blocks of s
!>          are reduced to positive diagonal form, i.e., if h(j+1,j) is
!>          non-zero, then t(j+1,j) = t(j,j+1) = 0, t(j,j) > 0, and
!>          t(j+1,j+1) > 0.
!>          if job = 'e', the diagonal blocks of t match those of p, but
!>          the rest of t is unspecified.
!> \endverbatim
!>
!> \param[in] ldt
!> \verbatim
!>          ldt is integer
!>          the leading dimension of the array t.  ldt >= max( 1, n ).
!> \endverbatim
!>
!> \param[out] alphar
!> \verbatim
!>          alphar is double precision array, dimension (n)
!>          the real parts of each scalar alpha defining an eigenvalue
!>          of gnep.
!> \endverbatim
!>
!> \param[out] alphai
!> \verbatim
!>          alphai is double precision array, dimension (n)
!>          the imaginary parts of each scalar alpha defining an
!>          eigenvalue of gnep.
!>          if alphai(j) is zero, then the j-th eigenvalue is real; if
!>          positive, then the j-th and (j+1)-st eigenvalues are a
!>          complex conjugate pair, with alphai(j+1) = -alphai(j).
!> \endverbatim
!>
!> \param[out] beta
!> \verbatim
!>          beta is double precision array, dimension (n)
!>          the scalars beta that define the eigenvalues of gnep.
!>          together, the quantities alpha = (alphar(j),alphai(j)) and
!>          beta = beta(j) represent the j-th eigenvalue of the matrix
!>          pair (a,b), in one of the forms lambda = alpha/beta or
!>          mu = beta/alpha.  since either lambda or mu may overflow,
!>          they should not, in general, be computed.
!> \endverbatim
!>
!> \param[in,out] q
!> \verbatim
!>          q is double precision array, dimension (ldq, n)
!>          on entry, if compq = 'v', the orthogonal matrix q1 used in
!>          the reduction of (a,b) to generalized hessenberg form.
!>          on exit, if compq = 'i', the orthogonal matrix of left schur
!>          vectors of (h,t), and if compq = 'v', the orthogonal matrix
!>          of left schur vectors of (a,b).
!>          not referenced if compq = 'n'.
!> \endverbatim
!>
!> \param[in] ldq
!> \verbatim
!>          ldq is integer
!>          the leading dimension of the array q.  ldq >= 1.
!>          if compq='v' or 'i', then ldq >= n.
!> \endverbatim
!>
!> \param[in,out] z
!> \verbatim
!>          z is double precision array, dimension (ldz, n)
!>          on entry, if compz = 'v', the orthogonal matrix z1 used in
!>          the reduction of (a,b) to generalized hessenberg form.
!>          on exit, if compz = 'i', the orthogonal matrix of
!>          right schur vectors of (h,t), and if compz = 'v', the
!>          orthogonal matrix of right schur vectors of (a,b).
!>          not referenced if compz = 'n'.
!> \endverbatim
!>
!> \param[in] ldz
!> \verbatim
!>          ldz is integer
!>          the leading dimension of the array z.  ldz >= 1.
!>          if compz='v' or 'i', then ldz >= n.
!> \endverbatim
!>
!> \param[out] work
!> \verbatim
!>          work is double precision array, dimension (max(1,lwork))
!>          on exit, if info >= 0, work(1) returns the optimal lwork.
!> \endverbatim
!>
!> \param[in] lwork
!> \verbatim
!>          lwork is integer
!>          the dimension of the array work.  lwork >= max(1,n).
!>
!>          if lwork = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the work array, returns
!>          this value as the first entry of the work array, and no error
!>          message related to lwork is issued by xerbla.
!> \endverbatim
!>
!> \param[out] info
!> \verbatim
!>          info is integer
!>          = 0: successful exit
!>          < 0: if info = -i, the i-th argument had an illegal value
!>          = 1,...,n: the qz iteration did not converge.  (h,t) is not
!>                     in schur form, but alphar(i), alphai(i), and
!>                     beta(i), i=info+1,...,n should be correct.
!>          = n+1,...,2*n: the shift calculation failed.  (h,t) is not
!>                     in schur form, but alphar(i), alphai(i), and
!>                     beta(i), i=info-n+1,...,n should be correct.
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup hgeqz
!
!> \par further details:
!  =====================
!>
!> \verbatim
!>
!>  iteration counters:
!>
!>  jiter  -- counts iterations.
!>  iiter  -- counts iterations run since ilast was last
!>            changed.  this is therefore reset only when a 1-by-1 or
!>            2-by-2 block deflates off the bottom.
!> \endverbatim
!>
!  =====================================================================
subroutine dhgeqz( job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, alphar, alphai, beta, q, ldq, z, ldz, &
    work, lwork, info )
implicit none
!
!  -- lapack computational routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
character, intent(in) :: compq !! accumulate left vectors: n none, i initialize, v update input basis.
character, intent(in) :: compz !! accumulate right vectors: n none, i initialize, v update input basis.
character, intent(in) :: job !! select eigenvalues only (e) or schur factors and eigenvalues (s).
integer, intent(in) :: ihi !! ilo and ihi mark the rows and columns of h which are in hessenberg
!! form.  it is assumed that a is already upper
integer, intent(in) :: ilo !! first row of the active hessenberg block; one for the full pencil.
integer, intent(out) :: info !! zero success; negative invalid argument or positive numerical failure.
integer, intent(in) :: ldh !! the leading dimension of the array h.  ldh >= max( 1, n ).
integer, intent(in) :: ldq !! the leading dimension of the array q.  ldq >= 1. if compq='v' or
!! 'i', then ldq >= n.
integer, intent(in) :: ldt !! the leading dimension of the array t.  ldt >= max( 1, n ).
integer, intent(in) :: ldz !! the leading dimension of the array z.  ldz >= 1. if compz='v' or
!! 'i', then ldz >= n.
integer, intent(in) :: lwork !! the dimension of the array work.  lwork >= max(1,n). if lwork =
!! -1, then a workspace query is assumed; the routine
integer, intent(in) :: n !! the order of the matrices h, t, q, and z.  n >= 0.
!     ..
!     .. array arguments ..
real(dp), intent(out) :: alphai( * ) !! the imaginary parts of each scalar alpha defining an
!! eigenvalue of gnep. if alphai(j) is zero, then the j-th
!! eigenvalue is real; if
real(dp), intent(out) :: alphar( * ) !! the real parts of each scalar alpha defining an eigenvalue
!! of gnep.
real(dp), intent(out) :: beta( * ) !! the scalars beta that define the eigenvalues of gnep.
!! together, the quantities alpha = (alphar(j),alphai(j)) and
real(dp), intent(inout) :: h( ldh, * ) !! on entry, the n-by-n upper hessenberg matrix h. on exit,
!! if job = 's', h contains the upper quasi-triangular
real(dp), intent(inout) :: q( ldq, * ) !! on entry, if compq = 'v', the orthogonal matrix q1 used
!! in the reduction of (a,b) to generalized hessenberg
!! form.
real(dp), intent(inout) :: t( ldt, * ) !! on entry, the n-by-n upper triangular matrix t. on exit,
!! if job = 's', t contains the upper triangular
real(dp), intent(out) :: work( * ) !! on exit, if info >= 0, work(1) returns the optimal lwork.
real(dp), intent(inout) :: z( ldz, * ) !! on entry, if compz = 'v', the orthogonal matrix z1 used
!! in the reduction of (a,b) to generalized hessenberg
!! form.
!     ..
!
!  =====================================================================
!
!     .. parameters ..
!    $                     safety = 1.0e+0 )
real(dp) :: half, zero, one, safety
parameter          ( half = 0.5d+0, zero = 0.0d+0, one = 1.0d+0, safety = 1.0d+2 )
!     ..
!     .. local scalars ..
logical :: ilazr2, ilazro, ilpivt, ilq, ilschr, ilz, lquery
integer :: icompq, icompz, ifirst, ifrstm, iiter, ilast, ilastm, in, ischur, istart, j, jc, jch, jiter, &
    jr, maxit
real(dp) :: a11, a12, a1i, a1r, a21, a22, a2i, a2r, ad11, ad11l, ad12, ad12l, ad21, ad21l, ad22, ad22l, &
    ad32l, an, anorm, ascale, atol, b11, b1a, b1i, b1r, b22, b2a, b2i, b2r, bn, bnorm, bscale, btol, c, &
    c11i, c11r, c12, c21, c22i, c22r, cl, cq, cr, cz, eshift, s, s1, s1inv, s2, safmax, safmin, scale, sl, &
    sqi, sqr, sr, szi, szr, t1, t2, t3, tau, temp, temp2, tempi, tempr, u1, u12, u12l, u2, ulp, vs, w11, &
    w12, w21, w22, wabs, wi, wr, wr2
!     ..
!     .. local arrays ..
real(dp) :: v( 3 )
!     ..
!     .. external functions ..
!     ..
!     .. external subroutines ..
!     ..
!     .. intrinsic functions ..
intrinsic          abs, dble, max, min, sqrt
!     ..
!     .. executable statements ..
!
!     decode job, compq, compz
!
if( lsame( job, 'e' ) ) then
ilschr = .false.
ischur = 1
else if( lsame( job, 's' ) ) then
ilschr = .true.
ischur = 2
else
ischur = 0
end if
!
if( lsame( compq, 'n' ) ) then
ilq = .false.
icompq = 1
else if( lsame( compq, 'v' ) ) then
ilq = .true.
icompq = 2
else if( lsame( compq, 'i' ) ) then
ilq = .true.
icompq = 3
else
icompq = 0
end if
!
if( lsame( compz, 'n' ) ) then
ilz = .false.
icompz = 1
else if( lsame( compz, 'v' ) ) then
ilz = .true.
icompz = 2
else if( lsame( compz, 'i' ) ) then
ilz = .true.
icompz = 3
else
icompz = 0
end if
!
!     check argument values
!
info = 0
work( 1 ) = max( 1, n )
lquery = ( lwork.eq.-1 )
if( ischur.eq.0 ) then
info = -1
else if( icompq.eq.0 ) then
info = -2
else if( icompz.eq.0 ) then
info = -3
else if( n.lt.0 ) then
info = -4
else if( ilo.lt.1 ) then
info = -5
else if( ihi.gt.n .or. ihi.lt.ilo-1 ) then
info = -6
else if( ldh.lt.n ) then
info = -8
else if( ldt.lt.n ) then
info = -10
else if( ldq.lt.1 .or. ( ilq .and. ldq.lt.n ) ) then
info = -15
else if( ldz.lt.1 .or. ( ilz .and. ldz.lt.n ) ) then
info = -17
else if( lwork.lt.max( 1, n ) .and. .not.lquery ) then
info = -19
end if
if( info.ne.0 ) then
call xerbla( 'dhgeqz', -info )
return
else if( lquery ) then
return
end if
!
!     quick return if possible
!
if( n.le.0 ) then
work( 1 ) = dble( 1 )
return
end if
!
!     initialize q and z
!
if( icompq.eq.3 ) call dlaset( 'full', n, n, zero, one, q, ldq )
if( icompz.eq.3 ) call dlaset( 'full', n, n, zero, one, z, ldz )
!
!     machine constants
!
in = ihi + 1 - ilo
safmin = dlamch( 's' )
safmax = one / safmin
ulp = dlamch( 'e' )*dlamch( 'b' )
anorm = dlanhs( 'f', in, h( ilo, ilo ), ldh, work )
bnorm = dlanhs( 'f', in, t( ilo, ilo ), ldt, work )
atol = max( safmin, ulp*anorm )
btol = max( safmin, ulp*bnorm )
ascale = one / max( safmin, anorm )
bscale = one / max( safmin, bnorm )
!
!     set eigenvalues ihi+1:n
!
do 30 j = ihi + 1, n
if( t( j, j ).lt.zero ) then
if( ilschr ) then
do 10 jr = 1, j
h( jr, j ) = -h( jr, j )
t( jr, j ) = -t( jr, j )
10 continue
else
h( j, j ) = -h( j, j )
t( j, j ) = -t( j, j )
end if
if( ilz ) then
do 20 jr = 1, n
z( jr, j ) = -z( jr, j )
20 continue
end if
end if
alphar( j ) = h( j, j )
alphai( j ) = zero
beta( j ) = t( j, j )
30 continue
!
!     if ihi < ilo, skip qz steps
!
if( ihi.lt.ilo ) go to 380
!
!     main qz iteration loop
!
!     initialize dynamic indices
!
!     eigenvalues ilast+1:n have been found.
!        column operations modify rows ifrstm:whatever.
!        row operations modify columns whatever:ilastm.
!
!     if only eigenvalues are being computed, then
!        ifrstm is the row of the last splitting row above row ilast;
!        this is always at least ilo.
!     iiter counts iterations since the last eigenvalue was found,
!        to tell when to use an extraordinary shift.
!     maxit is the maximum number of qz sweeps allowed.
!
ilast = ihi
if( ilschr ) then
ifrstm = 1
ilastm = n
else
ifrstm = ilo
ilastm = ihi
end if
iiter = 0
eshift = zero
maxit = 30*( ihi-ilo+1 )
!
do 360 jiter = 1, maxit
!
!        split the matrix if possible.
!
!        two tests:
!           1: h(j,j-1)=0  or  j=ilo
!           2: t(j,j)=0
!
if( ilast.eq.ilo ) then
!
!           special case: j=ilast
!
go to 80
else
if( abs( h( ilast, ilast-1 ) ).le.max( safmin, ulp*( abs( h( ilast, ilast ) ) + abs( h( ilast-1, ilast-1 &
    ) ) ) ) ) then
h( ilast, ilast-1 ) = zero
go to 80
end if
end if
!
if( abs( t( ilast, ilast ) ).le.btol ) then
t( ilast, ilast ) = zero
go to 70
end if
!
!        general case: j<ilast
!
do 60 j = ilast - 1, ilo, -1
!
!           test 1: for h(j,j-1)=0 or j=ilo
!
if( j.eq.ilo ) then
ilazro = .true.
else
if( abs( h( j, j-1 ) ).le.max( safmin, ulp*( abs( h( j, j ) ) + abs( h( j-1, j-1 ) ) ) ) ) then
h( j, j-1 ) = zero
ilazro = .true.
else
ilazro = .false.
end if
end if
!
!           test 2: for t(j,j)=0
!
if( abs( t( j, j ) ).lt.btol ) then
t( j, j ) = zero
!
!              test 1a: check for 2 consecutive small subdiagonals in a
!
ilazr2 = .false.
if( .not.ilazro ) then
temp = abs( h( j, j-1 ) )
temp2 = abs( h( j, j ) )
tempr = max( temp, temp2 )
if( tempr.lt.one .and. tempr.ne.zero ) then
temp = temp / tempr
temp2 = temp2 / tempr
end if
if( temp*( ascale*abs( h( j+1, j ) ) ).le.temp2* ( ascale*atol ) )ilazr2 = .true.
end if
!
!              if both tests pass (1 & 2), i.e., the leading diagonal
!              element of b in the block is zero, split a 1x1 block off
!              at the top. (i.e., at the j-th row/column) the leading
!              diagonal element of the remainder can also be zero, so
!              this may have to be done repeatedly.
!
if( ilazro .or. ilazr2 ) then
do 40 jch = j, ilast - 1
temp = h( jch, jch )
call dlartg( temp, h( jch+1, jch ), c, s, h( jch, jch ) )
h( jch+1, jch ) = zero
call drot( ilastm-jch, h( jch, jch+1 ), ldh, h( jch+1, jch+1 ), ldh, c, s )
call drot( ilastm-jch, t( jch, jch+1 ), ldt, t( jch+1, jch+1 ), ldt, c, s )
if( ilq ) call drot( n, q( 1, jch ), 1, q( 1, jch+1 ), 1, c, s )
if( ilazr2 ) h( jch, jch-1 ) = h( jch, jch-1 )*c
ilazr2 = .false.
if( abs( t( jch+1, jch+1 ) ).ge.btol ) then
if( jch+1.ge.ilast ) then
go to 80
else
ifirst = jch + 1
go to 110
end if
end if
t( jch+1, jch+1 ) = zero
40 continue
go to 70
else
!
!                 only test 2 passed -- chase the zero to t(ilast,ilast)
!                 then process as in the case t(ilast,ilast)=0
!
do 50 jch = j, ilast - 1
temp = t( jch, jch+1 )
call dlartg( temp, t( jch+1, jch+1 ), c, s, t( jch, jch+1 ) )
t( jch+1, jch+1 ) = zero
if( jch.lt.ilastm-1 ) call drot( ilastm-jch-1, t( jch, jch+2 ), ldt, t( jch+1, jch+2 ), ldt, c, s )
call drot( ilastm-jch+2, h( jch, jch-1 ), ldh, h( jch+1, jch-1 ), ldh, c, s )
if( ilq ) call drot( n, q( 1, jch ), 1, q( 1, jch+1 ), 1, c, s )
temp = h( jch+1, jch )
call dlartg( temp, h( jch+1, jch-1 ), c, s, h( jch+1, jch ) )
h( jch+1, jch-1 ) = zero
call drot( jch+1-ifrstm, h( ifrstm, jch ), 1, h( ifrstm, jch-1 ), 1, c, s )
call drot( jch-ifrstm, t( ifrstm, jch ), 1, t( ifrstm, jch-1 ), 1, c, s )
if( ilz ) call drot( n, z( 1, jch ), 1, z( 1, jch-1 ), 1, c, s )
50 continue
go to 70
end if
else if( ilazro ) then
!
!              only test 1 passed -- work on j:ilast
!
ifirst = j
go to 110
end if
!
!           neither test passed -- try next j
!
60 continue
!
!        (drop-through is "impossible")
!
info = n + 1
go to 420
!
!        t(ilast,ilast)=0 -- clear h(ilast,ilast-1) to split off a
!        1x1 block.
!
70 continue
temp = h( ilast, ilast )
call dlartg( temp, h( ilast, ilast-1 ), c, s, h( ilast, ilast ) )
h( ilast, ilast-1 ) = zero
call drot( ilast-ifrstm, h( ifrstm, ilast ), 1, h( ifrstm, ilast-1 ), 1, c, s )
call drot( ilast-ifrstm, t( ifrstm, ilast ), 1, t( ifrstm, ilast-1 ), 1, c, s )
if( ilz ) call drot( n, z( 1, ilast ), 1, z( 1, ilast-1 ), 1, c, s )
!
!        h(ilast,ilast-1)=0 -- standardize b, set alphar, alphai,
!                              and beta
!
80 continue
if( t( ilast, ilast ).lt.zero ) then
if( ilschr ) then
do 90 j = ifrstm, ilast
h( j, ilast ) = -h( j, ilast )
t( j, ilast ) = -t( j, ilast )
90 continue
else
h( ilast, ilast ) = -h( ilast, ilast )
t( ilast, ilast ) = -t( ilast, ilast )
end if
if( ilz ) then
do 100 j = 1, n
z( j, ilast ) = -z( j, ilast )
100 continue
end if
end if
alphar( ilast ) = h( ilast, ilast )
alphai( ilast ) = zero
beta( ilast ) = t( ilast, ilast )
!
!        go to next block -- exit if finished.
!
ilast = ilast - 1
if( ilast.lt.ilo ) go to 380
!
!        reset counters
!
iiter = 0
eshift = zero
if( .not.ilschr ) then
ilastm = ilast
if( ifrstm.gt.ilast ) ifrstm = ilo
end if
go to 350
!
!        qz step
!
!        this iteration only involves rows/columns ifirst:ilast. we
!        assume ifirst < ilast, and that the diagonal of b is non-zero.
!
110 continue
iiter = iiter + 1
if( .not.ilschr ) then
ifrstm = ifirst
end if
!
!        compute single shifts.
!
!        at this point, ifirst < ilast, and the diagonal elements of
!        t(ifirst:ilast,ifirst,ilast) are larger than btol (in
!        magnitude)
!
if( ( iiter / 10 )*10.eq.iiter ) then
!
!           exceptional shift.  chosen for no particularly good reason.
!           (single shift only.)
!
if( ( dble( maxit )*safmin )*abs( h( ilast, ilast-1 ) ).lt. abs( t( ilast-1, ilast-1 ) ) ) then
eshift = h( ilast, ilast-1 ) / t( ilast-1, ilast-1 )
else
eshift = eshift + one / ( safmin*dble( maxit ) )
end if
s1 = one
wr = eshift
!
else
!
!           shifts based on the generalized eigenvalues of the
!           bottom-right 2x2 block of a and b. the first eigenvalue
!           returned by dlag2 is the wilkinson shift (aep p.512),
!
call dlag2( h( ilast-1, ilast-1 ), ldh, t( ilast-1, ilast-1 ), ldt, safmin*safety, s1, s2, wr, wr2, wi )
!
if ( abs( (wr/s1)*t( ilast, ilast ) - h( ilast, ilast ) ) .gt. abs( (wr2/s2)*t( ilast, ilast ) - h( ilast, &
    ilast ) ) ) then
temp = wr
wr = wr2
wr2 = temp
temp = s1
s1 = s2
s2 = temp
end if
temp = max( s1, safmin*max( one, abs( wr ), abs( wi ) ) )
if( wi.ne.zero ) go to 200
end if
!
!        fiddle with shift to avoid overflow
!
temp = min( ascale, one )*( half*safmax )
if( s1.gt.temp ) then
scale = temp / s1
else
scale = one
end if
!
temp = min( bscale, one )*( half*safmax )
if( abs( wr ).gt.temp ) scale = min( scale, temp / abs( wr ) )
s1 = scale*s1
wr = scale*wr
!
!        now check for two consecutive small subdiagonals.
!
do 120 j = ilast - 1, ifirst + 1, -1
istart = j
temp = abs( s1*h( j, j-1 ) )
temp2 = abs( s1*h( j, j )-wr*t( j, j ) )
tempr = max( temp, temp2 )
if( tempr.lt.one .and. tempr.ne.zero ) then
temp = temp / tempr
temp2 = temp2 / tempr
end if
if( abs( ( ascale*h( j+1, j ) )*temp ).le.( ascale*atol )* temp2 )go to 130
120 continue
!
istart = ifirst
130 continue
!
!        do an implicit single-shift qz sweep.
!
!        initial q
!
temp = s1*h( istart, istart ) - wr*t( istart, istart )
temp2 = s1*h( istart+1, istart )
call dlartg( temp, temp2, c, s, tempr )
!
!        sweep
!
do 190 j = istart, ilast - 1
if( j.gt.istart ) then
temp = h( j, j-1 )
call dlartg( temp, h( j+1, j-1 ), c, s, h( j, j-1 ) )
h( j+1, j-1 ) = zero
end if
!
do 140 jc = j, ilastm
temp = c*h( j, jc ) + s*h( j+1, jc )
h( j+1, jc ) = -s*h( j, jc ) + c*h( j+1, jc )
h( j, jc ) = temp
temp2 = c*t( j, jc ) + s*t( j+1, jc )
t( j+1, jc ) = -s*t( j, jc ) + c*t( j+1, jc )
t( j, jc ) = temp2
140 continue
if( ilq ) then
do 150 jr = 1, n
temp = c*q( jr, j ) + s*q( jr, j+1 )
q( jr, j+1 ) = -s*q( jr, j ) + c*q( jr, j+1 )
q( jr, j ) = temp
150 continue
end if
!
temp = t( j+1, j+1 )
call dlartg( temp, t( j+1, j ), c, s, t( j+1, j+1 ) )
t( j+1, j ) = zero
!
do 160 jr = ifrstm, min( j+2, ilast )
temp = c*h( jr, j+1 ) + s*h( jr, j )
h( jr, j ) = -s*h( jr, j+1 ) + c*h( jr, j )
h( jr, j+1 ) = temp
160 continue
do 170 jr = ifrstm, j
temp = c*t( jr, j+1 ) + s*t( jr, j )
t( jr, j ) = -s*t( jr, j+1 ) + c*t( jr, j )
t( jr, j+1 ) = temp
170 continue
if( ilz ) then
do 180 jr = 1, n
temp = c*z( jr, j+1 ) + s*z( jr, j )
z( jr, j ) = -s*z( jr, j+1 ) + c*z( jr, j )
z( jr, j+1 ) = temp
180 continue
end if
190 continue
!
go to 350
!
!        use francis double-shift
!
!        note: the francis double-shift should work with real shifts,
!              but only if the block is at least 3x3.
!              this code may break if this point is reached with
!              a 2x2 block with real eigenvalues.
!
200 continue
if( ifirst+1.eq.ilast ) then
!
!           special case -- 2x2 block with complex eigenvectors
!
!           step 1: standardize, that is, rotate so that
!
!                       ( b11  0  )
!                   b = (         )  with b11 non-negative.
!                       (  0  b22 )
!
call dlasv2( t( ilast-1, ilast-1 ), t( ilast-1, ilast ), t( ilast, ilast ), b22, b11, sr, cr, sl, cl )
!
if( b11.lt.zero ) then
cr = -cr
sr = -sr
b11 = -b11
b22 = -b22
end if
!
call drot( ilastm+1-ifirst, h( ilast-1, ilast-1 ), ldh, h( ilast, ilast-1 ), ldh, cl, sl )
call drot( ilast+1-ifrstm, h( ifrstm, ilast-1 ), 1, h( ifrstm, ilast ), 1, cr, sr )
!
if( ilast.lt.ilastm ) call drot( ilastm-ilast, t( ilast-1, ilast+1 ), ldt, t( ilast, ilast+1 ), ldt, cl, sl )
if( ifrstm.lt.ilast-1 ) call drot( ifirst-ifrstm, t( ifrstm, ilast-1 ), 1, t( ifrstm, ilast ), 1, cr, sr )
!
if( ilq ) call drot( n, q( 1, ilast-1 ), 1, q( 1, ilast ), 1, cl, sl )
if( ilz ) call drot( n, z( 1, ilast-1 ), 1, z( 1, ilast ), 1, cr, sr )
!
t( ilast-1, ilast-1 ) = b11
t( ilast-1, ilast ) = zero
t( ilast, ilast-1 ) = zero
t( ilast, ilast ) = b22
!
!           if b22 is negative, negate column ilast
!
if( b22.lt.zero ) then
do 210 j = ifrstm, ilast
h( j, ilast ) = -h( j, ilast )
t( j, ilast ) = -t( j, ilast )
210 continue
!
if( ilz ) then
do 220 j = 1, n
z( j, ilast ) = -z( j, ilast )
220 continue
end if
b22 = -b22
end if
!
!           step 2: compute alphar, alphai, and beta (see refs.)
!
!           recompute shift
!
call dlag2( h( ilast-1, ilast-1 ), ldh, t( ilast-1, ilast-1 ), ldt, safmin*safety, s1, temp, wr, temp2, wi )
!
!           if standardization has perturbed the shift onto real line,
!           do another (real single-shift) qr step.
!
if( wi.eq.zero ) go to 350
s1inv = one / s1
!
!           do eispack (qzval) computation of alpha and beta
!
a11 = h( ilast-1, ilast-1 )
a21 = h( ilast, ilast-1 )
a12 = h( ilast-1, ilast )
a22 = h( ilast, ilast )
!
!           compute complex givens rotation on right
!           (assume some element of c = (sa - wb) > unfl )
!                            __
!           (sa - wb) ( cz   -sz )
!                     ( sz    cz )
!
c11r = s1*a11 - wr*b11
c11i = -wi*b11
c12 = s1*a12
c21 = s1*a21
c22r = s1*a22 - wr*b22
c22i = -wi*b22
!
if( abs( c11r )+abs( c11i )+abs( c12 ).gt.abs( c21 )+ abs( c22r )+abs( c22i ) ) then
t1 = dlapy3( c12, c11r, c11i )
cz = c12 / t1
szr = -c11r / t1
szi = -c11i / t1
else
cz = dlapy2( c22r, c22i )
if( cz.le.safmin ) then
cz = zero
szr = one
szi = zero
else
tempr = c22r / cz
tempi = c22i / cz
t1 = dlapy2( cz, c21 )
cz = cz / t1
szr = -c21*tempr / t1
szi = c21*tempi / t1
end if
end if
!
!           compute givens rotation on left
!
!           (  cq   sq )
!           (  __      )  a or b
!           ( -sq   cq )
!
an = abs( a11 ) + abs( a12 ) + abs( a21 ) + abs( a22 )
bn = abs( b11 ) + abs( b22 )
wabs = abs( wr ) + abs( wi )
if( s1*an.gt.wabs*bn ) then
cq = cz*b11
sqr = szr*b22
sqi = -szi*b22
else
a1r = cz*a11 + szr*a12
a1i = szi*a12
a2r = cz*a21 + szr*a22
a2i = szi*a22
cq = dlapy2( a1r, a1i )
if( cq.le.safmin ) then
cq = zero
sqr = one
sqi = zero
else
tempr = a1r / cq
tempi = a1i / cq
sqr = tempr*a2r + tempi*a2i
sqi = tempi*a2r - tempr*a2i
end if
end if
t1 = dlapy3( cq, sqr, sqi )
cq = cq / t1
sqr = sqr / t1
sqi = sqi / t1
!
!           compute diagonal elements of qbz
!
tempr = sqr*szr - sqi*szi
tempi = sqr*szi + sqi*szr
b1r = cq*cz*b11 + tempr*b22
b1i = tempi*b22
b1a = dlapy2( b1r, b1i )
b2r = cq*cz*b22 + tempr*b11
b2i = -tempi*b11
b2a = dlapy2( b2r, b2i )
!
!           normalize so beta > 0, and im( alpha1 ) > 0
!
beta( ilast-1 ) = b1a
beta( ilast ) = b2a
alphar( ilast-1 ) = ( wr*b1a )*s1inv
alphai( ilast-1 ) = ( wi*b1a )*s1inv
alphar( ilast ) = ( wr*b2a )*s1inv
alphai( ilast ) = -( wi*b2a )*s1inv
!
!           step 3: go to next block -- exit if finished.
!
ilast = ifirst - 1
if( ilast.lt.ilo ) go to 380
!
!           reset counters
!
iiter = 0
eshift = zero
if( .not.ilschr ) then
ilastm = ilast
if( ifrstm.gt.ilast ) ifrstm = ilo
end if
go to 350
else
!
!           usual case: 3x3 or larger block, using francis implicit
!                       double-shift
!
!                                    2
!           eigenvalue equation is  w  - c w + d = 0,
!
!                                         -1 2        -1
!           so compute 1st column of  (a b  )  - c a b   + d
!           using the formula in qzit (from eispack)
!
!           we assume that the block is at least 3x3
!
ad11 = ( ascale*h( ilast-1, ilast-1 ) ) / ( bscale*t( ilast-1, ilast-1 ) )
ad21 = ( ascale*h( ilast, ilast-1 ) ) / ( bscale*t( ilast-1, ilast-1 ) )
ad12 = ( ascale*h( ilast-1, ilast ) ) / ( bscale*t( ilast, ilast ) )
ad22 = ( ascale*h( ilast, ilast ) ) / ( bscale*t( ilast, ilast ) )
u12 = t( ilast-1, ilast ) / t( ilast, ilast )
ad11l = ( ascale*h( ifirst, ifirst ) ) / ( bscale*t( ifirst, ifirst ) )
ad21l = ( ascale*h( ifirst+1, ifirst ) ) / ( bscale*t( ifirst, ifirst ) )
ad12l = ( ascale*h( ifirst, ifirst+1 ) ) / ( bscale*t( ifirst+1, ifirst+1 ) )
ad22l = ( ascale*h( ifirst+1, ifirst+1 ) ) / ( bscale*t( ifirst+1, ifirst+1 ) )
ad32l = ( ascale*h( ifirst+2, ifirst+1 ) ) / ( bscale*t( ifirst+1, ifirst+1 ) )
u12l = t( ifirst, ifirst+1 ) / t( ifirst+1, ifirst+1 )
!
v( 1 ) = ( ad11-ad11l )*( ad22-ad11l ) - ad12*ad21 + ad21*u12*ad11l + ( ad12l-ad11l*u12l )*ad21l
v( 2 ) = ( ( ad22l-ad11l )-ad21l*u12l-( ad11-ad11l )- ( ad22-ad11l )+ad21*u12 )*ad21l
v( 3 ) = ad32l*ad21l
!
istart = ifirst
!
call dlarfg( 3, v( 1 ), v( 2 ), 1, tau )
v( 1 ) = one
!
!           sweep
!
do 290 j = istart, ilast - 2
!
!              all but last elements: use 3x3 householder transforms.
!
!              zero (j-1)st column of a
!
if( j.gt.istart ) then
v( 1 ) = h( j, j-1 )
v( 2 ) = h( j+1, j-1 )
v( 3 ) = h( j+2, j-1 )
!
call dlarfg( 3, h( j, j-1 ), v( 2 ), 1, tau )
v( 1 ) = one
h( j+1, j-1 ) = zero
h( j+2, j-1 ) = zero
end if
!
t2 = tau*v( 2 )
t3 = tau*v( 3 )
do 230 jc = j, ilastm
temp = h( j, jc )+v( 2 )*h( j+1, jc )+v( 3 )* h( j+2, jc )
h( j, jc ) = h( j, jc ) - temp*tau
h( j+1, jc ) = h( j+1, jc ) - temp*t2
h( j+2, jc ) = h( j+2, jc ) - temp*t3
temp2 = t( j, jc )+v( 2 )*t( j+1, jc )+v( 3 )* t( j+2, jc )
t( j, jc ) = t( j, jc ) - temp2*tau
t( j+1, jc ) = t( j+1, jc ) - temp2*t2
t( j+2, jc ) = t( j+2, jc ) - temp2*t3
230 continue
if( ilq ) then
do 240 jr = 1, n
temp = q( jr, j )+v( 2 )*q( jr, j+1 )+v( 3 )* q( jr, j+2 )
q( jr, j ) = q( jr, j ) - temp*tau
q( jr, j+1 ) = q( jr, j+1 ) - temp*t2
q( jr, j+2 ) = q( jr, j+2 ) - temp*t3
240 continue
end if
!
!              zero j-th column of b (see dlagbc for details)
!
!              swap rows to pivot
!
ilpivt = .false.
temp = max( abs( t( j+1, j+1 ) ), abs( t( j+1, j+2 ) ) )
temp2 = max( abs( t( j+2, j+1 ) ), abs( t( j+2, j+2 ) ) )
if( max( temp, temp2 ).lt.safmin ) then
scale = zero
u1 = one
u2 = zero
go to 250
else if( temp.ge.temp2 ) then
w11 = t( j+1, j+1 )
w21 = t( j+2, j+1 )
w12 = t( j+1, j+2 )
w22 = t( j+2, j+2 )
u1 = t( j+1, j )
u2 = t( j+2, j )
else
w21 = t( j+1, j+1 )
w11 = t( j+2, j+1 )
w22 = t( j+1, j+2 )
w12 = t( j+2, j+2 )
u2 = t( j+1, j )
u1 = t( j+2, j )
end if
!
!              swap columns if nec.
!
if( abs( w12 ).gt.abs( w11 ) ) then
ilpivt = .true.
temp = w12
temp2 = w22
w12 = w11
w22 = w21
w11 = temp
w21 = temp2
end if
!
!              lu-factor
!
temp = w21 / w11
u2 = u2 - temp*u1
w22 = w22 - temp*w12
w21 = zero
!
!              compute scale
!
scale = one
if( abs( w22 ).lt.safmin ) then
scale = zero
u2 = one
u1 = -w12 / w11
go to 250
end if
if( abs( w22 ).lt.abs( u2 ) ) scale = abs( w22 / u2 )
if( abs( w11 ).lt.abs( u1 ) ) scale = min( scale, abs( w11 / u1 ) )
!
!              solve
!
u2 = ( scale*u2 ) / w22
u1 = ( scale*u1-w12*u2 ) / w11
!
250 continue
if( ilpivt ) then
temp = u2
u2 = u1
u1 = temp
end if
!
!              compute householder vector
!
t1 = sqrt( scale**2+u1**2+u2**2 )
tau = one + scale / t1
vs = -one / ( scale+t1 )
v( 1 ) = one
v( 2 ) = vs*u1
v( 3 ) = vs*u2
!
!              apply transformations from the right.
!
t2 = tau*v(2)
t3 = tau*v(3)
do 260 jr = ifrstm, min( j+3, ilast )
temp = h( jr, j )+v( 2 )*h( jr, j+1 )+v( 3 )* h( jr, j+2 )
h( jr, j ) = h( jr, j ) - temp*tau
h( jr, j+1 ) = h( jr, j+1 ) - temp*t2
h( jr, j+2 ) = h( jr, j+2 ) - temp*t3
260 continue
do 270 jr = ifrstm, j + 2
temp = t( jr, j )+v( 2 )*t( jr, j+1 )+v( 3 )* t( jr, j+2 )
t( jr, j ) = t( jr, j ) - temp*tau
t( jr, j+1 ) = t( jr, j+1 ) - temp*t2
t( jr, j+2 ) = t( jr, j+2 ) - temp*t3
270 continue
if( ilz ) then
do 280 jr = 1, n
temp = z( jr, j )+v( 2 )*z( jr, j+1 )+v( 3 )* z( jr, j+2 )
z( jr, j ) = z( jr, j ) - temp*tau
z( jr, j+1 ) = z( jr, j+1 ) - temp*t2
z( jr, j+2 ) = z( jr, j+2 ) - temp*t3
280 continue
end if
t( j+1, j ) = zero
t( j+2, j ) = zero
290 continue
!
!           last elements: use givens rotations
!
!           rotations from the left
!
j = ilast - 1
temp = h( j, j-1 )
call dlartg( temp, h( j+1, j-1 ), c, s, h( j, j-1 ) )
h( j+1, j-1 ) = zero
!
do 300 jc = j, ilastm
temp = c*h( j, jc ) + s*h( j+1, jc )
h( j+1, jc ) = -s*h( j, jc ) + c*h( j+1, jc )
h( j, jc ) = temp
temp2 = c*t( j, jc ) + s*t( j+1, jc )
t( j+1, jc ) = -s*t( j, jc ) + c*t( j+1, jc )
t( j, jc ) = temp2
300 continue
if( ilq ) then
do 310 jr = 1, n
temp = c*q( jr, j ) + s*q( jr, j+1 )
q( jr, j+1 ) = -s*q( jr, j ) + c*q( jr, j+1 )
q( jr, j ) = temp
310 continue
end if
!
!           rotations from the right.
!
temp = t( j+1, j+1 )
call dlartg( temp, t( j+1, j ), c, s, t( j+1, j+1 ) )
t( j+1, j ) = zero
!
do 320 jr = ifrstm, ilast
temp = c*h( jr, j+1 ) + s*h( jr, j )
h( jr, j ) = -s*h( jr, j+1 ) + c*h( jr, j )
h( jr, j+1 ) = temp
320 continue
do 330 jr = ifrstm, ilast - 1
temp = c*t( jr, j+1 ) + s*t( jr, j )
t( jr, j ) = -s*t( jr, j+1 ) + c*t( jr, j )
t( jr, j+1 ) = temp
330 continue
if( ilz ) then
do 340 jr = 1, n
temp = c*z( jr, j+1 ) + s*z( jr, j )
z( jr, j ) = -s*z( jr, j+1 ) + c*z( jr, j )
z( jr, j+1 ) = temp
340 continue
end if
!
!           end of double-shift code
!
end if
!
go to 350
!
!        end of iteration loop
!
350 continue
360 continue
!
!     drop-through = non-convergence
!
info = ilast
go to 420
!
!     successful completion of all qz steps
!
380 continue
!
!     set eigenvalues 1:ilo-1
!
do 410 j = 1, ilo - 1
if( t( j, j ).lt.zero ) then
if( ilschr ) then
do 390 jr = 1, j
h( jr, j ) = -h( jr, j )
t( jr, j ) = -t( jr, j )
390 continue
else
h( j, j ) = -h( j, j )
t( j, j ) = -t( j, j )
end if
if( ilz ) then
do 400 jr = 1, n
z( jr, j ) = -z( jr, j )
400 continue
end if
end if
alphar( j ) = h( j, j )
alphai( j ) = zero
beta( j ) = t( j, j )
410 continue
!
!     normal termination
!
info = 0
!
!     exit (other than argument error) -- return optimal workspace size
!
420 continue
work( 1 ) = dble( n )
return
!
!     end of dhgeqz
!
end subroutine dhgeqz
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dgghrd.f
!> \brief \b dgghrd
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download dgghrd + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgghrd.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgghrd.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgghrd.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       subroutine dgghrd( compq, compz, n, ilo, ihi, a, lda, b, ldb, q,
!                          ldq, z, ldz, info )
!
!       .. scalar arguments ..
!       character          compq, compz
!       integer            ihi, ilo, info, lda, ldb, ldq, ldz, n
!       ..
!       .. array arguments ..
!       double precision   a( lda, * ), b( ldb, * ), q( ldq, * ),
!      $                   z( ldz, * )
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dgghrd reduces a pair of real matrices (a,b) to generalized upper
!> hessenberg form using orthogonal transformations, where a is a
!> general matrix and b is upper triangular.  the form of the
!> generalized eigenvalue problem is
!>    a*x = lambda*b*x,
!> and b is typically made upper triangular by computing its qr
!> factorization and moving the orthogonal matrix q to the left side
!> of the equation.
!>
!> this subroutine simultaneously reduces a to a hessenberg matrix h:
!>    q**t*a*z = h
!> and transforms b to another upper triangular matrix t:
!>    q**t*b*z = t
!> in order to reduce the problem to its standard form
!>    h*y = lambda*t*y
!> where y = z**t*x.
!>
!> the orthogonal matrices q and z are determined as products of givens
!> rotations.  they may either be formed explicitly, or they may be
!> postmultiplied into input matrices q1 and z1, so that
!>
!>      q1 * a * z1**t = (q1*q) * h * (z1*z)**t
!>
!>      q1 * b * z1**t = (q1*q) * t * (z1*z)**t
!>
!> if q1 is the orthogonal matrix from the qr factorization of b in the
!> original equation a*x = lambda*b*x, then dgghrd reduces the original
!> problem to generalized hessenberg form.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] compq
!> \verbatim
!>          compq is character*1
!>          = 'n': do not compute q;
!>          = 'i': q is initialized to the unit matrix, and the
!>                 orthogonal matrix q is returned;
!>          = 'v': q must contain an orthogonal matrix q1 on entry,
!>                 and the product q1*q is returned.
!> \endverbatim
!>
!> \param[in] compz
!> \verbatim
!>          compz is character*1
!>          = 'n': do not compute z;
!>          = 'i': z is initialized to the unit matrix, and the
!>                 orthogonal matrix z is returned;
!>          = 'v': z must contain an orthogonal matrix z1 on entry,
!>                 and the product z1*z is returned.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is integer
!>          the order of the matrices a and b.  n >= 0.
!> \endverbatim
!>
!> \param[in] ilo
!> \verbatim
!>          ilo is integer
!> \endverbatim
!>
!> \param[in] ihi
!> \verbatim
!>          ihi is integer
!>
!>          ilo and ihi mark the rows and columns of a which are to be
!>          reduced.  it is assumed that a is already upper triangular
!>          in rows and columns 1:ilo-1 and ihi+1:n.  ilo and ihi are
!>          normally set by a previous call to dggbal; otherwise they
!>          should be set to 1 and n respectively.
!>          1 <= ilo <= ihi <= n, if n > 0; ilo=1 and ihi=0, if n=0.
!> \endverbatim
!>
!> \param[in,out] a
!> \verbatim
!>          a is double precision array, dimension (lda, n)
!>          on entry, the n-by-n general matrix to be reduced.
!>          on exit, the upper triangle and the first subdiagonal of a
!>          are overwritten with the upper hessenberg matrix h, and the
!>          rest is set to zero.
!> \endverbatim
!>
!> \param[in] lda
!> \verbatim
!>          lda is integer
!>          the leading dimension of the array a.  lda >= max(1,n).
!> \endverbatim
!>
!> \param[in,out] b
!> \verbatim
!>          b is double precision array, dimension (ldb, n)
!>          on entry, the n-by-n upper triangular matrix b.
!>          on exit, the upper triangular matrix t = q**t b z.  the
!>          elements below the diagonal are set to zero.
!> \endverbatim
!>
!> \param[in] ldb
!> \verbatim
!>          ldb is integer
!>          the leading dimension of the array b.  ldb >= max(1,n).
!> \endverbatim
!>
!> \param[in,out] q
!> \verbatim
!>          q is double precision array, dimension (ldq, n)
!>          on entry, if compq = 'v', the orthogonal matrix q1,
!>          typically from the qr factorization of b.
!>          on exit, if compq='i', the orthogonal matrix q, and if
!>          compq = 'v', the product q1*q.
!>          not referenced if compq='n'.
!> \endverbatim
!>
!> \param[in] ldq
!> \verbatim
!>          ldq is integer
!>          the leading dimension of the array q.
!>          ldq >= n if compq='v' or 'i'; ldq >= 1 otherwise.
!> \endverbatim
!>
!> \param[in,out] z
!> \verbatim
!>          z is double precision array, dimension (ldz, n)
!>          on entry, if compz = 'v', the orthogonal matrix z1.
!>          on exit, if compz='i', the orthogonal matrix z, and if
!>          compz = 'v', the product z1*z.
!>          not referenced if compz='n'.
!> \endverbatim
!>
!> \param[in] ldz
!> \verbatim
!>          ldz is integer
!>          the leading dimension of the array z.
!>          ldz >= n if compz='v' or 'i'; ldz >= 1 otherwise.
!> \endverbatim
!>
!> \param[out] info
!> \verbatim
!>          info is integer
!>          = 0:  successful exit.
!>          < 0:  if info = -i, the i-th argument had an illegal value.
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup gghrd
!
!> \par further details:
!  =====================
!>
!> \verbatim
!>
!>  this routine reduces a to hessenberg and b to triangular form by
!>  an unblocked reduction, as described in _matrix_computations_,
!>  by golub and van loan (johns hopkins press.)
!> \endverbatim
!>
!  =====================================================================
subroutine dgghrd( compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, info )
implicit none
!
!  -- lapack computational routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
character, intent(in) :: compq !! accumulate left vectors: n none, i initialize, v update input basis.
character, intent(in) :: compz !! accumulate right vectors: n none, i initialize, v update input basis.
integer, intent(in) :: ihi !! ilo and ihi mark the rows and columns of a which are to be reduced.
!! it is assumed that a is already upper triangular
integer, intent(in) :: ilo !! first row of the active hessenberg block; one for the full pencil.
integer, intent(out) :: info !! zero success; negative invalid argument or positive numerical failure.
integer, intent(in) :: lda !! the leading dimension of the array a.  lda >= max(1,n).
integer, intent(in) :: ldb !! the leading dimension of the array b.  ldb >= max(1,n).
integer, intent(in) :: ldq !! the leading dimension of the array q. ldq >= n if compq='v' or 'i';
!! ldq >= 1 otherwise.
integer, intent(in) :: ldz !! the leading dimension of the array z. ldz >= n if compz='v' or 'i';
!! ldz >= 1 otherwise.
integer, intent(in) :: n !! the order of the matrices a and b.  n >= 0.
!     ..
!     .. array arguments ..
real(dp), intent(inout) :: a( lda, * ) !! on entry, the n-by-n general matrix to be reduced. on
!! exit, the upper triangle and the first subdiagonal of a
real(dp), intent(inout) :: b( ldb, * ) !! on entry, the n-by-n upper triangular matrix b. on exit,
!! the upper triangular matrix t = q**t b z.  the
real(dp), intent(inout) :: q( ldq, * ) !! on entry, if compq = 'v', the orthogonal matrix q1,
!! typically from the qr factorization of b.
real(dp), intent(inout) :: z( ldz, * ) !! on entry, if compz = 'v', the orthogonal matrix z1. on
!! exit, if compz='i', the orthogonal matrix z, and if
!     ..
!
!  =====================================================================
!
!     .. parameters ..
real(dp) :: one, zero
parameter          ( one = 1.0d+0, zero = 0.0d+0 )
!     ..
!     .. local scalars ..
logical :: ilq, ilz
integer :: icompq, icompz, jcol, jrow
real(dp) :: c, s, temp
!     ..
!     .. external functions ..
!     ..
!     .. external subroutines ..
!     ..
!     .. intrinsic functions ..
intrinsic          max
!     ..
!     .. executable statements ..
!
!     decode compq
!
if( lsame( compq, 'n' ) ) then
ilq = .false.
icompq = 1
else if( lsame( compq, 'v' ) ) then
ilq = .true.
icompq = 2
else if( lsame( compq, 'i' ) ) then
ilq = .true.
icompq = 3
else
icompq = 0
end if
!
!     decode compz
!
if( lsame( compz, 'n' ) ) then
ilz = .false.
icompz = 1
else if( lsame( compz, 'v' ) ) then
ilz = .true.
icompz = 2
else if( lsame( compz, 'i' ) ) then
ilz = .true.
icompz = 3
else
icompz = 0
end if
!
!     test the input parameters.
!
info = 0
if( icompq.le.0 ) then
info = -1
else if( icompz.le.0 ) then
info = -2
else if( n.lt.0 ) then
info = -3
else if( ilo.lt.1 ) then
info = -4
else if( ihi.gt.n .or. ihi.lt.ilo-1 ) then
info = -5
else if( lda.lt.max( 1, n ) ) then
info = -7
else if( ldb.lt.max( 1, n ) ) then
info = -9
else if( ( ilq .and. ldq.lt.n ) .or. ldq.lt.1 ) then
info = -11
else if( ( ilz .and. ldz.lt.n ) .or. ldz.lt.1 ) then
info = -13
end if
if( info.ne.0 ) then
call xerbla( 'dgghrd', -info )
return
end if
!
!     initialize q and z if desired.
!
if( icompq.eq.3 ) call dlaset( 'full', n, n, zero, one, q, ldq )
if( icompz.eq.3 ) call dlaset( 'full', n, n, zero, one, z, ldz )
!
!     quick return if possible
!
if( n.le.1 ) return
!
!     zero out lower triangle of b
!
do 20 jcol = 1, n - 1
do 10 jrow = jcol + 1, n
b( jrow, jcol ) = zero
10 continue
20 continue
!
!     reduce a and b
!
do 40 jcol = ilo, ihi - 2
!
do 30 jrow = ihi, jcol + 2, -1
!
!           step 1: rotate rows jrow-1, jrow to kill a(jrow,jcol)
!
temp = a( jrow-1, jcol )
call dlartg( temp, a( jrow, jcol ), c, s, a( jrow-1, jcol ) )
a( jrow, jcol ) = zero
call drot( n-jcol, a( jrow-1, jcol+1 ), lda, a( jrow, jcol+1 ), lda, c, s )
call drot( n+2-jrow, b( jrow-1, jrow-1 ), ldb, b( jrow, jrow-1 ), ldb, c, s )
if( ilq ) call drot( n, q( 1, jrow-1 ), 1, q( 1, jrow ), 1, c, s )
!
!           step 2: rotate columns jrow, jrow-1 to kill b(jrow,jrow-1)
!
temp = b( jrow, jrow )
call dlartg( temp, b( jrow, jrow-1 ), c, s, b( jrow, jrow ) )
b( jrow, jrow-1 ) = zero
call drot( ihi, a( 1, jrow ), 1, a( 1, jrow-1 ), 1, c, s )
call drot( jrow-1, b( 1, jrow ), 1, b( 1, jrow-1 ), 1, c, s )
if( ilz ) call drot( n, z( 1, jrow ), 1, z( 1, jrow-1 ), 1, c, s )
30 continue
40 continue
!
return
!
!     end of dgghrd
!
end subroutine dgghrd
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dtgexc.f
!> \brief \b dtgexc
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download dtgexc + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtgexc.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtgexc.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtgexc.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       subroutine dtgexc( wantq, wantz, n, a, lda, b, ldb, q, ldq, z,
!                          ldz, ifst, ilst, work, lwork, info )
!
!       .. scalar arguments ..
!       logical            wantq, wantz
!       integer            ifst, ilst, info, lda, ldb, ldq, ldz, lwork, n
!       ..
!       .. array arguments ..
!       double precision   a( lda, * ), b( ldb, * ), q( ldq, * ),
!      $                   work( * ), z( ldz, * )
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dtgexc reorders the generalized real schur decomposition of a real
!> matrix pair (a,b) using an orthogonal equivalence transformation
!>
!>                (a, b) = q * (a, b) * z**t,
!>
!> so that the diagonal block of (a, b) with row index ifst is moved
!> to row ilst.
!>
!> (a, b) must be in generalized real schur canonical form (as returned
!> by dgges), i.e. a is block upper triangular with 1-by-1 and 2-by-2
!> diagonal blocks. b is upper triangular.
!>
!> optionally, the matrices q and z of generalized schur vectors are
!> updated.
!>
!>        q(in) * a(in) * z(in)**t = q(out) * a(out) * z(out)**t
!>        q(in) * b(in) * z(in)**t = q(out) * b(out) * z(out)**t
!>
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] wantq
!> \verbatim
!>          wantq is logical
!>          .true. : update the left transformation matrix q;
!>          .false.: do not update q.
!> \endverbatim
!>
!> \param[in] wantz
!> \verbatim
!>          wantz is logical
!>          .true. : update the right transformation matrix z;
!>          .false.: do not update z.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is integer
!>          the order of the matrices a and b. n >= 0.
!> \endverbatim
!>
!> \param[in,out] a
!> \verbatim
!>          a is double precision array, dimension (lda,n)
!>          on entry, the matrix a in generalized real schur canonical
!>          form.
!>          on exit, the updated matrix a, again in generalized
!>          real schur canonical form.
!> \endverbatim
!>
!> \param[in] lda
!> \verbatim
!>          lda is integer
!>          the leading dimension of the array a. lda >= max(1,n).
!> \endverbatim
!>
!> \param[in,out] b
!> \verbatim
!>          b is double precision array, dimension (ldb,n)
!>          on entry, the matrix b in generalized real schur canonical
!>          form (a,b).
!>          on exit, the updated matrix b, again in generalized
!>          real schur canonical form (a,b).
!> \endverbatim
!>
!> \param[in] ldb
!> \verbatim
!>          ldb is integer
!>          the leading dimension of the array b. ldb >= max(1,n).
!> \endverbatim
!>
!> \param[in,out] q
!> \verbatim
!>          q is double precision array, dimension (ldq,n)
!>          on entry, if wantq = .true., the orthogonal matrix q.
!>          on exit, the updated matrix q.
!>          if wantq = .false., q is not referenced.
!> \endverbatim
!>
!> \param[in] ldq
!> \verbatim
!>          ldq is integer
!>          the leading dimension of the array q. ldq >= 1.
!>          if wantq = .true., ldq >= n.
!> \endverbatim
!>
!> \param[in,out] z
!> \verbatim
!>          z is double precision array, dimension (ldz,n)
!>          on entry, if wantz = .true., the orthogonal matrix z.
!>          on exit, the updated matrix z.
!>          if wantz = .false., z is not referenced.
!> \endverbatim
!>
!> \param[in] ldz
!> \verbatim
!>          ldz is integer
!>          the leading dimension of the array z. ldz >= 1.
!>          if wantz = .true., ldz >= n.
!> \endverbatim
!>
!> \param[in,out] ifst
!> \verbatim
!>          ifst is integer
!> \endverbatim
!>
!> \param[in,out] ilst
!> \verbatim
!>          ilst is integer
!>          specify the reordering of the diagonal blocks of (a, b).
!>          the block with row index ifst is moved to row ilst, by a
!>          sequence of swapping between adjacent blocks.
!>          on exit, if ifst pointed on entry to the second row of
!>          a 2-by-2 block, it is changed to point to the first row;
!>          ilst always points to the first row of the block in its
!>          final position (which may differ from its input value by
!>          +1 or -1). 1 <= ifst, ilst <= n.
!> \endverbatim
!>
!> \param[out] work
!> \verbatim
!>          work is double precision array, dimension (max(1,lwork))
!>          on exit, if info = 0, work(1) returns the optimal lwork.
!> \endverbatim
!>
!> \param[in] lwork
!> \verbatim
!>          lwork is integer
!>          the dimension of the array work.
!>          lwork >= 1 when n <= 1, otherwise lwork >= 4*n + 16.
!>
!>          if lwork = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the work array, returns
!>          this value as the first entry of the work array, and no error
!>          message related to lwork is issued by xerbla.
!> \endverbatim
!>
!> \param[out] info
!> \verbatim
!>          info is integer
!>           =0:  successful exit.
!>           <0:  if info = -i, the i-th argument had an illegal value.
!>           =1:  the transformed matrix pair (a, b) would be too far
!>                from generalized schur form; the problem is ill-
!>                conditioned. (a, b) may have been partially reordered,
!>                and ilst points to the first row of the current
!>                position of the block being moved.
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup tgexc
!
!> \par contributors:
!  ==================
!>
!>     bo kagstrom and peter poromaa, department of computing science,
!>     umea university, s-901 87 umea, sweden.
!
!> \par references:
!  ================
!>
!> \verbatim
!>
!>  [1] b. kagstrom; a direct method for reordering eigenvalues in the
!>      generalized real schur form of a regular matrix pair (a, b), in
!>      m.s. moonen et al (eds), linear algebra for large scale and
!>      real-time applications, kluwer academic publ. 1993, pp 195-218.
!> \endverbatim
!>
!  =====================================================================
subroutine dtgexc( wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, ifst, ilst, work, lwork, info )
implicit none
!
!  -- lapack computational routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
logical, intent(in) :: wantq !! .true. : update the left transformation matrix q; .false.: do not
!! update q.
logical, intent(in) :: wantz !! .true. : update the right transformation matrix z; .false.: do not
!! update z.
integer, intent(inout) :: ifst !! first row of the block to move; updated to the actual block start.
integer, intent(inout) :: ilst !! specify the reordering of the diagonal blocks of (a, b). the
!! block with row index ifst is moved to row ilst, by a
integer, intent(out) :: info !! zero success; negative invalid argument or positive numerical failure.
integer, intent(in) :: lda !! the leading dimension of the array a. lda >= max(1,n).
integer, intent(in) :: ldb !! the leading dimension of the array b. ldb >= max(1,n).
integer, intent(in) :: ldq !! the leading dimension of the array q. ldq >= 1. if wantq = .true.,
!! ldq >= n.
integer, intent(in) :: ldz !! the leading dimension of the array z. ldz >= 1. if wantz = .true.,
!! ldz >= n.
integer, intent(in) :: lwork !! the dimension of the array work. lwork >= 1 when n <= 1, otherwise
!! lwork >= 4*n + 16. if lwork = -1, then a workspace query is
!! assumed; the routine
integer, intent(in) :: n !! the order of the matrices a and b. n >= 0.
!     ..
!     .. array arguments ..
real(dp), intent(inout) :: a( lda, * ) !! on entry, the matrix a in generalized real schur
!! canonical form. on exit, the updated matrix a, again in
!! generalized
real(dp), intent(inout) :: b( ldb, * ) !! on entry, the matrix b in generalized real schur
!! canonical form (a,b). on exit, the updated matrix b,
!! again in generalized
real(dp), intent(inout) :: q( ldq, * ) !! on entry, if wantq = .true., the orthogonal matrix q. on
!! exit, the updated matrix q. if wantq = .false., q is not
!! referenced.
real(dp), intent(out) :: work( * ) !! on exit, if info = 0, work(1) returns the optimal lwork.
real(dp), intent(inout) :: z( ldz, * ) !! on entry, if wantz = .true., the orthogonal matrix z. on
!! exit, the updated matrix z. if wantz = .false., z is not
!! referenced.
!     ..
!
!  =====================================================================
!
!     .. parameters ..
real(dp) :: zero
parameter          ( zero = 0.0d+0 )
!     ..
!     .. local scalars ..
logical :: lquery
integer :: here, lwmin, nbf, nbl, nbnext
!     ..
!     .. external subroutines ..
!     ..
!     .. intrinsic functions ..
intrinsic          max
!     ..
!     .. executable statements ..
!
!     decode and test input arguments.
!
info = 0
lquery = ( lwork.eq.-1 )
if( n.lt.0 ) then
info = -3
else if( lda.lt.max( 1, n ) ) then
info = -5
else if( ldb.lt.max( 1, n ) ) then
info = -7
else if( ldq.lt.1 .or. wantq .and. ( ldq.lt.max( 1, n ) ) ) then
info = -9
else if( ldz.lt.1 .or. wantz .and. ( ldz.lt.max( 1, n ) ) ) then
info = -11
else if( ifst.lt.1 .or. ifst.gt.n ) then
info = -12
else if( ilst.lt.1 .or. ilst.gt.n ) then
info = -13
end if
!
if( info.eq.0 ) then
if( n.le.1 ) then
lwmin = 1
else
lwmin = 4*n + 16
end if
work(1) = lwmin
!
if (lwork.lt.lwmin .and. .not.lquery) then
info = -15
end if
end if
!
if( info.ne.0 ) then
call xerbla( 'dtgexc', -info )
return
else if( lquery ) then
return
end if
!
!     quick return if possible
!
if( n.le.1 ) return
!
!     determine the first row of the specified block and find out
!     if it is 1-by-1 or 2-by-2.
!
if( ifst.gt.1 ) then
if( a( ifst, ifst-1 ).ne.zero ) ifst = ifst - 1
end if
nbf = 1
if( ifst.lt.n ) then
if( a( ifst+1, ifst ).ne.zero ) nbf = 2
end if
!
!     determine the first row of the final block
!     and find out if it is 1-by-1 or 2-by-2.
!
if( ilst.gt.1 ) then
if( a( ilst, ilst-1 ).ne.zero ) ilst = ilst - 1
end if
nbl = 1
if( ilst.lt.n ) then
if( a( ilst+1, ilst ).ne.zero ) nbl = 2
end if
if( ifst.eq.ilst ) return
!
if( ifst.lt.ilst ) then
!
!        update ilst.
!
if( nbf.eq.2 .and. nbl.eq.1 ) ilst = ilst - 1
if( nbf.eq.1 .and. nbl.eq.2 ) ilst = ilst + 1
!
here = ifst
!
10 continue
!
!        swap with next one below.
!
if( nbf.eq.1 .or. nbf.eq.2 ) then
!
!           current block either 1-by-1 or 2-by-2.
!
nbnext = 1
if( here+nbf+1.le.n ) then
if( a( here+nbf+1, here+nbf ).ne.zero ) nbnext = 2
end if
call dtgex2( wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, here, nbf, nbnext, work, lwork, info )
if( info.ne.0 ) then
ilst = here
return
end if
here = here + nbnext
!
!           test if 2-by-2 block breaks into two 1-by-1 blocks.
!
if( nbf.eq.2 ) then
if( a( here+1, here ).eq.zero ) nbf = 3
end if
!
else
!
!           current block consists of two 1-by-1 blocks, each of which
!           must be swapped individually.
!
nbnext = 1
if( here+3.le.n ) then
if( a( here+3, here+2 ).ne.zero ) nbnext = 2
end if
call dtgex2( wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, here+1, 1, nbnext, work, lwork, info )
if( info.ne.0 ) then
ilst = here
return
end if
if( nbnext.eq.1 ) then
!
!              swap two 1-by-1 blocks.
!
call dtgex2( wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, here, 1, 1, work, lwork, info )
if( info.ne.0 ) then
ilst = here
return
end if
here = here + 1
!
else
!
!              recompute nbnext in case of 2-by-2 split.
!
if( a( here+2, here+1 ).eq.zero ) nbnext = 1
if( nbnext.eq.2 ) then
!
!                 2-by-2 block did not split.
!
call dtgex2( wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, here, 1, nbnext, work, lwork, info )
if( info.ne.0 ) then
ilst = here
return
end if
here = here + 2
else
!
!                 2-by-2 block did split.
!
call dtgex2( wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, here, 1, 1, work, lwork, info )
if( info.ne.0 ) then
ilst = here
return
end if
here = here + 1
call dtgex2( wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, here, 1, 1, work, lwork, info )
if( info.ne.0 ) then
ilst = here
return
end if
here = here + 1
end if
!
end if
end if
if( here.lt.ilst ) go to 10
else
here = ifst
!
20 continue
!
!        swap with next one below.
!
if( nbf.eq.1 .or. nbf.eq.2 ) then
!
!           current block either 1-by-1 or 2-by-2.
!
nbnext = 1
if( here.ge.3 ) then
if( a( here-1, here-2 ).ne.zero ) nbnext = 2
end if
call dtgex2( wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, here-nbnext, nbnext, nbf, work, lwork, info )
if( info.ne.0 ) then
ilst = here
return
end if
here = here - nbnext
!
!           test if 2-by-2 block breaks into two 1-by-1 blocks.
!
if( nbf.eq.2 ) then
if( a( here+1, here ).eq.zero ) nbf = 3
end if
!
else
!
!           current block consists of two 1-by-1 blocks, each of which
!           must be swapped individually.
!
nbnext = 1
if( here.ge.3 ) then
if( a( here-1, here-2 ).ne.zero ) nbnext = 2
end if
call dtgex2( wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, here-nbnext, nbnext, 1, work, lwork, info )
if( info.ne.0 ) then
ilst = here
return
end if
if( nbnext.eq.1 ) then
!
!              swap two 1-by-1 blocks.
!
call dtgex2( wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, here, nbnext, 1, work, lwork, info )
if( info.ne.0 ) then
ilst = here
return
end if
here = here - 1
else
!
!             recompute nbnext in case of 2-by-2 split.
!
if( a( here, here-1 ).eq.zero ) nbnext = 1
if( nbnext.eq.2 ) then
!
!                 2-by-2 block did not split.
!
call dtgex2( wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, here-1, 2, 1, work, lwork, info )
if( info.ne.0 ) then
ilst = here
return
end if
here = here - 2
else
!
!                 2-by-2 block did split.
!
call dtgex2( wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, here, 1, 1, work, lwork, info )
if( info.ne.0 ) then
ilst = here
return
end if
here = here - 1
call dtgex2( wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, here, 1, 1, work, lwork, info )
if( info.ne.0 ) then
ilst = here
return
end if
here = here - 1
end if
end if
end if
if( here.gt.ilst ) go to 20
end if
ilst = here
work( 1 ) = lwmin
return
!
!     end of dtgexc
!
end subroutine dtgexc
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/INSTALL/dlamch.f
!> \brief \b dlamch
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  definition:
!  ===========
!
!      double precision function dlamch( cmach )
!
!     .. scalar arguments ..
!     character          cmach
!     ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dlamch determines double precision machine parameters.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] cmach
!> \verbatim
!>          cmach is character*1
!>          specifies the value to be returned by dlamch:
!>          = 'e' or 'e',   dlamch := eps
!>          = 's' or 's ,   dlamch := sfmin
!>          = 'b' or 'b',   dlamch := base
!>          = 'p' or 'p',   dlamch := eps*base
!>          = 'n' or 'n',   dlamch := t
!>          = 'r' or 'r',   dlamch := rnd
!>          = 'm' or 'm',   dlamch := emin
!>          = 'u' or 'u',   dlamch := rmin
!>          = 'l' or 'l',   dlamch := emax
!>          = 'o' or 'o',   dlamch := rmax
!>          where
!>          eps   = relative machine precision
!>          sfmin = safe minimum, such that 1/sfmin does not overflow
!>          base  = base of the machine
!>          prec  = eps*base
!>          t     = number of (base) digits in the mantissa
!>          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
!>          emin  = minimum exponent before (gradual) underflow
!>          rmin  = underflow threshold - base**(emin-1)
!>          emax  = largest exponent before overflow
!>          rmax  = overflow threshold  - (base**emax)*(1-eps)
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup lamch
!
!  =====================================================================
real(dp) function dlamch( cmach ) result(lapack_result)
implicit none
!
!  -- lapack auxiliary routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
character, intent(in) :: cmach !! specifies the value to be returned by dlamch:
!     ..
!
! =====================================================================
!
!     .. parameters ..
real(dp) :: one, zero
parameter          ( one = 1.0d+0, zero = 0.0d+0 )
!     ..
!     .. local scalars ..
real(dp) :: rnd, eps, sfmin, small, rmach
!     ..
!     .. external functions ..
!     ..
!     .. intrinsic functions ..
intrinsic          digits, epsilon, huge, maxexponent, minexponent, radix, tiny
!     ..
!     .. executable statements ..
!
!
!     assume rounding, not chopping. always.
!
rnd = one
!
if( one.eq.rnd ) then
eps = epsilon(zero) * 0.5
else
eps = epsilon(zero)
end if
!
if( lsame( cmach, 'e' ) ) then
rmach = eps
else if( lsame( cmach, 's' ) ) then
sfmin = tiny(zero)
small = one / huge(zero)
if( small.ge.sfmin ) then
!
!           use small plus a bit, to avoid the possibility of rounding
!           causing overflow when computing  1/sfmin.
!
sfmin = small*( one+eps )
end if
rmach = sfmin
else if( lsame( cmach, 'b' ) ) then
rmach = radix(zero)
else if( lsame( cmach, 'p' ) ) then
rmach = eps * radix(zero)
else if( lsame( cmach, 'n' ) ) then
rmach = digits(zero)
else if( lsame( cmach, 'r' ) ) then
rmach = rnd
else if( lsame( cmach, 'm' ) ) then
rmach = minexponent(zero)
else if( lsame( cmach, 'u' ) ) then
rmach = tiny(zero)
else if( lsame( cmach, 'l' ) ) then
rmach = maxexponent(zero)
else if( lsame( cmach, 'o' ) ) then
rmach = huge(zero)
else
rmach = zero
end if
!
lapack_result = rmach
return
!
!     end of dlamch
!
end function dlamch
!***********************************************************************
!> \brief \b dlamc3
!> \details
!> \b purpose:
!> \verbatim
!> dlamc3  is intended to force  a  and  b  to be stored prior to doing
!> the addition of  a  and  b ,  for use in situations where optimizers
!> might hold one of these in a register.
!> \endverbatim
!> \author lapack is a software package provided by univ. of tennessee, univ. of california
!berkeley, univ. of colorado denver and nag ltd..
!> \param[in] a
!> \verbatim
!>          a is a double precision
!> \endverbatim
!>
!> \param[in] b
!> \verbatim
!>          b is a double precision
!>          the values a and b.
!> \endverbatim
!>
!> \ingroup lamc3
!>
real(dp) function dlamc3( a, b ) result(lapack_result)
implicit none
!
!  -- lapack auxiliary routine --
!     univ. of tennessee, univ. of california berkeley and nag ltd..
!
!     .. scalar arguments ..
real(dp), intent(in) :: a !! first scalar addend, retained in the original evaluation order.
real(dp), intent(in) :: b !! the values a and b.
!     ..
! =====================================================================
!
!     .. executable statements ..
!
lapack_result = a + b
!
return
!
!     end of dlamc3
!
end function dlamc3
!
!***********************************************************************
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlanhs.f
!> \brief \b dlanhs returns the value of the 1-norm, frobenius norm, infinity-norm, or the
!largest absolute value of any element of an upper hessenberg matrix.
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download dlanhs + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlanhs.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlanhs.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlanhs.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       double precision function dlanhs( norm, n, a, lda, work )
!
!       .. scalar arguments ..
!       character          norm
!       integer            lda, n
!       ..
!       .. array arguments ..
!       double precision   a( lda, * ), work( * )
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dlanhs  returns the value of the one norm,  or the frobenius norm, or
!> the  infinity norm,  or the  element of  largest absolute value  of a
!> hessenberg matrix a.
!> \endverbatim
!>
!> \return dlanhs
!> \verbatim
!>
!>    dlanhs = ( max(abs(a(i,j))), norm = 'm' or 'm'
!>             (
!>             ( norm1(a),         norm = '1', 'o' or 'o'
!>             (
!>             ( normi(a),         norm = 'i' or 'i'
!>             (
!>             ( normf(a),         norm = 'f', 'f', 'e' or 'e'
!>
!> where  norm1  denotes the  one norm of a matrix (maximum column sum),
!> normi  denotes the  infinity norm  of a matrix  (maximum row sum) and
!> normf  denotes the  frobenius norm of a matrix (square root of sum of
!> squares).  note that  max(abs(a(i,j)))  is not a consistent matrix norm.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] norm
!> \verbatim
!>          norm is character*1
!>          specifies the value to be returned in dlanhs as described
!>          above.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is integer
!>          the order of the matrix a.  n >= 0.  when n = 0, dlanhs is
!>          set to zero.
!> \endverbatim
!>
!> \param[in] a
!> \verbatim
!>          a is double precision array, dimension (lda,n)
!>          the n by n upper hessenberg matrix a; the part of a below the
!>          first sub-diagonal is not referenced.
!> \endverbatim
!>
!> \param[in] lda
!> \verbatim
!>          lda is integer
!>          the leading dimension of the array a.  lda >= max(n,1).
!> \endverbatim
!>
!> \param[out] work
!> \verbatim
!>          work is double precision array, dimension (max(1,lwork)),
!>          where lwork >= n when norm = 'i'; otherwise, work is not
!>          referenced.
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup lanhs
!
!  =====================================================================
real(dp) function dlanhs( norm, n, a, lda, work ) result(lapack_result)
implicit none
!
!  -- lapack auxiliary routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
character, intent(in) :: norm !! specifies the value to be returned in dlanhs as described above.
integer, intent(in) :: lda !! the leading dimension of the array a.  lda >= max(n,1).
integer, intent(in) :: n !! the order of the matrix a.  n >= 0.  when n = 0, dlanhs is set to zero.
!     ..
!     .. array arguments ..
real(dp), intent(in) :: a( lda, * ) !! the n by n upper hessenberg matrix a; the part of a below
!! the first sub-diagonal is not referenced.
real(dp), intent(out) :: work( * ) !! where lwork >= n when norm = 'i'; otherwise, work is not
!! referenced.
!     ..
!
! =====================================================================
!
!     .. parameters ..
real(dp) :: one, zero
parameter          ( one = 1.0d+0, zero = 0.0d+0 )
!     ..
!     .. local scalars ..
integer :: i, j
real(dp) :: scale, sum, value
!     ..
!     .. external subroutines ..
!     ..
!     .. external functions ..
!     ..
!     .. intrinsic functions ..
intrinsic          abs, min, sqrt
!     ..
!     .. executable statements ..
!
if( n.eq.0 ) then
value = zero
else if( lsame( norm, 'm' ) ) then
!
!        find max(abs(a(i,j))).
!
value = zero
do 20 j = 1, n
do 10 i = 1, min( n, j+1 )
sum = abs( a( i, j ) )
if( value .lt. sum .or. disnan( sum ) ) value = sum
10 continue
20 continue
else if( ( lsame( norm, 'o' ) ) .or. ( norm.eq.'1' ) ) then
!
!        find norm1(a).
!
value = zero
do 40 j = 1, n
sum = zero
do 30 i = 1, min( n, j+1 )
sum = sum + abs( a( i, j ) )
30 continue
if( value .lt. sum .or. disnan( sum ) ) value = sum
40 continue
else if( lsame( norm, 'i' ) ) then
!
!        find normi(a).
!
do 50 i = 1, n
work( i ) = zero
50 continue
do 70 j = 1, n
do 60 i = 1, min( n, j+1 )
work( i ) = work( i ) + abs( a( i, j ) )
60 continue
70 continue
value = zero
do 80 i = 1, n
sum = work( i )
if( value .lt. sum .or. disnan( sum ) ) value = sum
80 continue
else if( ( lsame( norm, 'f' ) ) .or. ( lsame( norm, 'e' ) ) ) then
!
!        find normf(a).
!
scale = zero
sum = one
do 90 j = 1, n
call dlassq( min( n, j+1 ), a( 1, j ), 1, scale, sum )
90 continue
value = scale*sqrt( sum )
end if
!
lapack_result = value
return
!
!     end of dlanhs
!
end function dlanhs
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlapy2.f
!> \brief \b dlapy2 returns sqrt(x2+y2).
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download dlapy2 + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlapy2.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlapy2.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlapy2.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       double precision function dlapy2( x, y )
!
!       .. scalar arguments ..
!       double precision   x, y
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dlapy2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
!> overflow and unnecessary underflow.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] x
!> \verbatim
!>          x is double precision
!> \endverbatim
!>
!> \param[in] y
!> \verbatim
!>          y is double precision
!>          x and y specify the values x and y.
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup lapy2
!
!  =====================================================================
real(dp) function dlapy2( x, y ) result(lapack_result)
implicit none
!
!  -- lapack auxiliary routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
real(dp), intent(in) :: x !! input component or strided vector used by this norm or nan check.
real(dp), intent(in) :: y !! x and y specify the values x and y.
!     ..
!
!  =====================================================================
!
!     .. parameters ..
real(dp) :: zero
parameter          ( zero = 0.0d0 )
real(dp) :: one
parameter          ( one = 1.0d0 )
!     ..
!     .. local scalars ..
real(dp) :: w, xabs, yabs, z, hugeval
logical :: x_is_nan, y_is_nan
!     ..
!     .. external functions ..
!     ..
!     .. external subroutines ..
!     ..
!     .. intrinsic functions ..
intrinsic          abs, max, min, sqrt
!     ..
!     .. executable statements ..
!
x_is_nan = disnan( x )
y_is_nan = disnan( y )
if ( x_is_nan ) lapack_result = x
if ( y_is_nan ) lapack_result = y
hugeval = dlamch( 'overflow' )
!
if ( .not.( x_is_nan.or.y_is_nan ) ) then
xabs = abs( x )
yabs = abs( y )
w = max( xabs, yabs )
z = min( xabs, yabs )
if( z.eq.zero .or. w.gt.hugeval ) then
lapack_result = w
else
lapack_result = w*sqrt( one+( z / w )**2 )
end if
end if
return
!
!     end of dlapy2
!
end function dlapy2
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlapy3.f
!> \brief \b dlapy3 returns sqrt(x2+y2+z2).
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download dlapy3 + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlapy3.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlapy3.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlapy3.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       double precision function dlapy3( x, y, z )
!
!       .. scalar arguments ..
!       double precision   x, y, z
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dlapy3 returns sqrt(x**2+y**2+z**2), taking care not to cause
!> unnecessary overflow and unnecessary underflow.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] x
!> \verbatim
!>          x is double precision
!> \endverbatim
!>
!> \param[in] y
!> \verbatim
!>          y is double precision
!> \endverbatim
!>
!> \param[in] z
!> \verbatim
!>          z is double precision
!>          x, y and z specify the values x, y and z.
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup lapy3
!
!  =====================================================================
real(dp) function dlapy3( x, y, z ) result(lapack_result)
implicit none
!
!  -- lapack auxiliary routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
real(dp), intent(in) :: x !! input component or strided vector used by this norm or nan check.
real(dp), intent(in) :: y !! second input component used by this norm or nan check.
real(dp), intent(in) :: z !! x, y and z specify the values x, y and z.
!     ..
!
!  =====================================================================
!
!     .. parameters ..
real(dp) :: zero
parameter          ( zero = 0.0d0 )
!     ..
!     .. local scalars ..
real(dp) :: w, xabs, yabs, zabs, hugeval
!     ..
!     .. external subroutines ..
!     ..
!     .. intrinsic functions ..
intrinsic          abs, max, sqrt
!     ..
!     .. executable statements ..
!
hugeval = dlamch( 'overflow' )
xabs = abs( x )
yabs = abs( y )
zabs = abs( z )
w = max( xabs, yabs, zabs )
if( w.eq.zero .or. w.gt.hugeval ) then
!     w can be zero for max(0,nan,0)
!     adding all three entries together will make sure
!     nan will not disappear.
lapack_result =  xabs + yabs + zabs
else
lapack_result = w*sqrt( ( xabs / w )**2+( yabs / w )**2+ ( zabs / w )**2 )
end if
return
!
!     end of dlapy3
!
end function dlapy3
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlag2.f
!> \brief \b dlag2 computes the eigenvalues of a 2-by-2 generalized eigenvalue problem, with
!scaling as necessary to avoid over-/underflow.
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download dlag2 + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlag2.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlag2.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlag2.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       subroutine dlag2( a, lda, b, ldb, safmin, scale1, scale2, wr1,
!                         wr2, wi )
!
!       .. scalar arguments ..
!       integer            lda, ldb
!       double precision   safmin, scale1, scale2, wi, wr1, wr2
!       ..
!       .. array arguments ..
!       double precision   a( lda, * ), b( ldb, * )
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dlag2 computes the eigenvalues of a 2 x 2 generalized eigenvalue
!> problem  a - w b, with scaling as necessary to avoid over-/underflow.
!>
!> the scaling factor "s" results in a modified eigenvalue equation
!>
!>     s a - w b
!>
!> where  s  is a non-negative scaling factor chosen so that  w,  w b,
!> and  s a  do not overflow and, if possible, do not underflow, either.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] a
!> \verbatim
!>          a is double precision array, dimension (lda, 2)
!>          on entry, the 2 x 2 matrix a.  it is assumed that its 1-norm
!>          is less than 1/safmin.  entries less than
!>          sqrt(safmin)*norm(a) are subject to being treated as zero.
!> \endverbatim
!>
!> \param[in] lda
!> \verbatim
!>          lda is integer
!>          the leading dimension of the array a.  lda >= 2.
!> \endverbatim
!>
!> \param[in] b
!> \verbatim
!>          b is double precision array, dimension (ldb, 2)
!>          on entry, the 2 x 2 upper triangular matrix b.  it is
!>          assumed that the one-norm of b is less than 1/safmin.  the
!>          diagonals should be at least sqrt(safmin) times the largest
!>          element of b (in absolute value); if a diagonal is smaller
!>          than that, then  +/- sqrt(safmin) will be used instead of
!>          that diagonal.
!> \endverbatim
!>
!> \param[in] ldb
!> \verbatim
!>          ldb is integer
!>          the leading dimension of the array b.  ldb >= 2.
!> \endverbatim
!>
!> \param[in] safmin
!> \verbatim
!>          safmin is double precision
!>          the smallest positive number s.t. 1/safmin does not
!>          overflow.  (this should always be dlamch('s') -- it is an
!>          argument in order to avoid having to call dlamch frequently.)
!> \endverbatim
!>
!> \param[out] scale1
!> \verbatim
!>          scale1 is double precision
!>          a scaling factor used to avoid over-/underflow in the
!>          eigenvalue equation which defines the first eigenvalue.  if
!>          the eigenvalues are complex, then the eigenvalues are
!>          ( wr1  +/-  wi i ) / scale1  (which may lie outside the
!>          exponent range of the machine), scale1=scale2, and scale1
!>          will always be positive.  if the eigenvalues are real, then
!>          the first (real) eigenvalue is  wr1 / scale1 , but this may
!>          overflow or underflow, and in fact, scale1 may be zero or
!>          less than the underflow threshold if the exact eigenvalue
!>          is sufficiently large.
!> \endverbatim
!>
!> \param[out] scale2
!> \verbatim
!>          scale2 is double precision
!>          a scaling factor used to avoid over-/underflow in the
!>          eigenvalue equation which defines the second eigenvalue.  if
!>          the eigenvalues are complex, then scale2=scale1.  if the
!>          eigenvalues are real, then the second (real) eigenvalue is
!>          wr2 / scale2 , but this may overflow or underflow, and in
!>          fact, scale2 may be zero or less than the underflow
!>          threshold if the exact eigenvalue is sufficiently large.
!> \endverbatim
!>
!> \param[out] wr1
!> \verbatim
!>          wr1 is double precision
!>          if the eigenvalue is real, then wr1 is scale1 times the
!>          eigenvalue closest to the (2,2) element of a b**(-1).  if the
!>          eigenvalue is complex, then wr1=wr2 is scale1 times the real
!>          part of the eigenvalues.
!> \endverbatim
!>
!> \param[out] wr2
!> \verbatim
!>          wr2 is double precision
!>          if the eigenvalue is real, then wr2 is scale2 times the
!>          other eigenvalue.  if the eigenvalue is complex, then
!>          wr1=wr2 is scale1 times the real part of the eigenvalues.
!> \endverbatim
!>
!> \param[out] wi
!> \verbatim
!>          wi is double precision
!>          if the eigenvalue is real, then wi is zero.  if the
!>          eigenvalue is complex, then wi is scale1 times the imaginary
!>          part of the eigenvalues.  wi will always be non-negative.
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup lag2
!
!  =====================================================================
subroutine dlag2( a, lda, b, ldb, safmin, scale1, scale2, wr1, wr2, wi )
implicit none
!
!  -- lapack auxiliary routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
integer, intent(in) :: lda !! the leading dimension of the array a.  lda >= 2.
integer, intent(in) :: ldb !! the leading dimension of the array b.  ldb >= 2.
real(dp), intent(in) :: safmin !! the smallest positive number s.t. 1/safmin does not overflow.
!! (this should always be dlamch('s') -- it is an
real(dp), intent(out) :: scale1 !! a scaling factor used to avoid over-/underflow in the
!! eigenvalue equation which defines the first eigenvalue.  if
real(dp), intent(out) :: scale2 !! a scaling factor used to avoid over-/underflow in the
!! eigenvalue equation which defines the second eigenvalue.  if
real(dp), intent(out) :: wi !! if the eigenvalue is real, then wi is zero.  if the eigenvalue is
!! complex, then wi is scale1 times the imaginary
real(dp), intent(out) :: wr1 !! if the eigenvalue is real, then wr1 is scale1 times the eigenvalue
!! closest to the (2,2) element of a b**(-1).  if the
real(dp), intent(out) :: wr2 !! if the eigenvalue is real, then wr2 is scale2 times the other
!! eigenvalue.  if the eigenvalue is complex, then
!     ..
!     .. array arguments ..
real(dp), intent(in) :: a( lda, * ) !! on entry, the 2 x 2 matrix a.  it is assumed that its
!! 1-norm is less than 1/safmin.  entries less than
real(dp), intent(in) :: b( ldb, * ) !! on entry, the 2 x 2 upper triangular matrix b.  it is
!! assumed that the one-norm of b is less than 1/safmin.  the
!     ..
!
!  =====================================================================
!
!     .. parameters ..
real(dp) :: zero, one, two
parameter          ( zero = 0.0d+0, one = 1.0d+0, two = 2.0d+0 )
real(dp) :: half
parameter          ( half = one / two )
real(dp) :: fuzzy1
parameter          ( fuzzy1 = one+1.0d-5 )
!     ..
!     .. local scalars ..
real(dp) :: a11, a12, a21, a22, abi22, anorm, as11, as12, as22, ascale, b11, b12, b22, binv11, binv22, &
    bmin, bnorm, bscale, bsize, c1, c2, c3, c4, c5, diff, discr, pp, qq, r, rtmax, rtmin, s1, s2, safmax, &
    shift, ss, sum, wabs, wbig, wdet, wscale, wsize, wsmall
!     ..
!     .. intrinsic functions ..
intrinsic          abs, max, min, sign, sqrt
!     ..
!     .. executable statements ..
!
rtmin = sqrt( safmin )
rtmax = one / rtmin
safmax = one / safmin
!
!     scale a
!
anorm = max( abs( a( 1, 1 ) )+abs( a( 2, 1 ) ), abs( a( 1, 2 ) )+abs( a( 2, 2 ) ), safmin )
ascale = one / anorm
a11 = ascale*a( 1, 1 )
a21 = ascale*a( 2, 1 )
a12 = ascale*a( 1, 2 )
a22 = ascale*a( 2, 2 )
!
!     perturb b if necessary to insure non-singularity
!
b11 = b( 1, 1 )
b12 = b( 1, 2 )
b22 = b( 2, 2 )
bmin = rtmin*max( abs( b11 ), abs( b12 ), abs( b22 ), rtmin )
if( abs( b11 ).lt.bmin ) b11 = sign( bmin, b11 )
if( abs( b22 ).lt.bmin ) b22 = sign( bmin, b22 )
!
!     scale b
!
bnorm = max( abs( b11 ), abs( b12 )+abs( b22 ), safmin )
bsize = max( abs( b11 ), abs( b22 ) )
bscale = one / bsize
b11 = b11*bscale
b12 = b12*bscale
b22 = b22*bscale
!
!     compute larger eigenvalue by method described by c. van loan
!
!     ( as is a shifted by -shift*b )
!
binv11 = one / b11
binv22 = one / b22
s1 = a11*binv11
s2 = a22*binv22
if( abs( s1 ).le.abs( s2 ) ) then
as12 = a12 - s1*b12
as22 = a22 - s1*b22
ss = a21*( binv11*binv22 )
abi22 = as22*binv22 - ss*b12
pp = half*abi22
shift = s1
else
as12 = a12 - s2*b12
as11 = a11 - s2*b11
ss = a21*( binv11*binv22 )
abi22 = -ss*b12
pp = half*( as11*binv11+abi22 )
shift = s2
end if
qq = ss*as12
if( abs( pp*rtmin ).ge.one ) then
discr = ( rtmin*pp )**2 + qq*safmin
r = sqrt( abs( discr ) )*rtmax
else
if( pp**2+abs( qq ).le.safmin ) then
discr = ( rtmax*pp )**2 + qq*safmax
r = sqrt( abs( discr ) )*rtmin
else
discr = pp**2 + qq
r = sqrt( abs( discr ) )
end if
end if
!
!     note: the test of r in the following if is to cover the case when
!           discr is small and negative and is flushed to zero during
!           the calculation of r.  on machines which have a consistent
!           flush-to-zero threshold and handle numbers above that
!           threshold correctly, it would not be necessary.
!
if( discr.ge.zero .or. r.eq.zero ) then
sum = pp + sign( r, pp )
diff = pp - sign( r, pp )
wbig = shift + sum
!
!        compute smaller eigenvalue
!
wsmall = shift + diff
if( half*abs( wbig ).gt.max( abs( wsmall ), safmin ) ) then
wdet = ( a11*a22-a12*a21 )*( binv11*binv22 )
wsmall = wdet / wbig
end if
!
!        choose (real) eigenvalue closest to 2,2 element of a*b**(-1)
!        for wr1.
!
if( pp.gt.abi22 ) then
wr1 = min( wbig, wsmall )
wr2 = max( wbig, wsmall )
else
wr1 = max( wbig, wsmall )
wr2 = min( wbig, wsmall )
end if
wi = zero
else
!
!        complex eigenvalues
!
wr1 = shift + pp
wr2 = wr1
wi = r
end if
!
!     further scaling to avoid underflow and overflow in computing
!     scale1 and overflow in computing w*b.
!
!     this scale factor (wscale) is bounded from above using c1 and c2,
!     and from below using c3 and c4.
!        c1 implements the condition  s a  must never overflow.
!        c2 implements the condition  w b  must never overflow.
!        c3, with c2,
!           implement the condition that s a - w b must never overflow.
!        c4 implements the condition  s    should not underflow.
!        c5 implements the condition  max(s,|w|) should be at least 2.
!
c1 = bsize*( safmin*max( one, ascale ) )
c2 = safmin*max( one, bnorm )
c3 = bsize*safmin
if( ascale.le.one .and. bsize.le.one ) then
c4 = min( one, ( ascale / safmin )*bsize )
else
c4 = one
end if
if( ascale.le.one .or. bsize.le.one ) then
c5 = min( one, ascale*bsize )
else
c5 = one
end if
!
!     scale first eigenvalue
!
wabs = abs( wr1 ) + abs( wi )
wsize = max( safmin, c1, fuzzy1*( wabs*c2+c3 ), min( c4, half*max( wabs, c5 ) ) )
if( wsize.ne.one ) then
wscale = one / wsize
if( wsize.gt.one ) then
scale1 = ( max( ascale, bsize )*wscale )* min( ascale, bsize )
else
scale1 = ( min( ascale, bsize )*wscale )* max( ascale, bsize )
end if
wr1 = wr1*wscale
if( wi.ne.zero ) then
wi = wi*wscale
wr2 = wr1
scale2 = scale1
end if
else
scale1 = ascale*bsize
scale2 = scale1
end if
!
!     scale second eigenvalue (if real)
!
if( wi.eq.zero ) then
wsize = max( safmin, c1, fuzzy1*( abs( wr2 )*c2+c3 ), min( c4, half*max( abs( wr2 ), c5 ) ) )
if( wsize.ne.one ) then
wscale = one / wsize
if( wsize.gt.one ) then
scale2 = ( max( ascale, bsize )*wscale )* min( ascale, bsize )
else
scale2 = ( min( ascale, bsize )*wscale )* max( ascale, bsize )
end if
wr2 = wr2*wscale
else
scale2 = ascale*bsize
end if
end if
!
!     end of dlag2
!
return
end subroutine dlag2
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlarfg.f
!> \brief \b dlarfg generates an elementary reflector (householder matrix).
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download dlarfg + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarfg.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarfg.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarfg.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       subroutine dlarfg( n, alpha, x, incx, tau )
!
!       .. scalar arguments ..
!       integer            incx, n
!       double precision   alpha, tau
!       ..
!       .. array arguments ..
!       double precision   x( * )
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dlarfg generates a real elementary reflector h of order n, such
!> that
!>
!>       h * ( alpha ) = ( beta ),   h**t * h = i.
!>           (   x   )   (   0  )
!>
!> where alpha and beta are scalars, and x is an (n-1)-element real
!> vector. h is represented in the form
!>
!>       h = i - tau * ( 1 ) * ( 1 v**t ) ,
!>                     ( v )
!>
!> where tau is a real scalar and v is a real (n-1)-element
!> vector.
!>
!> if the elements of x are all zero, then tau = 0 and h is taken to be
!> the unit matrix.
!>
!> otherwise  1 <= tau <= 2.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] n
!> \verbatim
!>          n is integer
!>          the order of the elementary reflector.
!> \endverbatim
!>
!> \param[in,out] alpha
!> \verbatim
!>          alpha is double precision
!>          on entry, the value alpha.
!>          on exit, it is overwritten with the value beta.
!> \endverbatim
!>
!> \param[in,out] x
!> \verbatim
!>          x is double precision array, dimension
!>                         (1+(n-2)*abs(incx))
!>          on entry, the vector x.
!>          on exit, it is overwritten with the vector v.
!> \endverbatim
!>
!> \param[in] incx
!> \verbatim
!>          incx is integer
!>          the increment between elements of x. incx > 0.
!> \endverbatim
!>
!> \param[out] tau
!> \verbatim
!>          tau is double precision
!>          the value tau.
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup larfg
!
!  =====================================================================
subroutine dlarfg( n, alpha, x, incx, tau )
implicit none
!
!  -- lapack auxiliary routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
integer, intent(in) :: incx !! the increment between elements of x. incx > 0.
integer, intent(in) :: n !! the order of the elementary reflector.
real(dp), intent(inout) :: alpha !! on entry, the value alpha. on exit, it is overwritten with the
!! value beta.
real(dp), intent(out) :: tau !! the value tau.
!     ..
!     .. array arguments ..
real(dp), intent(inout) :: x( * ) !! (1+(n-2)*abs(incx)) on entry, the vector x. on exit, it is
!! overwritten with the vector v.
!     ..
!
!  =====================================================================
!
!     .. parameters ..
real(dp) :: one, zero
parameter          ( one = 1.0d+0, zero = 0.0d+0 )
!     ..
!     .. local scalars ..
integer :: j, knt
real(dp) :: beta, rsafmn, safmin, xnorm
!     ..
!     .. external functions ..
!     ..
!     .. intrinsic functions ..
intrinsic          abs, sign
!     ..
!     .. external subroutines ..
!     ..
!     .. executable statements ..
!
if( n.le.1 ) then
tau = zero
return
end if
!
xnorm = dnrm2( n-1, x, incx )
!
if( xnorm.eq.zero ) then
!
!        h  =  i
!
tau = zero
else
!
!        general case
!
beta = -sign( dlapy2( alpha, xnorm ), alpha )
safmin = dlamch( 's' ) / dlamch( 'e' )
knt = 0
if( abs( beta ).lt.safmin ) then
!
!           xnorm, beta may be inaccurate; scale x and recompute them
!
rsafmn = one / safmin
10 continue
knt = knt + 1
call dscal( n-1, rsafmn, x, incx )
beta = beta*rsafmn
alpha = alpha*rsafmn
if( (abs( beta ).lt.safmin) .and. (knt .lt. 20) ) go to 10
!
!           new beta is at most 1, at least safmin
!
xnorm = dnrm2( n-1, x, incx )
beta = -sign( dlapy2( alpha, xnorm ), alpha )
end if
tau = ( beta-alpha ) / beta
call dscal( n-1, one / ( alpha-beta ), x, incx )
!
!        if alpha is subnormal, it may lose relative accuracy
!
do 20 j = 1, knt
beta = beta*safmin
20 continue
alpha = beta
end if
!
return
!
!     end of dlarfg
!
end subroutine dlarfg
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlartg.f90
!> \brief \b dlartg generates a plane rotation with real cosine and real sine.
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  definition:
!  ===========
!
!       subroutine dlartg( f, g, c, s, r )
!
!       .. scalar arguments ..
!       real(wp)          c, f, g, r, s
!       ..
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dlartg generates a plane rotation so that
!>
!>    [  c  s  ]  .  [ f ]  =  [ r ]
!>    [ -s  c  ]     [ g ]     [ 0 ]
!>
!> where c**2 + s**2 = 1.
!>
!> the mathematical formulas used for c and s are
!>    r = sign(f) * sqrt(f**2 + g**2)
!>    c = f / r
!>    s = g / r
!> hence c >= 0. the algorithm used to compute these quantities
!> incorporates scaling to avoid overflow or underflow in computing the
!> square root of the sum of squares.
!>
!> this version is discontinuous in r at f = 0 but it returns the same
!> c and s as zlartg for complex inputs (f,0) and (g,0).
!>
!> this is a more accurate version of the blas1 routine drotg,
!> with the following other differences:
!>    f and g are unchanged on return.
!>    if g=0, then c=1 and s=0.
!>    if f=0 and (g .ne. 0), then c=0 and s=sign(1,g) without doing any
!>       floating point operations (saves work in dbdsqr when
!>       there are zeros on the diagonal).
!>
!> below, wp=>dp stands for double precision from la_constants module.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] f
!> \verbatim
!>          f is real(wp)
!>          the first component of vector to be rotated.
!> \endverbatim
!>
!> \param[in] g
!> \verbatim
!>          g is real(wp)
!>          the second component of vector to be rotated.
!> \endverbatim
!>
!> \param[out] c
!> \verbatim
!>          c is real(wp)
!>          the cosine of the rotation.
!> \endverbatim
!>
!> \param[out] s
!> \verbatim
!>          s is real(wp)
!>          the sine of the rotation.
!> \endverbatim
!>
!> \param[out] r
!> \verbatim
!>          r is real(wp)
!>          the nonzero component of the rotated vector.
!> \endverbatim
!
!  authors:
!  ========
!
!> \author edward anderson, lockheed martin
!
!> \date july 2016
!
!> \ingroup lartg
!
!> \par contributors:
!  ==================
!>
!> weslley pereira, university of colorado denver, usa
!
!> \par further details:
!  =====================
!>
!> \verbatim
!>
!>  anderson e. (2017)
!>  algorithm 978: safe scaling in the level 1 blas
!>  acm trans math softw 44:1--28
!>  https://doi.org/10.1145/3061665
!>
!> \endverbatim
!
subroutine dlartg( f, g, c, s, r )
implicit none
!
!  -- lapack auxiliary routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!     february 2021
!
!  .. scalar arguments ..
real(dp), intent(out) :: c !! the cosine of the rotation.
real(dp), intent(in) :: f !! the first component of vector to be rotated.
real(dp), intent(in) :: g !! the second component of vector to be rotated.
real(dp), intent(out) :: r !! the nonzero component of the rotated vector.
real(dp), intent(out) :: s !! the sine of the rotation.
!  ..
!  .. local scalars ..
real(dp) :: d, f1, fs, g1, gs, u, rtmin, rtmax
!  ..
!  .. intrinsic functions ..
intrinsic :: abs, sign, sqrt
!  ..
!  .. constants ..
rtmin = sqrt( safmin )
rtmax = sqrt( safmax/2 )
!  ..
!  .. executable statements ..
!
f1 = abs( f )
g1 = abs( g )
if( g == zero ) then
c = one
s = zero
r = f
else if( f == zero ) then
c = zero
s = sign( one, g )
r = g1
else if( f1 > rtmin .and. f1 < rtmax .and. g1 > rtmin .and. g1 < rtmax ) then
d = sqrt( f*f + g*g )
c = f1 / d
r = sign( d, f )
s = g / r
else
u = min( safmax, max( safmin, f1, g1 ) )
fs = f / u
gs = g / u
d = sqrt( fs*fs + gs*gs )
c = abs( fs ) / d
r = sign( d, f )
s = gs / r
r = r*u
end if
return
end subroutine dlartg
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlaset.f
!> \brief \b dlaset initializes the off-diagonal elements and the diagonal elements of a matrix
!to given values.
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download dlaset + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaset.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaset.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaset.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       subroutine dlaset( uplo, m, n, alpha, beta, a, lda )
!
!       .. scalar arguments ..
!       character          uplo
!       integer            lda, m, n
!       double precision   alpha, beta
!       ..
!       .. array arguments ..
!       double precision   a( lda, * )
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dlaset initializes an m-by-n matrix a to beta on the diagonal and
!> alpha on the offdiagonals.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] uplo
!> \verbatim
!>          uplo is character*1
!>          specifies the part of the matrix a to be set.
!>          = 'u':      upper triangular part is set; the strictly lower
!>                      triangular part of a is not changed.
!>          = 'l':      lower triangular part is set; the strictly upper
!>                      triangular part of a is not changed.
!>          otherwise:  all of the matrix a is set.
!> \endverbatim
!>
!> \param[in] m
!> \verbatim
!>          m is integer
!>          the number of rows of the matrix a.  m >= 0.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is integer
!>          the number of columns of the matrix a.  n >= 0.
!> \endverbatim
!>
!> \param[in] alpha
!> \verbatim
!>          alpha is double precision
!>          the constant to which the offdiagonal elements are to be set.
!> \endverbatim
!>
!> \param[in] beta
!> \verbatim
!>          beta is double precision
!>          the constant to which the diagonal elements are to be set.
!> \endverbatim
!>
!> \param[out] a
!> \verbatim
!>          a is double precision array, dimension (lda,n)
!>          on exit, the leading m-by-n submatrix of a is set as follows:
!>
!>          if uplo = 'u', a(i,j) = alpha, 1<=i<=j-1, 1<=j<=n,
!>          if uplo = 'l', a(i,j) = alpha, j+1<=i<=m, 1<=j<=n,
!>          otherwise,     a(i,j) = alpha, 1<=i<=m, 1<=j<=n, i.ne.j,
!>
!>          and, for all uplo, a(i,i) = beta, 1<=i<=min(m,n).
!> \endverbatim
!>
!> \param[in] lda
!> \verbatim
!>          lda is integer
!>          the leading dimension of the array a.  lda >= max(1,m).
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup laset
!
!  =====================================================================
subroutine dlaset( uplo, m, n, alpha, beta, a, lda )
implicit none
!
!  -- lapack auxiliary routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
character, intent(in) :: uplo !! specifies the part of the matrix a to be set.
integer, intent(in) :: lda !! the leading dimension of the array a.  lda >= max(1,m).
integer, intent(in) :: m !! the number of rows of the matrix a.  m >= 0.
integer, intent(in) :: n !! the number of columns of the matrix a.  n >= 0.
real(dp), intent(in) :: alpha !! the constant to which the offdiagonal elements are to be set.
real(dp), intent(in) :: beta !! the constant to which the diagonal elements are to be set.
!     ..
!     .. array arguments ..
real(dp), intent(out) :: a( lda, * ) !! on exit, the leading m-by-n submatrix of a is set as
!! follows: if uplo = 'u', a(i,j) = alpha, 1<=i<=j-1,
!! 1<=j<=n,
!     ..
!
! =====================================================================
!
!     .. local scalars ..
integer :: i, j
!     ..
!     .. external functions ..
!     ..
!     .. intrinsic functions ..
intrinsic          min
!     ..
!     .. executable statements ..
!
if( lsame( uplo, 'u' ) ) then
!
!        set the strictly upper triangular or trapezoidal part of the
!        array to alpha.
!
do 20 j = 2, n
do 10 i = 1, min( j-1, m )
a( i, j ) = alpha
10 continue
20 continue
!
else if( lsame( uplo, 'l' ) ) then
!
!        set the strictly lower triangular or trapezoidal part of the
!        array to alpha.
!
do 40 j = 1, min( m, n )
do 30 i = j + 1, m
a( i, j ) = alpha
30 continue
40 continue
!
else
!
!        set the leading m-by-n submatrix to alpha.
!
do 60 j = 1, n
do 50 i = 1, m
a( i, j ) = alpha
50 continue
60 continue
end if
!
!     set the first min(m,n) diagonal elements to beta.
!
do 70 i = 1, min( m, n )
a( i, i ) = beta
70 continue
!
return
!
!     end of dlaset
!
end subroutine dlaset
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlasv2.f
!> \brief \b dlasv2 computes the singular value decomposition of a 2-by-2 triangular matrix.
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download dlasv2 + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasv2.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasv2.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasv2.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       subroutine dlasv2( f, g, h, ssmin, ssmax, snr, csr, snl, csl )
!
!       .. scalar arguments ..
!       double precision   csl, csr, f, g, h, snl, snr, ssmax, ssmin
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dlasv2 computes the singular value decomposition of a 2-by-2
!> triangular matrix
!>    [  f   g  ]
!>    [  0   h  ].
!> on return, abs(ssmax) is the larger singular value, abs(ssmin) is the
!> smaller singular value, and (csl,snl) and (csr,snr) are the left and
!> right singular vectors for abs(ssmax), giving the decomposition
!>
!>    [ csl  snl ] [  f   g  ] [ csr -snr ]  =  [ ssmax   0   ]
!>    [-snl  csl ] [  0   h  ] [ snr  csr ]     [  0    ssmin ].
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] f
!> \verbatim
!>          f is double precision
!>          the (1,1) element of the 2-by-2 matrix.
!> \endverbatim
!>
!> \param[in] g
!> \verbatim
!>          g is double precision
!>          the (1,2) element of the 2-by-2 matrix.
!> \endverbatim
!>
!> \param[in] h
!> \verbatim
!>          h is double precision
!>          the (2,2) element of the 2-by-2 matrix.
!> \endverbatim
!>
!> \param[out] ssmin
!> \verbatim
!>          ssmin is double precision
!>          abs(ssmin) is the smaller singular value.
!> \endverbatim
!>
!> \param[out] ssmax
!> \verbatim
!>          ssmax is double precision
!>          abs(ssmax) is the larger singular value.
!> \endverbatim
!>
!> \param[out] snl
!> \verbatim
!>          snl is double precision
!> \endverbatim
!>
!> \param[out] csl
!> \verbatim
!>          csl is double precision
!>          the vector (csl, snl) is a unit left singular vector for the
!>          singular value abs(ssmax).
!> \endverbatim
!>
!> \param[out] snr
!> \verbatim
!>          snr is double precision
!> \endverbatim
!>
!> \param[out] csr
!> \verbatim
!>          csr is double precision
!>          the vector (csr, snr) is a unit right singular vector for the
!>          singular value abs(ssmax).
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup lasv2
!
!> \par further details:
!  =====================
!>
!> \verbatim
!>
!>  any input parameter may be aliased with any output parameter.
!>
!>  barring over/underflow and assuming a guard digit in subtraction, all
!>  output quantities are correct to within a few units in the last
!>  place (ulps).
!>
!>  in ieee arithmetic, the code works correctly if one matrix element is
!>  infinite.
!>
!>  overflow will not occur unless the largest singular value itself
!>  overflows or is within a few ulps of overflow.
!>
!>  underflow is harmless if underflow is gradual. otherwise, results
!>  may correspond to a matrix modified by perturbations of size near
!>  the underflow threshold.
!> \endverbatim
!>
!  =====================================================================
subroutine dlasv2( f, g, h, ssmin, ssmax, snr, csr, snl, csl )
implicit none
!
!  -- lapack auxiliary routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
real(dp), intent(out) :: csl !! the vector (csl, snl) is a unit left singular vector for the
!! singular value abs(ssmax).
real(dp), intent(out) :: csr !! the vector (csr, snr) is a unit right singular vector for the
!! singular value abs(ssmax).
real(dp), intent(in) :: f !! the (1,1) element of the 2-by-2 matrix.
real(dp), intent(in) :: g !! the (1,2) element of the 2-by-2 matrix.
real(dp), intent(in) :: h !! the (2,2) element of the 2-by-2 matrix.
real(dp), intent(out) :: snl !! sine of the left real plane rotation.
real(dp), intent(out) :: snr !! sine of the right real plane rotation.
real(dp), intent(out) :: ssmax !! abs(ssmax) is the larger singular value.
real(dp), intent(out) :: ssmin !! abs(ssmin) is the smaller singular value.
!     ..
!
! =====================================================================
!
!     .. parameters ..
real(dp) :: zero
parameter          ( zero = 0.0d0 )
real(dp) :: half
parameter          ( half = 0.5d0 )
real(dp) :: one
parameter          ( one = 1.0d0 )
real(dp) :: two
parameter          ( two = 2.0d0 )
real(dp) :: four
parameter          ( four = 4.0d0 )
!     ..
!     .. local scalars ..
logical :: gasmal, swap
integer :: pmax
real(dp) :: a, clt, crt, d, fa, ft, ga, gt, ha, ht, l, m, mm, r, s, slt, srt, t, temp, tsign, tt
!     ..
!     .. intrinsic functions ..
intrinsic          abs, sign, sqrt
!     ..
!     .. external functions ..
!     ..
!     .. executable statements ..
!
ft = f
fa = abs( ft )
ht = h
ha = abs( h )
!
!     pmax points to the maximum absolute element of matrix
!       pmax = 1 if f largest in absolute values
!       pmax = 2 if g largest in absolute values
!       pmax = 3 if h largest in absolute values
!
pmax = 1
swap = ( ha.gt.fa )
if( swap ) then
pmax = 3
temp = ft
ft = ht
ht = temp
temp = fa
fa = ha
ha = temp
!
!        now fa .ge. ha
!
end if
gt = g
ga = abs( gt )
if( ga.eq.zero ) then
!
!        diagonal matrix
!
ssmin = ha
ssmax = fa
clt = one
crt = one
slt = zero
srt = zero
else
gasmal = .true.
if( ga.gt.fa ) then
pmax = 2
if( ( fa / ga ).lt.dlamch( 'eps' ) ) then
!
!              case of very large ga
!
gasmal = .false.
ssmax = ga
if( ha.gt.one ) then
ssmin = fa / ( ga / ha )
else
ssmin = ( fa / ga )*ha
end if
clt = one
slt = ht / gt
srt = one
crt = ft / gt
end if
end if
if( gasmal ) then
!
!           normal case
!
d = fa - ha
if( d.eq.fa ) then
!
!              copes with infinite f or h
!
l = one
else
l = d / fa
end if
!
!           note that 0 .le. l .le. 1
!
m = gt / ft
!
!           note that abs(m) .le. 1/macheps
!
t = two - l
!
!           note that t .ge. 1
!
mm = m*m
tt = t*t
s = sqrt( tt+mm )
!
!           note that 1 .le. s .le. 1 + 1/macheps
!
if( l.eq.zero ) then
r = abs( m )
else
r = sqrt( l*l+mm )
end if
!
!           note that 0 .le. r .le. 1 + 1/macheps
!
a = half*( s+r )
!
!           note that 1 .le. a .le. 1 + abs(m)
!
ssmin = ha / a
ssmax = fa*a
if( mm.eq.zero ) then
!
!              note that m is very tiny
!
if( l.eq.zero ) then
t = sign( two, ft )*sign( one, gt )
else
t = gt / sign( d, ft ) + m / t
end if
else
t = ( m / ( s+t )+m / ( r+l ) )*( one+a )
end if
l = sqrt( t*t+four )
crt = two / l
srt = t / l
clt = ( crt+srt*m ) / a
slt = ( ht / ft )*srt / a
end if
end if
if( swap ) then
csl = srt
snl = crt
csr = slt
snr = clt
else
csl = clt
snl = slt
csr = crt
snr = srt
end if
!
!     correct signs of ssmax and ssmin
!
if( pmax.eq.1 ) tsign = sign( one, csr )*sign( one, csl )*sign( one, f )
if( pmax.eq.2 ) tsign = sign( one, snr )*sign( one, csl )*sign( one, g )
if( pmax.eq.3 ) tsign = sign( one, snr )*sign( one, snl )*sign( one, h )
ssmax = sign( ssmax, tsign )
ssmin = sign( ssmin, tsign*sign( one, f )*sign( one, h ) )
return
!
!     end of dlasv2
!
end subroutine dlasv2
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/BLAS/SRC/drot.f
!> \brief \b drot
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  definition:
!  ===========
!
!       subroutine drot(n,dx,incx,dy,incy,c,s)
!
!       .. scalar arguments ..
!       double precision c,s
!       integer incx,incy,n
!       ..
!       .. array arguments ..
!       double precision dx(*),dy(*)
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!>    drot applies a plane rotation.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] n
!> \verbatim
!>          n is integer
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in,out] dx
!> \verbatim
!>          dx is double precision array, dimension ( 1 + ( n - 1 )*abs( incx ) )
!> \endverbatim
!>
!> \param[in] incx
!> \verbatim
!>          incx is integer
!>         storage spacing between elements of dx
!> \endverbatim
!>
!> \param[in,out] dy
!> \verbatim
!>          dy is double precision array, dimension ( 1 + ( n - 1 )*abs( incy ) )
!> \endverbatim
!>
!> \param[in] incy
!> \verbatim
!>          incy is integer
!>         storage spacing between elements of dy
!> \endverbatim
!>
!> \param[in] c
!> \verbatim
!>          c is double precision
!> \endverbatim
!>
!> \param[in] s
!> \verbatim
!>          s is double precision
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup rot
!
!> \par further details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
subroutine drot(n,dx,incx,dy,incy,c,s)
implicit none
!
!  -- reference blas level1 routine --
!  -- reference blas is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
real(dp), intent(in) :: c !! cosine of the real plane rotation.
real(dp), intent(in) :: s !! sine of the real plane rotation.
integer, intent(in) :: incx !! storage spacing between elements of dx
integer, intent(in) :: incy !! storage spacing between elements of dy
integer, intent(in) :: n !! number of elements in input vector(s)
!     ..
!     .. array arguments ..
real(dp), intent(inout) :: dx(*) !! first strided vector, updated in place when requested by the operation.
real(dp), intent(inout) :: dy(*) !! second strided vector, updated in place when requested by the operation.
!     ..
!
!  =====================================================================
!
!     .. local scalars ..
real(dp) :: dtemp
integer :: i, ix, iy
!     ..
if (n.le.0) return
if (incx.eq.1 .and. incy.eq.1) then
!
!       code for both increments equal to 1
!
do i = 1,n
dtemp = c*dx(i) + s*dy(i)
dy(i) = c*dy(i) - s*dx(i)
dx(i) = dtemp
end do
else
!
!       code for unequal increments or equal increments not equal
!         to 1
!
ix = 1
iy = 1
if (incx.lt.0) ix = (-n+1)*incx + 1
if (incy.lt.0) iy = (-n+1)*incy + 1
do i = 1,n
dtemp = c*dx(ix) + s*dy(iy)
dy(iy) = c*dy(iy) - s*dx(ix)
dx(ix) = dtemp
ix = ix + incx
iy = iy + incy
end do
end if
return
!
!     end of drot
!
end subroutine drot
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dtgex2.f
!> \brief \b dtgex2 swaps adjacent diagonal blocks in an upper (quasi) triangular matrix pair by
!an orthogonal equivalence transformation.
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download dtgex2 + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtgex2.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtgex2.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtgex2.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       subroutine dtgex2( wantq, wantz, n, a, lda, b, ldb, q, ldq, z,
!                          ldz, j1, n1, n2, work, lwork, info )
!
!       .. scalar arguments ..
!       logical            wantq, wantz
!       integer            info, j1, lda, ldb, ldq, ldz, lwork, n, n1, n2
!       ..
!       .. array arguments ..
!       double precision   a( lda, * ), b( ldb, * ), q( ldq, * ),
!      $                   work( * ), z( ldz, * )
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dtgex2 swaps adjacent diagonal blocks (a11, b11) and (a22, b22)
!> of size 1-by-1 or 2-by-2 in an upper (quasi) triangular matrix pair
!> (a, b) by an orthogonal equivalence transformation.
!>
!> (a, b) must be in generalized real schur canonical form (as returned
!> by dgges), i.e. a is block upper triangular with 1-by-1 and 2-by-2
!> diagonal blocks. b is upper triangular.
!>
!> optionally, the matrices q and z of generalized schur vectors are
!> updated.
!>
!>        q(in) * a(in) * z(in)**t = q(out) * a(out) * z(out)**t
!>        q(in) * b(in) * z(in)**t = q(out) * b(out) * z(out)**t
!>
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] wantq
!> \verbatim
!>          wantq is logical
!>          .true. : update the left transformation matrix q;
!>          .false.: do not update q.
!> \endverbatim
!>
!> \param[in] wantz
!> \verbatim
!>          wantz is logical
!>          .true. : update the right transformation matrix z;
!>          .false.: do not update z.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is integer
!>          the order of the matrices a and b. n >= 0.
!> \endverbatim
!>
!> \param[in,out] a
!> \verbatim
!>          a is double precision array, dimensions (lda,n)
!>          on entry, the matrix a in the pair (a, b).
!>          on exit, the updated matrix a.
!> \endverbatim
!>
!> \param[in] lda
!> \verbatim
!>          lda is integer
!>          the leading dimension of the array a. lda >= max(1,n).
!> \endverbatim
!>
!> \param[in,out] b
!> \verbatim
!>          b is double precision array, dimensions (ldb,n)
!>          on entry, the matrix b in the pair (a, b).
!>          on exit, the updated matrix b.
!> \endverbatim
!>
!> \param[in] ldb
!> \verbatim
!>          ldb is integer
!>          the leading dimension of the array b. ldb >= max(1,n).
!> \endverbatim
!>
!> \param[in,out] q
!> \verbatim
!>          q is double precision array, dimension (ldq,n)
!>          on entry, if wantq = .true., the orthogonal matrix q.
!>          on exit, the updated matrix q.
!>          not referenced if wantq = .false..
!> \endverbatim
!>
!> \param[in] ldq
!> \verbatim
!>          ldq is integer
!>          the leading dimension of the array q. ldq >= 1.
!>          if wantq = .true., ldq >= n.
!> \endverbatim
!>
!> \param[in,out] z
!> \verbatim
!>          z is double precision array, dimension (ldz,n)
!>          on entry, if wantz =.true., the orthogonal matrix z.
!>          on exit, the updated matrix z.
!>          not referenced if wantz = .false..
!> \endverbatim
!>
!> \param[in] ldz
!> \verbatim
!>          ldz is integer
!>          the leading dimension of the array z. ldz >= 1.
!>          if wantz = .true., ldz >= n.
!> \endverbatim
!>
!> \param[in] j1
!> \verbatim
!>          j1 is integer
!>          the index to the first block (a11, b11). 1 <= j1 <= n.
!> \endverbatim
!>
!> \param[in] n1
!> \verbatim
!>          n1 is integer
!>          the order of the first block (a11, b11). n1 = 0, 1 or 2.
!> \endverbatim
!>
!> \param[in] n2
!> \verbatim
!>          n2 is integer
!>          the order of the second block (a22, b22). n2 = 0, 1 or 2.
!> \endverbatim
!>
!> \param[out] work
!> \verbatim
!>          work is double precision array, dimension (max(1,lwork)).
!> \endverbatim
!>
!> \param[in] lwork
!> \verbatim
!>          lwork is integer
!>          the dimension of the array work.
!>          lwork >=  max( 1, n*(n2+n1), (n2+n1)*(n2+n1)*2 )
!> \endverbatim
!>
!> \param[out] info
!> \verbatim
!>          info is integer
!>            =0: successful exit
!>            >0: if info = 1, the transformed matrix (a, b) would be
!>                too far from generalized schur form; the blocks are
!>                not swapped and (a, b) and (q, z) are unchanged.
!>                the problem of swapping is too ill-conditioned.
!>            <0: if info = -16: lwork is too small. appropriate value
!>                for lwork is returned in work(1).
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup tgex2
!
!> \par further details:
!  =====================
!>
!>  in the current code both weak and strong stability tests are
!>  performed. the user can omit the strong stability test by changing
!>  the internal logical parameter wands to .false.. see ref. [2] for
!>  details.
!
!> \par contributors:
!  ==================
!>
!>     bo kagstrom and peter poromaa, department of computing science,
!>     umea university, s-901 87 umea, sweden.
!
!> \par references:
!  ================
!>
!> \verbatim
!>
!>  [1] b. kagstrom; a direct method for reordering eigenvalues in the
!>      generalized real schur form of a regular matrix pair (a, b), in
!>      m.s. moonen et al (eds), linear algebra for large scale and
!>      real-time applications, kluwer academic publ. 1993, pp 195-218.
!>
!>  [2] b. kagstrom and p. poromaa; computing eigenspaces with specified
!>      eigenvalues of a regular matrix pair (a, b) and condition
!>      estimation: theory, algorithms and software,
!>      report uminf - 94.04, department of computing science, umea
!>      university, s-901 87 umea, sweden, 1994. also as lapack working
!>      note 87. to appear in numerical algorithms, 1996.
!> \endverbatim
!>
!  =====================================================================
subroutine dtgex2( wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, j1, n1, n2, work, lwork, info )
implicit none
!
!  -- lapack auxiliary routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
logical, intent(in) :: wantq !! .true. : update the left transformation matrix q; .false.: do not
!! update q.
logical, intent(in) :: wantz !! .true. : update the right transformation matrix z; .false.: do not
!! update z.
integer, intent(out) :: info !! zero success; negative invalid argument or positive numerical failure.
integer, intent(in) :: j1 !! the index to the first block (a11, b11). 1 <= j1 <= n.
integer, intent(in) :: lda !! the leading dimension of the array a. lda >= max(1,n).
integer, intent(in) :: ldb !! the leading dimension of the array b. ldb >= max(1,n).
integer, intent(in) :: ldq !! the leading dimension of the array q. ldq >= 1. if wantq = .true.,
!! ldq >= n.
integer, intent(in) :: ldz !! the leading dimension of the array z. ldz >= 1. if wantz = .true.,
!! ldz >= n.
integer, intent(in) :: lwork !! the dimension of the array work. lwork >=  max( 1, n*(n2+n1),
!! (n2+n1)*(n2+n1)*2 )
integer, intent(in) :: n !! the order of the matrices a and b. n >= 0.
integer, intent(in) :: n1 !! the order of the first block (a11, b11). n1 = 0, 1 or 2.
integer, intent(in) :: n2 !! the order of the second block (a22, b22). n2 = 0, 1 or 2.
!     ..
!     .. array arguments ..
real(dp), intent(inout) :: a( lda, * ) !! on entry, the matrix a in the pair (a, b). on exit, the
!! updated matrix a.
real(dp), intent(inout) :: b( ldb, * ) !! on entry, the matrix b in the pair (a, b). on exit, the
!! updated matrix b.
real(dp), intent(inout) :: q( ldq, * ) !! on entry, if wantq = .true., the orthogonal matrix q. on
!! exit, the updated matrix q. not referenced if wantq =
!! .false..
real(dp), intent(out) :: work( * ) !! caller-owned scratch array, sized according to the retained contract.
real(dp), intent(inout) :: z( ldz, * ) !! on entry, if wantz =.true., the orthogonal matrix z. on
!! exit, the updated matrix z. not referenced if wantz =
!! .false..
!     ..
!
!  =====================================================================
!  replaced various illegal calls to dcopy by calls to dlaset, or by do
!  loops. sven hammarling, 1/5/02.
!
!     .. parameters ..
real(dp) :: zero, one
parameter          ( zero = 0.0d+0, one = 1.0d+0 )
real(dp) :: twenty
parameter          ( twenty = 2.0d+01 )
integer :: ldst
parameter          ( ldst = 4 )
logical :: wands
parameter          ( wands = .true. )
!     ..
!     .. local scalars ..
logical :: strong, weak
integer :: i, idum, linfo, m
real(dp) :: bqra21, brqa21, ddum, dnorma, dnormb, dscale, dsum, eps, f, g, sa, sb, scale, smlnum, thresha, &
    threshb
!     ..
!     .. local arrays ..
integer :: iwork( ldst + 2 )
real(dp) :: ai( 2 ), ar( 2 ), be( 2 ), ir( ldst, ldst ), ircop( ldst, ldst ), li( ldst, ldst ), licop( &
    ldst, ldst ), s( ldst, ldst ), scpy( ldst, ldst ), t( ldst, ldst ), taul( ldst ), taur( ldst ), tcpy( &
    ldst, ldst )
!     ..
!     .. external functions ..
!     ..
!     .. external subroutines ..
!     ..
!     .. intrinsic functions ..
intrinsic          abs, max, sqrt
!     ..
!     .. executable statements ..
!
info = 0
!
!     quick return if possible
!
if( n.le.1 .or. n1.le.0 .or. n2.le.0 ) return
if( n1.gt.n .or. ( j1+n1 ).gt.n ) return
m = n1 + n2
if( lwork.lt.max( 1, n*m, m*m*2 ) ) then
info = -16
work( 1 ) = max( 1, n*m, m*m*2 )
return
end if
!
weak = .false.
strong = .false.
!
!     make a local copy of selected block
!
call dlaset( 'full', ldst, ldst, zero, zero, li, ldst )
call dlaset( 'full', ldst, ldst, zero, zero, ir, ldst )
call dlacpy( 'full', m, m, a( j1, j1 ), lda, s, ldst )
call dlacpy( 'full', m, m, b( j1, j1 ), ldb, t, ldst )
!
!     compute threshold for testing acceptance of swapping.
!
eps = dlamch( 'p' )
smlnum = dlamch( 's' ) / eps
dscale = zero
dsum = one
call dlacpy( 'full', m, m, s, ldst, work, m )
call dlassq( m*m, work, 1, dscale, dsum )
dnorma = dscale*sqrt( dsum )
dscale = zero
dsum = one
call dlacpy( 'full', m, m, t, ldst, work, m )
call dlassq( m*m, work, 1, dscale, dsum )
dnormb = dscale*sqrt( dsum )
!
!     thres has been changed from
!        thresh = max( ten*eps*sa, smlnum )
!     to
!        thresh = max( twenty*eps*sa, smlnum )
!     on 04/01/10.
!     "bug" reported by ondra kamenik, confirmed by julie langou, fixed by
!     jim demmel and guillaume revy. see forum post 1783.
!
thresha = max( twenty*eps*dnorma, smlnum )
threshb = max( twenty*eps*dnormb, smlnum )
!
if( m.eq.2 ) then
!
!        case 1: swap 1-by-1 and 1-by-1 blocks.
!
!        compute orthogonal ql and rq that swap 1-by-1 and 1-by-1 blocks
!        using givens rotations and perform the swap tentatively.
!
f = s( 2, 2 )*t( 1, 1 ) - t( 2, 2 )*s( 1, 1 )
g = s( 2, 2 )*t( 1, 2 ) - t( 2, 2 )*s( 1, 2 )
sa = abs( s( 2, 2 ) ) * abs( t( 1, 1 ) )
sb = abs( s( 1, 1 ) ) * abs( t( 2, 2 ) )
call dlartg( f, g, ir( 1, 2 ), ir( 1, 1 ), ddum )
ir( 2, 1 ) = -ir( 1, 2 )
ir( 2, 2 ) = ir( 1, 1 )
call drot( 2, s( 1, 1 ), 1, s( 1, 2 ), 1, ir( 1, 1 ), ir( 2, 1 ) )
call drot( 2, t( 1, 1 ), 1, t( 1, 2 ), 1, ir( 1, 1 ), ir( 2, 1 ) )
if( sa.ge.sb ) then
call dlartg( s( 1, 1 ), s( 2, 1 ), li( 1, 1 ), li( 2, 1 ), ddum )
else
call dlartg( t( 1, 1 ), t( 2, 1 ), li( 1, 1 ), li( 2, 1 ), ddum )
end if
call drot( 2, s( 1, 1 ), ldst, s( 2, 1 ), ldst, li( 1, 1 ), li( 2, 1 ) )
call drot( 2, t( 1, 1 ), ldst, t( 2, 1 ), ldst, li( 1, 1 ), li( 2, 1 ) )
li( 2, 2 ) = li( 1, 1 )
li( 1, 2 ) = -li( 2, 1 )
!
!        weak stability test: |s21| <= o(eps f-norm((a)))
!                           and  |t21| <= o(eps f-norm((b)))
!
weak = abs( s( 2, 1 ) ) .le. thresha .and. abs( t( 2, 1 ) ) .le. threshb
if( .not.weak ) go to 70
!
if( wands ) then
!
!           strong stability test:
!               f-norm((a-ql**h*s*qr)) <= o(eps*f-norm((a)))
!               and
!               f-norm((b-ql**h*t*qr)) <= o(eps*f-norm((b)))
!
call dlacpy( 'full', m, m, a( j1, j1 ), lda, work( m*m+1 ), m )
call dgemm( 'n', 'n', m, m, m, one, li, ldst, s, ldst, zero, work, m )
call dgemm( 'n', 't', m, m, m, -one, work, m, ir, ldst, one, work( m*m+1 ), m )
dscale = zero
dsum = one
call dlassq( m*m, work( m*m+1 ), 1, dscale, dsum )
sa = dscale*sqrt( dsum )
!
call dlacpy( 'full', m, m, b( j1, j1 ), ldb, work( m*m+1 ), m )
call dgemm( 'n', 'n', m, m, m, one, li, ldst, t, ldst, zero, work, m )
call dgemm( 'n', 't', m, m, m, -one, work, m, ir, ldst, one, work( m*m+1 ), m )
dscale = zero
dsum = one
call dlassq( m*m, work( m*m+1 ), 1, dscale, dsum )
sb = dscale*sqrt( dsum )
strong = sa.le.thresha .and. sb.le.threshb
if( .not.strong ) go to 70
end if
!
!        update (a(j1:j1+m-1, m+j1:n), b(j1:j1+m-1, m+j1:n)) and
!               (a(1:j1-1, j1:j1+m), b(1:j1-1, j1:j1+m)).
!
call drot( j1+1, a( 1, j1 ), 1, a( 1, j1+1 ), 1, ir( 1, 1 ), ir( 2, 1 ) )
call drot( j1+1, b( 1, j1 ), 1, b( 1, j1+1 ), 1, ir( 1, 1 ), ir( 2, 1 ) )
call drot( n-j1+1, a( j1, j1 ), lda, a( j1+1, j1 ), lda, li( 1, 1 ), li( 2, 1 ) )
call drot( n-j1+1, b( j1, j1 ), ldb, b( j1+1, j1 ), ldb, li( 1, 1 ), li( 2, 1 ) )
!
!        set  n1-by-n2 (2,1) - blocks to zero.
!
a( j1+1, j1 ) = zero
b( j1+1, j1 ) = zero
!
!        accumulate transformations into q and z if requested.
!
if( wantz ) call drot( n, z( 1, j1 ), 1, z( 1, j1+1 ), 1, ir( 1, 1 ), ir( 2, 1 ) )
if( wantq ) call drot( n, q( 1, j1 ), 1, q( 1, j1+1 ), 1, li( 1, 1 ), li( 2, 1 ) )
!
!        exit with info = 0 if swap was successfully performed.
!
return
!
else
!
!        case 2: swap 1-by-1 and 2-by-2 blocks, or 2-by-2
!                and 2-by-2 blocks.
!
!        solve the generalized sylvester equation
!                 s11 * r - l * s22 = scale * s12
!                 t11 * r - l * t22 = scale * t12
!        for r and l. solutions in li and ir.
!
call dlacpy( 'full', n1, n2, t( 1, n1+1 ), ldst, li, ldst )
call dlacpy( 'full', n1, n2, s( 1, n1+1 ), ldst, ir( n2+1, n1+1 ), ldst )
call dtgsy2( 'n', 0, n1, n2, s, ldst, s( n1+1, n1+1 ), ldst, ir( n2+1, n1+1 ), ldst, t, ldst, t( n1+1, &
    n1+1 ), ldst, li, ldst, scale, dsum, dscale, iwork, idum, linfo )
if( linfo.ne.0 ) go to 70
!
!        compute orthogonal matrix ql:
!
!                    ql**t * li = [ tl ]
!                                 [ 0  ]
!        where
!                    li =  [      -l              ]
!                          [ scale * identity(n2) ]
!
do 10 i = 1, n2
call dscal( n1, -one, li( 1, i ), 1 )
li( n1+i, i ) = scale
10 continue
call dgeqr2( m, n2, li, ldst, taul, work, linfo )
if( linfo.ne.0 ) go to 70
call dorg2r( m, m, n2, li, ldst, taul, work, linfo )
if( linfo.ne.0 ) go to 70
!
!        compute orthogonal matrix rq:
!
!                    ir * rq**t =   [ 0  tr],
!
!         where ir = [ scale * identity(n1), r ]
!
do 20 i = 1, n1
ir( n2+i, i ) = scale
20 continue
call dgerq2( n1, m, ir( n2+1, 1 ), ldst, taur, work, linfo )
if( linfo.ne.0 ) go to 70
call dorgr2( m, m, n1, ir, ldst, taur, work, linfo )
if( linfo.ne.0 ) go to 70
!
!        perform the swapping tentatively:
!
call dgemm( 't', 'n', m, m, m, one, li, ldst, s, ldst, zero, work, m )
call dgemm( 'n', 't', m, m, m, one, work, m, ir, ldst, zero, s, ldst )
call dgemm( 't', 'n', m, m, m, one, li, ldst, t, ldst, zero, work, m )
call dgemm( 'n', 't', m, m, m, one, work, m, ir, ldst, zero, t, ldst )
call dlacpy( 'f', m, m, s, ldst, scpy, ldst )
call dlacpy( 'f', m, m, t, ldst, tcpy, ldst )
call dlacpy( 'f', m, m, ir, ldst, ircop, ldst )
call dlacpy( 'f', m, m, li, ldst, licop, ldst )
!
!        triangularize the b-part by an rq factorization.
!        apply transformation (from left) to a-part, giving s.
!
call dgerq2( m, m, t, ldst, taur, work, linfo )
if( linfo.ne.0 ) go to 70
call dormr2( 'r', 't', m, m, m, t, ldst, taur, s, ldst, work, linfo )
if( linfo.ne.0 ) go to 70
call dormr2( 'l', 'n', m, m, m, t, ldst, taur, ir, ldst, work, linfo )
if( linfo.ne.0 ) go to 70
!
!        compute f-norm(s21) in brqa21. (t21 is 0.)
!
dscale = zero
dsum = one
do 30 i = 1, n2
call dlassq( n1, s( n2+1, i ), 1, dscale, dsum )
30 continue
brqa21 = dscale*sqrt( dsum )
!
!        triangularize the b-part by a qr factorization.
!        apply transformation (from right) to a-part, giving s.
!
call dgeqr2( m, m, tcpy, ldst, taul, work, linfo )
if( linfo.ne.0 ) go to 70
call dorm2r( 'l', 't', m, m, m, tcpy, ldst, taul, scpy, ldst, work, info )
call dorm2r( 'r', 'n', m, m, m, tcpy, ldst, taul, licop, ldst, work, info )
if( linfo.ne.0 ) go to 70
!
!        compute f-norm(s21) in bqra21. (t21 is 0.)
!
dscale = zero
dsum = one
do 40 i = 1, n2
call dlassq( n1, scpy( n2+1, i ), 1, dscale, dsum )
40 continue
bqra21 = dscale*sqrt( dsum )
!
!        decide which method to use.
!          weak stability test:
!             f-norm(s21) <= o(eps * f-norm((s)))
!
if( bqra21.le.brqa21 .and. bqra21.le.thresha ) then
call dlacpy( 'f', m, m, scpy, ldst, s, ldst )
call dlacpy( 'f', m, m, tcpy, ldst, t, ldst )
call dlacpy( 'f', m, m, ircop, ldst, ir, ldst )
call dlacpy( 'f', m, m, licop, ldst, li, ldst )
else if( brqa21.ge.thresha ) then
go to 70
end if
!
!        set lower triangle of b-part to zero
!
call dlaset( 'lower', m-1, m-1, zero, zero, t(2,1), ldst )
!
if( wands ) then
!
!           strong stability test:
!               f-norm((a-ql**h*s*qr)) <= o(eps*f-norm((a)))
!               and
!               f-norm((b-ql**h*t*qr)) <= o(eps*f-norm((b)))
!
call dlacpy( 'full', m, m, a( j1, j1 ), lda, work( m*m+1 ), m )
call dgemm( 'n', 'n', m, m, m, one, li, ldst, s, ldst, zero, work, m )
call dgemm( 'n', 'n', m, m, m, -one, work, m, ir, ldst, one, work( m*m+1 ), m )
dscale = zero
dsum = one
call dlassq( m*m, work( m*m+1 ), 1, dscale, dsum )
sa = dscale*sqrt( dsum )
!
call dlacpy( 'full', m, m, b( j1, j1 ), ldb, work( m*m+1 ), m )
call dgemm( 'n', 'n', m, m, m, one, li, ldst, t, ldst, zero, work, m )
call dgemm( 'n', 'n', m, m, m, -one, work, m, ir, ldst, one, work( m*m+1 ), m )
dscale = zero
dsum = one
call dlassq( m*m, work( m*m+1 ), 1, dscale, dsum )
sb = dscale*sqrt( dsum )
strong = sa.le.thresha .and. sb.le.threshb
if( .not.strong ) go to 70
!
end if
!
!        if the swap is accepted ("weakly" and "strongly"), apply the
!        transformations and set n1-by-n2 (2,1)-block to zero.
!
call dlaset( 'full', n1, n2, zero, zero, s(n2+1,1), ldst )
!
!        copy back m-by-m diagonal block starting at index j1 of (a, b)
!
call dlacpy( 'f', m, m, s, ldst, a( j1, j1 ), lda )
call dlacpy( 'f', m, m, t, ldst, b( j1, j1 ), ldb )
call dlaset( 'full', ldst, ldst, zero, zero, t, ldst )
!
!        standardize existing 2-by-2 blocks.
!
call dlaset( 'full', m, m, zero, zero, work, m )
work( 1 ) = one
t( 1, 1 ) = one
idum = lwork - m*m - 2
if( n2.gt.1 ) then
call dlagv2( a( j1, j1 ), lda, b( j1, j1 ), ldb, ar, ai, be, work( 1 ), work( 2 ), t( 1, 1 ), t( 2, 1 ) )
work( m+1 ) = -work( 2 )
work( m+2 ) = work( 1 )
t( n2, n2 ) = t( 1, 1 )
t( 1, 2 ) = -t( 2, 1 )
end if
work( m*m ) = one
t( m, m ) = one
!
if( n1.gt.1 ) then
call dlagv2( a( j1+n2, j1+n2 ), lda, b( j1+n2, j1+n2 ), ldb, taur, taul, work( m*m+1 ), work( n2*m+n2+1 ), &
    work( n2*m+n2+2 ), t( n2+1, n2+1 ), t( m, m-1 ) )
work( m*m ) = work( n2*m+n2+1 )
work( m*m-1 ) = -work( n2*m+n2+2 )
t( m, m ) = t( n2+1, n2+1 )
t( m-1, m ) = -t( m, m-1 )
end if
call dgemm( 't', 'n', n2, n1, n2, one, work, m, a( j1, j1+n2 ), lda, zero, work( m*m+1 ), n2 )
call dlacpy( 'full', n2, n1, work( m*m+1 ), n2, a( j1, j1+n2 ), lda )
call dgemm( 't', 'n', n2, n1, n2, one, work, m, b( j1, j1+n2 ), ldb, zero, work( m*m+1 ), n2 )
call dlacpy( 'full', n2, n1, work( m*m+1 ), n2, b( j1, j1+n2 ), ldb )
call dgemm( 'n', 'n', m, m, m, one, li, ldst, work, m, zero, work( m*m+1 ), m )
call dlacpy( 'full', m, m, work( m*m+1 ), m, li, ldst )
call dgemm( 'n', 'n', n2, n1, n1, one, a( j1, j1+n2 ), lda, t( n2+1, n2+1 ), ldst, zero, work, n2 )
call dlacpy( 'full', n2, n1, work, n2, a( j1, j1+n2 ), lda )
call dgemm( 'n', 'n', n2, n1, n1, one, b( j1, j1+n2 ), ldb, t( n2+1, n2+1 ), ldst, zero, work, n2 )
call dlacpy( 'full', n2, n1, work, n2, b( j1, j1+n2 ), ldb )
call dgemm( 't', 'n', m, m, m, one, ir, ldst, t, ldst, zero, work, m )
call dlacpy( 'full', m, m, work, m, ir, ldst )
!
!        accumulate transformations into q and z if requested.
!
if( wantq ) then
call dgemm( 'n', 'n', n, m, m, one, q( 1, j1 ), ldq, li, ldst, zero, work, n )
call dlacpy( 'full', n, m, work, n, q( 1, j1 ), ldq )
!
end if
!
if( wantz ) then
call dgemm( 'n', 'n', n, m, m, one, z( 1, j1 ), ldz, ir, ldst, zero, work, n )
call dlacpy( 'full', n, m, work, n, z( 1, j1 ), ldz )
!
end if
!
!        update (a(j1:j1+m-1, m+j1:n), b(j1:j1+m-1, m+j1:n)) and
!                (a(1:j1-1, j1:j1+m), b(1:j1-1, j1:j1+m)).
!
i = j1 + m
if( i.le.n ) then
call dgemm( 't', 'n', m, n-i+1, m, one, li, ldst, a( j1, i ), lda, zero, work, m )
call dlacpy( 'full', m, n-i+1, work, m, a( j1, i ), lda )
call dgemm( 't', 'n', m, n-i+1, m, one, li, ldst, b( j1, i ), ldb, zero, work, m )
call dlacpy( 'full', m, n-i+1, work, m, b( j1, i ), ldb )
end if
i = j1 - 1
if( i.gt.0 ) then
call dgemm( 'n', 'n', i, m, m, one, a( 1, j1 ), lda, ir, ldst, zero, work, i )
call dlacpy( 'full', i, m, work, i, a( 1, j1 ), lda )
call dgemm( 'n', 'n', i, m, m, one, b( 1, j1 ), ldb, ir, ldst, zero, work, i )
call dlacpy( 'full', i, m, work, i, b( 1, j1 ), ldb )
end if
!
!        exit with info = 0 if swap was successfully performed.
!
return
!
end if
!
!     exit with info = 1 if swap was rejected.
!
70 continue
!
info = 1
return
!
!     end of dtgex2
!
end subroutine dtgex2
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlassq.f90
!> \brief \b dlassq updates a sum of squares represented in scaled form.
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download dlassq + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlassq.f90">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlassq.f90">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlassq.f90">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       subroutine dlassq( n, x, incx, scale, sumsq )
!
!       .. scalar arguments ..
!       integer            incx, n
!       double precision   scale, sumsq
!       ..
!       .. array arguments ..
!       double precision   x( * )
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dlassq returns the values scale_out and sumsq_out such that
!>
!>    (scale_out**2)*sumsq_out = x( 1 )**2 +...+ x( n )**2 + (scale**2)*sumsq,
!>
!> where x( i ) = x( 1 + ( i - 1 )*incx ). the value of sumsq is
!> assumed to be non-negative.
!>
!> scale and sumsq must be supplied in scale and sumsq and
!> scale_out and sumsq_out are overwritten on scale and sumsq respectively.
!>
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] n
!> \verbatim
!>          n is integer
!>          the number of elements to be used from the vector x.
!> \endverbatim
!>
!> \param[in] x
!> \verbatim
!>          x is double precision array, dimension (1+(n-1)*abs(incx))
!>          the vector for which a scaled sum of squares is computed.
!>             x( i ) = x( 1 + ( i - 1 )*incx ), 1 <= i <= n.
!> \endverbatim
!>
!> \param[in] incx
!> \verbatim
!>          incx is integer
!>          the increment between successive values of the vector x.
!>          if incx > 0, x(1+(i-1)*incx) = x(i) for 1 <= i <= n
!>          if incx < 0, x(1-(n-i)*incx) = x(i) for 1 <= i <= n
!>          if incx = 0, x isn't a vector so there is no need to call
!>          this subroutine. if you call it anyway, it will count x(1)
!>          in the vector norm n times.
!> \endverbatim
!>
!> \param[in,out] scale
!> \verbatim
!>          scale is double precision
!>          on entry, the value scale in the equation above.
!>          on exit, scale is overwritten by scale_out, the scaling factor
!>          for the sum of squares.
!> \endverbatim
!>
!> \param[in,out] sumsq
!> \verbatim
!>          sumsq is double precision
!>          on entry, the value sumsq in the equation above.
!>          on exit, sumsq is overwritten by sumsq_out, the basic sum of
!>          squares from which scale_out has been factored out.
!> \endverbatim
!
!  authors:
!  ========
!
!> \author edward anderson, lockheed martin
!
!> \par contributors:
!  ==================
!>
!> weslley pereira, university of colorado denver, usa
!> nick papior, technical university of denmark, dk
!
!> \par further details:
!  =====================
!>
!> \verbatim
!>
!>  anderson e. (2017)
!>  algorithm 978: safe scaling in the level 1 blas
!>  acm trans math softw 44:1--28
!>  https://doi.org/10.1145/3061665
!>
!>  blue, james l. (1978)
!>  a portable fortran program to find the euclidean norm of a vector
!>  acm trans math softw 4:15--23
!>  https://doi.org/10.1145/355769.355771
!>
!> \endverbatim
!
!> \ingroup lassq
!
!  =====================================================================
subroutine dlassq( n, x, incx, scale, sumsq )
implicit none
!
!  -- lapack auxiliary routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!  .. scalar arguments ..
integer, intent(in) :: incx !! the increment between successive values of the vector x. if incx >
!! 0, x(1+(i-1)*incx) = x(i) for 1 <= i <= n
integer, intent(in) :: n !! the number of elements to be used from the vector x.
real(dp), intent(inout) :: scale !! on entry, the value scale in the equation above. on exit,
!! scale is overwritten by scale_out, the scaling factor
real(dp), intent(inout) :: sumsq !! on entry, the value sumsq in the equation above. on exit,
!! sumsq is overwritten by sumsq_out, the basic sum of
!  ..
!  .. array arguments ..
real(dp), intent(in) :: x(*) !! the vector for which a scaled sum of squares is computed. x( i ) =
!! x( 1 + ( i - 1 )*incx ), 1 <= i <= n.
!  ..
!  .. local scalars ..
integer :: i, ix
logical :: notbig
real(dp) :: abig, amed, asml, ax, ymax, ymin
!  ..
!
!  quick return if possible
!
if( ieee_is_nan(scale) .or. ieee_is_nan(sumsq) ) return
if( sumsq == zero ) scale = one
if( scale == zero ) then
scale = one
sumsq = zero
end if
if (n <= 0) then
return
end if
!
!  compute the sum of squares in 3 accumulators:
!     abig -- sums of squares scaled down to avoid overflow
!     asml -- sums of squares scaled up to avoid underflow
!     amed -- sums of squares that do not require scaling
!  the thresholds and multipliers are
!     tbig -- values bigger than this are scaled down by sbig
!     tsml -- values smaller than this are scaled up by ssml
!
notbig = .true.
asml = zero
amed = zero
abig = zero
ix = 1
if( incx < 0 ) ix = 1 - (n-1)*incx
do i = 1, n
ax = abs(x(ix))
if (ax > tbig) then
abig = abig + (ax*sbig)**2
notbig = .false.
else if (ax < tsml) then
if (notbig) asml = asml + (ax*ssml)**2
else
amed = amed + ax**2
end if
ix = ix + incx
end do
!
!  put the existing sum of squares into one of the accumulators
!
if( sumsq > zero ) then
ax = scale*sqrt( sumsq )
if (ax > tbig) then
if (scale > one) then
scale = scale * sbig
abig = abig + scale * (scale * sumsq)
else
! sumsq > tbig^2 => (sbig * (sbig * sumsq)) is representable
abig = abig + scale * (scale * (sbig * (sbig * sumsq)))
end if
else if (ax < tsml) then
if (notbig) then
if (scale < one) then
scale = scale * ssml
asml = asml + scale * (scale * sumsq)
else
! sumsq < tsml^2 => (ssml * (ssml * sumsq)) is representable
asml = asml + scale * (scale * (ssml * (ssml * sumsq)))
end if
end if
else
amed = amed + scale * (scale * sumsq)
end if
end if
!
!  combine abig and amed or amed and asml if more than one
!  accumulator was used.
!
if (abig > zero) then
!
!     combine abig and amed if abig > 0.
!
if (amed > zero .or. ieee_is_nan(amed)) then
abig = abig + (amed*sbig)*sbig
end if
scale = one / sbig
sumsq = abig
else if (asml > zero) then
!
!     combine amed and asml if asml > 0.
!
if (amed > zero .or. ieee_is_nan(amed)) then
amed = sqrt(amed)
asml = sqrt(asml) / ssml
if (asml > amed) then
ymin = amed
ymax = asml
else
ymin = asml
ymax = amed
end if
scale = one
sumsq = ymax**2*( one + (ymin/ymax)**2 )
else
scale = one / ssml
sumsq = asml
end if
else
!
!     otherwise all values are mid-range or zero
!
scale = one
sumsq = amed
end if
return
end subroutine dlassq
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/disnan.f
!> \brief \b disnan tests input for nan.
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download disnan + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/disnan.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/disnan.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/disnan.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       logical function disnan( din )
!
!       .. scalar arguments ..
!       double precision, intent(in) :: din
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> disnan returns .true. if its argument is nan, and .false.
!> otherwise.  to be replaced by the fortran 2003 intrinsic in the
!> future.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] din
!> \verbatim
!>          din is double precision
!>          input to test for nan.
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup isnan
!
!  =====================================================================
logical function disnan( din ) result(lapack_result)
implicit none
!
!  -- lapack auxiliary routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
real(dp), intent(in) :: din
!     ..
!
!  =====================================================================
!
!  .. external functions ..
!  ..
!  .. executable statements ..
lapack_result = dlaisnan(din,din)
return
end function disnan
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/BLAS/SRC/dnrm2.f90
!> \brief \b dnrm2
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  definition:
!  ===========
!
!       double precision function dnrm2(n,x,incx)
!
!       .. scalar arguments ..
!       integer incx,n
!       ..
!       .. array arguments ..
!       double precision x(*)
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dnrm2 returns the euclidean norm of a vector via the function
!> name, so that
!>
!>    dnrm2 := sqrt( x'*x )
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] n
!> \verbatim
!>          n is integer
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in] x
!> \verbatim
!>          x is double precision array, dimension ( 1 + ( n - 1 )*abs( incx ) )
!> \endverbatim
!>
!> \param[in] incx
!> \verbatim
!>          incx is integer, storage spacing between elements of x
!>          if incx > 0, x(1+(i-1)*incx) = x(i) for 1 <= i <= n
!>          if incx < 0, x(1-(n-i)*incx) = x(i) for 1 <= i <= n
!>          if incx = 0, x isn't a vector so there is no need to call
!>          this subroutine.  if you call it anyway, it will count x(1)
!>          in the vector norm n times.
!> \endverbatim
!
!  authors:
!  ========
!
!> \author edward anderson, lockheed martin
!
!> \date august 2016
!
!> \ingroup nrm2
!
!> \par contributors:
!  ==================
!>
!> weslley pereira, university of colorado denver, usa
!
!> \par further details:
!  =====================
!>
!> \verbatim
!>
!>  anderson e. (2017)
!>  algorithm 978: safe scaling in the level 1 blas
!>  acm trans math softw 44:1--28
!>  https://doi.org/10.1145/3061665
!>
!>  blue, james l. (1978)
!>  a portable fortran program to find the euclidean norm of a vector
!>  acm trans math softw 4:15--23
!>  https://doi.org/10.1145/355769.355771
!>
!> \endverbatim
!>
!  =====================================================================
function dnrm2( n, x, incx ) result(lapack_result)
implicit none
real(dp) :: lapack_result
!
!  -- reference blas level1 routine (version 3.9.1) --
!  -- reference blas is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!     march 2021
!
!  .. constants ..
real(dp), parameter :: zero = 0.0_dp
real(dp), parameter :: one  = 1.0_dp
real(dp), parameter :: maxn = huge(0.0_dp)
!  ..
!  .. blue's scaling constants ..
real(dp), parameter :: tsml = real(radix(0._dp), dp)**ceiling( (minexponent(0._dp) - 1) * 0.5_dp)
real(dp), parameter :: tbig = real(radix(0._dp), dp)**floor( (maxexponent(0._dp) - digits(0._dp) + 1) * &
    0.5_dp)
real(dp), parameter :: ssml = real(radix(0._dp), dp)**( - floor( (minexponent(0._dp) - digits(0._dp)) * &
    0.5_dp))
real(dp), parameter :: sbig = real(radix(0._dp), dp)**( - ceiling( (maxexponent(0._dp) + digits(0._dp) - &
    1) * 0.5_dp))
!  ..
!  .. scalar arguments ..
integer, intent(in) :: incx !! if incx > 0, x(1+(i-1)*incx) = x(i) for 1 <= i <= n if incx < 0,
!! x(1-(n-i)*incx) = x(i) for 1 <= i <= n
integer, intent(in) :: n !! number of elements in input vector(s)
!  ..
!  .. array arguments ..
real(dp), intent(in) :: x(*) !! input component or strided vector used by this norm or nan check.
!  ..
!  .. local scalars ..
integer :: i, ix
logical :: notbig
real(dp) :: abig, amed, asml, ax, scl, sumsq, ymax, ymin
!
!  quick return if possible
!
lapack_result = zero
if( n <= 0 ) return
!
scl = one
sumsq = zero
!
!  compute the sum of squares in 3 accumulators:
!     abig -- sums of squares scaled down to avoid overflow
!     asml -- sums of squares scaled up to avoid underflow
!     amed -- sums of squares that do not require scaling
!  the thresholds and multipliers are
!     tbig -- values bigger than this are scaled down by sbig
!     tsml -- values smaller than this are scaled up by ssml
!
notbig = .true.
asml = zero
amed = zero
abig = zero
ix = 1
if( incx < 0 ) ix = 1 - (n-1)*incx
do i = 1, n
ax = abs(x(ix))
if (ax > tbig) then
abig = abig + (ax*sbig)**2
notbig = .false.
else if (ax < tsml) then
if (notbig) asml = asml + (ax*ssml)**2
else
amed = amed + ax**2
end if
ix = ix + incx
end do
!
!  combine abig and amed or amed and asml if more than one
!  accumulator was used.
!
if (abig > zero) then
!
!     combine abig and amed if abig > 0.
!
if ( (amed > zero) .or. (amed > maxn) .or. (amed /= amed) ) then
abig = abig + (amed*sbig)*sbig
end if
scl = one / sbig
sumsq = abig
else if (asml > zero) then
!
!     combine amed and asml if asml > 0.
!
if ( (amed > zero) .or. (amed > maxn) .or. (amed /= amed) ) then
amed = sqrt(amed)
asml = sqrt(asml) / ssml
if (asml > amed) then
ymin = amed
ymax = asml
else
ymin = asml
ymax = amed
end if
scl = one
sumsq = ymax**2*( one + (ymin/ymax)**2 )
else
scl = one / ssml
sumsq = asml
end if
else
!
!     otherwise all values are mid-range
!
scl = one
sumsq = amed
end if
lapack_result = scl*sqrt( sumsq )
return
end function dnrm2
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/BLAS/SRC/dscal.f
!> \brief \b dscal
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  definition:
!  ===========
!
!       subroutine dscal(n,da,dx,incx)
!
!       .. scalar arguments ..
!       double precision da
!       integer incx,n
!       ..
!       .. array arguments ..
!       double precision dx(*)
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!>    dscal scales a vector by a constant.
!>    uses unrolled loops for increment equal to 1.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] n
!> \verbatim
!>          n is integer
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in] da
!> \verbatim
!>          da is double precision
!>           on entry, da specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in,out] dx
!> \verbatim
!>          dx is double precision array, dimension ( 1 + ( n - 1 )*abs( incx ) )
!> \endverbatim
!>
!> \param[in] incx
!> \verbatim
!>          incx is integer
!>         storage spacing between elements of dx
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup scal
!
!> \par further details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 3/93 to return if incx .le. 0.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
subroutine dscal(n,da,dx,incx)
implicit none
!
!  -- reference blas level1 routine --
!  -- reference blas is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
real(dp), intent(in) :: da !! on entry, da specifies the scalar alpha.
integer, intent(in) :: incx !! storage spacing between elements of dx
integer, intent(in) :: n !! number of elements in input vector(s)
!     ..
!     .. array arguments ..
real(dp), intent(inout) :: dx(*) !! first strided vector, updated in place when requested by the operation.
!     ..
!
!  =====================================================================
!
!     .. local scalars ..
integer :: i, m, mp1, nincx
!     .. parameters ..
real(dp) :: one
parameter (one=1.0d+0)
!     ..
!     .. intrinsic functions ..
intrinsic mod
!     ..
if (n.le.0 .or. incx.le.0 .or. da.eq.one) return
if (incx.eq.1) then
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
m = mod(n,5)
if (m.ne.0) then
do i = 1,m
dx(i) = da*dx(i)
end do
if (n.lt.5) return
end if
mp1 = m + 1
do i = mp1,n,5
dx(i) = da*dx(i)
dx(i+1) = da*dx(i+1)
dx(i+2) = da*dx(i+2)
dx(i+3) = da*dx(i+3)
dx(i+4) = da*dx(i+4)
end do
else
!
!        code for increment not equal to 1
!
nincx = n*incx
do i = 1,nincx,incx
dx(i) = da*dx(i)
end do
end if
return
!
!     end of dscal
!
end subroutine dscal
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/BLAS/SRC/dgemm.f
!> \brief \b dgemm
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  definition:
!  ===========
!
!       subroutine dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
!
!       .. scalar arguments ..
!       double precision alpha,beta
!       integer k,lda,ldb,ldc,m,n
!       character transa,transb
!       ..
!       .. array arguments ..
!       double precision a(lda,*),b(ldb,*),c(ldc,*)
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dgemm  performs one of the matrix-matrix operations
!>
!>    c := alpha*op( a )*op( b ) + beta*c,
!>
!> where  op( x ) is one of
!>
!>    op( x ) = x   or   op( x ) = x**t,
!>
!> alpha and beta are scalars, and a, b and c are matrices, with op( a )
!> an m by k matrix,  op( b )  a  k by n matrix and  c an m by n matrix.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] transa
!> \verbatim
!>          transa is character*1
!>           on entry, transa specifies the form of op( a ) to be used in
!>           the matrix multiplication as follows:
!>
!>              transa = 'n' or 'n',  op( a ) = a.
!>
!>              transa = 't' or 't',  op( a ) = a**t.
!>
!>              transa = 'c' or 'c',  op( a ) = a**t.
!> \endverbatim
!>
!> \param[in] transb
!> \verbatim
!>          transb is character*1
!>           on entry, transb specifies the form of op( b ) to be used in
!>           the matrix multiplication as follows:
!>
!>              transb = 'n' or 'n',  op( b ) = b.
!>
!>              transb = 't' or 't',  op( b ) = b**t.
!>
!>              transb = 'c' or 'c',  op( b ) = b**t.
!> \endverbatim
!>
!> \param[in] m
!> \verbatim
!>          m is integer
!>           on entry,  m  specifies  the number  of rows  of the  matrix
!>           op( a )  and of the  matrix  c.  m  must  be at least  zero.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is integer
!>           on entry,  n  specifies the number  of columns of the matrix
!>           op( b ) and the number of columns of the matrix c. n must be
!>           at least zero.
!> \endverbatim
!>
!> \param[in] k
!> \verbatim
!>          k is integer
!>           on entry,  k  specifies  the number of columns of the matrix
!>           op( a ) and the number of rows of the matrix op( b ). k must
!>           be at least  zero.
!> \endverbatim
!>
!> \param[in] alpha
!> \verbatim
!>          alpha is double precision.
!>           on entry, alpha specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] a
!> \verbatim
!>          a is double precision array, dimension ( lda, ka ), where ka is
!>           k  when  transa = 'n' or 'n',  and is  m  otherwise.
!>           before entry with  transa = 'n' or 'n',  the leading  m by k
!>           part of the array  a  must contain the matrix  a,  otherwise
!>           the leading  k by m  part of the array  a  must contain  the
!>           matrix a.
!> \endverbatim
!>
!> \param[in] lda
!> \verbatim
!>          lda is integer
!>           on entry, lda specifies the first dimension of a as declared
!>           in the calling (sub) program. when  transa = 'n' or 'n' then
!>           lda must be at least  max( 1, m ), otherwise  lda must be at
!>           least  max( 1, k ).
!> \endverbatim
!>
!> \param[in] b
!> \verbatim
!>          b is double precision array, dimension ( ldb, kb ), where kb is
!>           n  when  transb = 'n' or 'n',  and is  k  otherwise.
!>           before entry with  transb = 'n' or 'n',  the leading  k by n
!>           part of the array  b  must contain the matrix  b,  otherwise
!>           the leading  n by k  part of the array  b  must contain  the
!>           matrix b.
!> \endverbatim
!>
!> \param[in] ldb
!> \verbatim
!>          ldb is integer
!>           on entry, ldb specifies the first dimension of b as declared
!>           in the calling (sub) program. when  transb = 'n' or 'n' then
!>           ldb must be at least  max( 1, k ), otherwise  ldb must be at
!>           least  max( 1, n ).
!> \endverbatim
!>
!> \param[in] beta
!> \verbatim
!>          beta is double precision.
!>           on entry,  beta  specifies the scalar  beta.  when  beta  is
!>           supplied as zero then c need not be set on input.
!> \endverbatim
!>
!> \param[in,out] c
!> \verbatim
!>          c is double precision array, dimension ( ldc, n )
!>           before entry, the leading  m by n  part of the array  c must
!>           contain the matrix  c,  except when  beta  is zero, in which
!>           case c need not be set on entry.
!>           on exit, the array  c  is overwritten by the  m by n  matrix
!>           ( alpha*op( a )*op( b ) + beta*c ).
!> \endverbatim
!>
!> \param[in] ldc
!> \verbatim
!>          ldc is integer
!>           on entry, ldc specifies the first dimension of c as declared
!>           in  the  calling  (sub)  program.   ldc  must  be  at  least
!>           max( 1, m ).
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup gemm
!
!> \par further details:
!  =====================
!>
!> \verbatim
!>
!>  level 3 blas routine.
!>
!>  -- written on 8-february-1989.
!>     jack dongarra, argonne national laboratory.
!>     iain duff, aere harwell.
!>     jeremy du croz, numerical algorithms group ltd.
!>     sven hammarling, numerical algorithms group ltd.
!> \endverbatim
!>
!  =====================================================================
subroutine dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb, beta,c,ldc)
implicit none
!
!  -- reference blas level3 routine --
!  -- reference blas is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
real(dp), intent(in) :: alpha !! on entry, alpha specifies the scalar alpha.
real(dp), intent(in) :: beta !! on entry,  beta  specifies the scalar  beta.  when  beta  is
!! supplied as zero then c need not be set on input.
integer, intent(in) :: k !! on entry,  k  specifies  the number of columns of the matrix op( a )
!! and the number of rows of the matrix op( b ). k must
integer, intent(in) :: lda !! on entry, lda specifies the first dimension of a as declared in the
!! calling (sub) program. when  transa = 'n' or 'n' then
integer, intent(in) :: ldb !! on entry, ldb specifies the first dimension of b as declared in the
!! calling (sub) program. when  transb = 'n' or 'n' then
integer, intent(in) :: ldc !! on entry, ldc specifies the first dimension of c as declared in  the
!! calling  (sub)  program.   ldc  must  be  at  least
integer, intent(in) :: m !! on entry,  m  specifies  the number  of rows  of the  matrix op( a )
!! and of the  matrix  c.  m  must  be at least  zero.
integer, intent(in) :: n !! on entry,  n  specifies the number  of columns of the matrix op( b )
!! and the number of columns of the matrix c. n must be
character, intent(in) :: transa !! on entry, transa specifies the form of op( a ) to be used in
!! the matrix multiplication as follows:
character, intent(in) :: transb !! on entry, transb specifies the form of op( b ) to be used in
!! the matrix multiplication as follows:
!     ..
!     .. array arguments ..
real(dp), intent(in) :: a(lda,*) !! k  when  transa = 'n' or 'n',  and is  m  otherwise. before
!! entry with  transa = 'n' or 'n',  the leading  m by k
real(dp), intent(in) :: b(ldb,*) !! n  when  transb = 'n' or 'n',  and is  k  otherwise. before
!! entry with  transb = 'n' or 'n',  the leading  k by n
real(dp), intent(inout) :: c(ldc,*) !! before entry, the leading  m by n  part of the array  c
!! must contain the matrix  c,  except when  beta  is zero, in
!! which
!     ..
!
!  =====================================================================
!
!     .. external functions ..
!     ..
!     .. external subroutines ..
!     ..
!     .. intrinsic functions ..
intrinsic max
!     ..
!     .. local scalars ..
real(dp) :: temp
integer :: i, info, j, l, nrowa, nrowb
logical :: nota, notb
!     ..
!     .. parameters ..
real(dp) :: one, zero
parameter (one=1.0d+0,zero=0.0d+0)
!     ..
!
!     set  nota  and  notb  as  true if  a  and  b  respectively are not
!     transposed and set  nrowa and nrowb  as the number of rows of  a
!     and  b  respectively.
!
nota = lsame(transa,'n')
notb = lsame(transb,'n')
if (nota) then
nrowa = m
else
nrowa = k
end if
if (notb) then
nrowb = k
else
nrowb = n
end if
!
!     test the input parameters.
!
info = 0
if ((.not.nota) .and. (.not.lsame(transa,'c')) .and. (.not.lsame(transa,'t'))) then
info = 1
else if ((.not.notb) .and. (.not.lsame(transb,'c')) .and. (.not.lsame(transb,'t'))) then
info = 2
else if (m.lt.0) then
info = 3
else if (n.lt.0) then
info = 4
else if (k.lt.0) then
info = 5
else if (lda.lt.max(1,nrowa)) then
info = 8
else if (ldb.lt.max(1,nrowb)) then
info = 10
else if (ldc.lt.max(1,m)) then
info = 13
end if
if (info.ne.0) then
call xerbla('dgemm ',info)
return
end if
!
!     quick return if possible.
!
if ((m.eq.0) .or. (n.eq.0) .or. (((alpha.eq.zero).or. (k.eq.0)).and. (beta.eq.one))) return
!
!     and if  alpha.eq.zero.
!
if (alpha.eq.zero) then
if (beta.eq.zero) then
do 20 j = 1,n
do 10 i = 1,m
c(i,j) = zero
10 continue
20 continue
else
do 40 j = 1,n
do 30 i = 1,m
c(i,j) = beta*c(i,j)
30 continue
40 continue
end if
return
end if
!
!     start the operations.
!
if (notb) then
if (nota) then
!
!           form  c := alpha*a*b + beta*c.
!
do 90 j = 1,n
if (beta.eq.zero) then
do 50 i = 1,m
c(i,j) = zero
50 continue
else if (beta.ne.one) then
do 60 i = 1,m
c(i,j) = beta*c(i,j)
60 continue
end if
do 80 l = 1,k
temp = alpha*b(l,j)
do 70 i = 1,m
c(i,j) = c(i,j) + temp*a(i,l)
70 continue
80 continue
90 continue
else
!
!           form  c := alpha*a**t*b + beta*c
!
do 120 j = 1,n
do 110 i = 1,m
temp = zero
do 100 l = 1,k
temp = temp + a(l,i)*b(l,j)
100 continue
if (beta.eq.zero) then
c(i,j) = alpha*temp
else
c(i,j) = alpha*temp + beta*c(i,j)
end if
110 continue
120 continue
end if
else
if (nota) then
!
!           form  c := alpha*a*b**t + beta*c
!
do 170 j = 1,n
if (beta.eq.zero) then
do 130 i = 1,m
c(i,j) = zero
130 continue
else if (beta.ne.one) then
do 140 i = 1,m
c(i,j) = beta*c(i,j)
140 continue
end if
do 160 l = 1,k
temp = alpha*b(j,l)
do 150 i = 1,m
c(i,j) = c(i,j) + temp*a(i,l)
150 continue
160 continue
170 continue
else
!
!           form  c := alpha*a**t*b**t + beta*c
!
do 200 j = 1,n
do 190 i = 1,m
temp = zero
do 180 l = 1,k
temp = temp + a(l,i)*b(j,l)
180 continue
if (beta.eq.zero) then
c(i,j) = alpha*temp
else
c(i,j) = alpha*temp + beta*c(i,j)
end if
190 continue
200 continue
end if
end if
!
return
!
!     end of dgemm
!
end subroutine dgemm
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dgeqr2.f
!> \brief \b dgeqr2 computes the qr factorization of a general rectangular matrix using an
!unblocked algorithm.
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download dgeqr2 + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeqr2.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeqr2.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeqr2.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       subroutine dgeqr2( m, n, a, lda, tau, work, info )
!
!       .. scalar arguments ..
!       integer            info, lda, m, n
!       ..
!       .. array arguments ..
!       double precision   a( lda, * ), tau( * ), work( * )
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dgeqr2 computes a qr factorization of a real m-by-n matrix a:
!>
!>    a = q * ( r ),
!>            ( 0 )
!>
!> where:
!>
!>    q is a m-by-m orthogonal matrix;
!>    r is an upper-triangular n-by-n matrix;
!>    0 is a (m-n)-by-n zero matrix, if m > n.
!>
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] m
!> \verbatim
!>          m is integer
!>          the number of rows of the matrix a.  m >= 0.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is integer
!>          the number of columns of the matrix a.  n >= 0.
!> \endverbatim
!>
!> \param[in,out] a
!> \verbatim
!>          a is double precision array, dimension (lda,n)
!>          on entry, the m by n matrix a.
!>          on exit, the elements on and above the diagonal of the array
!>          contain the min(m,n) by n upper trapezoidal matrix r (r is
!>          upper triangular if m >= n); the elements below the diagonal,
!>          with the array tau, represent the orthogonal matrix q as a
!>          product of elementary reflectors (see further details).
!> \endverbatim
!>
!> \param[in] lda
!> \verbatim
!>          lda is integer
!>          the leading dimension of the array a.  lda >= max(1,m).
!> \endverbatim
!>
!> \param[out] tau
!> \verbatim
!>          tau is double precision array, dimension (min(m,n))
!>          the scalar factors of the elementary reflectors (see further
!>          details).
!> \endverbatim
!>
!> \param[out] work
!> \verbatim
!>          work is double precision array, dimension (n)
!> \endverbatim
!>
!> \param[out] info
!> \verbatim
!>          info is integer
!>          = 0: successful exit
!>          < 0: if info = -i, the i-th argument had an illegal value
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup geqr2
!
!> \par further details:
!  =====================
!>
!> \verbatim
!>
!>  the matrix q is represented as a product of elementary reflectors
!>
!>     q = h(1) h(2) . . . h(k), where k = min(m,n).
!>
!>  each h(i) has the form
!>
!>     h(i) = i - tau * v * v**t
!>
!>  where tau is a real scalar, and v is a real vector with
!>  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in a(i+1:m,i),
!>  and tau in tau(i).
!> \endverbatim
!>
!  =====================================================================
subroutine dgeqr2( m, n, a, lda, tau, work, info )
implicit none
!
!  -- lapack computational routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
integer, intent(out) :: info !! zero success; negative invalid argument or positive numerical failure.
integer, intent(in) :: lda !! the leading dimension of the array a.  lda >= max(1,m).
integer, intent(in) :: m !! the number of rows of the matrix a.  m >= 0.
integer, intent(in) :: n !! the number of columns of the matrix a.  n >= 0.
!     ..
!     .. array arguments ..
real(dp), intent(inout) :: a( lda, * ) !! on entry, the m by n matrix a. on exit, the elements on
!! and above the diagonal of the array
real(dp), intent(out) :: tau( * ) !! the scalar factors of the elementary reflectors (see further
!! details).
real(dp), intent(out) :: work( * ) !! caller-owned scratch array, sized according to the retained contract.
!     ..
!
!  =====================================================================
!
!     .. parameters ..
real(dp) :: one
parameter          ( one = 1.0d+0 )
!     ..
!     .. local scalars ..
integer :: i, k
!     ..
!     .. external subroutines ..
!     ..
!     .. intrinsic functions ..
intrinsic          max, min
!     ..
!     .. executable statements ..
!
!     test the input arguments
!
info = 0
if( m.lt.0 ) then
info = -1
else if( n.lt.0 ) then
info = -2
else if( lda.lt.max( 1, m ) ) then
info = -4
end if
if( info.ne.0 ) then
call xerbla( 'dgeqr2', -info )
return
end if
!
k = min( m, n )
!
do 10 i = 1, k
!
!        generate elementary reflector h(i) to annihilate a(i+1:m,i)
!
call dlarfg( m-i+1, a( i, i ), a( min( i+1, m ), i ), 1, tau( i ) )
if( i.lt.n ) then
!
!           apply h(i) to a(i:m,i+1:n) from the left
!
call dlarf1f( 'left', m-i+1, n-i, a( i, i ), 1, tau( i ), a( i, i+1 ), lda, work )
end if
10 continue
return
!
!     end of dgeqr2
!
end subroutine dgeqr2
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dgerq2.f
!> \brief \b dgerq2 computes the rq factorization of a general rectangular matrix using an
!unblocked algorithm.
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download dgerq2 + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgerq2.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgerq2.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgerq2.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       subroutine dgerq2( m, n, a, lda, tau, work, info )
!
!       .. scalar arguments ..
!       integer            info, lda, m, n
!       ..
!       .. array arguments ..
!       double precision   a( lda, * ), tau( * ), work( * )
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dgerq2 computes an rq factorization of a real m by n matrix a:
!> a = r * q.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] m
!> \verbatim
!>          m is integer
!>          the number of rows of the matrix a.  m >= 0.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is integer
!>          the number of columns of the matrix a.  n >= 0.
!> \endverbatim
!>
!> \param[in,out] a
!> \verbatim
!>          a is double precision array, dimension (lda,n)
!>          on entry, the m by n matrix a.
!>          on exit, if m <= n, the upper triangle of the subarray
!>          a(1:m,n-m+1:n) contains the m by m upper triangular matrix r;
!>          if m >= n, the elements on and above the (m-n)-th subdiagonal
!>          contain the m by n upper trapezoidal matrix r; the remaining
!>          elements, with the array tau, represent the orthogonal matrix
!>          q as a product of elementary reflectors (see further
!>          details).
!> \endverbatim
!>
!> \param[in] lda
!> \verbatim
!>          lda is integer
!>          the leading dimension of the array a.  lda >= max(1,m).
!> \endverbatim
!>
!> \param[out] tau
!> \verbatim
!>          tau is double precision array, dimension (min(m,n))
!>          the scalar factors of the elementary reflectors (see further
!>          details).
!> \endverbatim
!>
!> \param[out] work
!> \verbatim
!>          work is double precision array, dimension (m)
!> \endverbatim
!>
!> \param[out] info
!> \verbatim
!>          info is integer
!>          = 0: successful exit
!>          < 0: if info = -i, the i-th argument had an illegal value
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup gerq2
!
!> \par further details:
!  =====================
!>
!> \verbatim
!>
!>  the matrix q is represented as a product of elementary reflectors
!>
!>     q = h(1) h(2) . . . h(k), where k = min(m,n).
!>
!>  each h(i) has the form
!>
!>     h(i) = i - tau * v * v**t
!>
!>  where tau is a real scalar, and v is a real vector with
!>  v(n-k+i+1:n) = 0 and v(n-k+i) = 1; v(1:n-k+i-1) is stored on exit in
!>  a(m-k+i,1:n-k+i-1), and tau in tau(i).
!> \endverbatim
!>
!  =====================================================================
subroutine dgerq2( m, n, a, lda, tau, work, info )
implicit none
!
!  -- lapack computational routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
integer, intent(out) :: info !! zero success; negative invalid argument or positive numerical failure.
integer, intent(in) :: lda !! the leading dimension of the array a.  lda >= max(1,m).
integer, intent(in) :: m !! the number of rows of the matrix a.  m >= 0.
integer, intent(in) :: n !! the number of columns of the matrix a.  n >= 0.
!     ..
!     .. array arguments ..
real(dp), intent(inout) :: a( lda, * ) !! on entry, the m by n matrix a. on exit, if m <= n, the
!! upper triangle of the subarray a(1:m,n-m+1:n) contains
!! the m by m upper triangular matrix r;
real(dp), intent(out) :: tau( * ) !! the scalar factors of the elementary reflectors (see further
!! details).
real(dp), intent(out) :: work( * ) !! caller-owned scratch array, sized according to the retained contract.
!     ..
!
!  =====================================================================
!
!     .. parameters ..
real(dp) :: one
parameter          ( one = 1.0d+0 )
!     ..
!     .. local scalars ..
integer :: i, k
!     ..
!     .. external subroutines ..
!     ..
!     .. intrinsic functions ..
intrinsic          max, min
!     ..
!     .. executable statements ..
!
!     test the input arguments
!
info = 0
if( m.lt.0 ) then
info = -1
else if( n.lt.0 ) then
info = -2
else if( lda.lt.max( 1, m ) ) then
info = -4
end if
if( info.ne.0 ) then
call xerbla( 'dgerq2', -info )
return
end if
!
k = min( m, n )
!
do 10 i = k, 1, -1
!
!        generate elementary reflector h(i) to annihilate
!        a(m-k+i,1:n-k+i-1)
!
call dlarfg( n-k+i, a( m-k+i, n-k+i ), a( m-k+i, 1 ), lda, tau( i ) )
!
!        apply h(i) to a(1:m-k+i-1,1:n-k+i) from the right
!
call dlarf1l( 'right', m-k+i-1, n-k+i, a( m-k+i, 1 ), lda, tau( i ), a, lda, work )
10 continue
return
!
!     end of dgerq2
!
end subroutine dgerq2
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlacpy.f
!> \brief \b dlacpy copies all or part of one two-dimensional array to another.
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download dlacpy + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlacpy.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlacpy.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlacpy.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       subroutine dlacpy( uplo, m, n, a, lda, b, ldb )
!
!       .. scalar arguments ..
!       character          uplo
!       integer            lda, ldb, m, n
!       ..
!       .. array arguments ..
!       double precision   a( lda, * ), b( ldb, * )
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dlacpy copies all or part of a two-dimensional matrix a to another
!> matrix b.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] uplo
!> \verbatim
!>          uplo is character*1
!>          specifies the part of the matrix a to be copied to b.
!>          = 'u':      upper triangular part
!>          = 'l':      lower triangular part
!>          otherwise:  all of the matrix a
!> \endverbatim
!>
!> \param[in] m
!> \verbatim
!>          m is integer
!>          the number of rows of the matrix a.  m >= 0.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is integer
!>          the number of columns of the matrix a.  n >= 0.
!> \endverbatim
!>
!> \param[in] a
!> \verbatim
!>          a is double precision array, dimension (lda,n)
!>          the m by n matrix a.  if uplo = 'u', only the upper triangle
!>          or trapezoid is accessed; if uplo = 'l', only the lower
!>          triangle or trapezoid is accessed.
!> \endverbatim
!>
!> \param[in] lda
!> \verbatim
!>          lda is integer
!>          the leading dimension of the array a.  lda >= max(1,m).
!> \endverbatim
!>
!> \param[out] b
!> \verbatim
!>          b is double precision array, dimension (ldb,n)
!>          on exit, b = a in the locations specified by uplo.
!> \endverbatim
!>
!> \param[in] ldb
!> \verbatim
!>          ldb is integer
!>          the leading dimension of the array b.  ldb >= max(1,m).
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup lacpy
!
!  =====================================================================
subroutine dlacpy( uplo, m, n, a, lda, b, ldb )
implicit none
!
!  -- lapack auxiliary routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
character, intent(in) :: uplo !! specifies the part of the matrix a to be copied to b.
integer, intent(in) :: lda !! the leading dimension of the array a.  lda >= max(1,m).
integer, intent(in) :: ldb !! the leading dimension of the array b.  ldb >= max(1,m).
integer, intent(in) :: m !! the number of rows of the matrix a.  m >= 0.
integer, intent(in) :: n !! the number of columns of the matrix a.  n >= 0.
!     ..
!     .. array arguments ..
real(dp), intent(in) :: a( lda, * ) !! the m by n matrix a.  if uplo = 'u', only the upper
!! triangle or trapezoid is accessed; if uplo = 'l', only the
!! lower
real(dp), intent(out) :: b( ldb, * ) !! on exit, b = a in the locations specified by uplo.
!     ..
!
!  =====================================================================
!
!     .. local scalars ..
integer :: i, j
!     ..
!     .. external functions ..
!     ..
!     .. intrinsic functions ..
intrinsic          min
!     ..
!     .. executable statements ..
!
if( lsame( uplo, 'u' ) ) then
do 20 j = 1, n
do 10 i = 1, min( j, m )
b( i, j ) = a( i, j )
10 continue
20 continue
else if( lsame( uplo, 'l' ) ) then
do 40 j = 1, n
do 30 i = j, m
b( i, j ) = a( i, j )
30 continue
40 continue
else
do 60 j = 1, n
do 50 i = 1, m
b( i, j ) = a( i, j )
50 continue
60 continue
end if
return
!
!     end of dlacpy
!
end subroutine dlacpy
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlagv2.f
!> \brief \b dlagv2 computes the generalized schur factorization of a real 2-by-2 matrix pencil
!(a,b) where b is upper triangular.
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download dlagv2 + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlagv2.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlagv2.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlagv2.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       subroutine dlagv2( a, lda, b, ldb, alphar, alphai, beta, csl, snl,
!                          csr, snr )
!
!       .. scalar arguments ..
!       integer            lda, ldb
!       double precision   csl, csr, snl, snr
!       ..
!       .. array arguments ..
!       double precision   a( lda, * ), alphai( 2 ), alphar( 2 ),
!      $                   b( ldb, * ), beta( 2 )
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dlagv2 computes the generalized schur factorization of a real 2-by-2
!> matrix pencil (a,b) where b is upper triangular. this routine
!> computes orthogonal (rotation) matrices given by csl, snl and csr,
!> snr such that
!>
!> 1) if the pencil (a,b) has two real eigenvalues (include 0/0 or 1/0
!>    types), then
!>
!>    [ a11 a12 ] := [  csl  snl ] [ a11 a12 ] [  csr -snr ]
!>    [  0  a22 ]    [ -snl  csl ] [ a21 a22 ] [  snr  csr ]
!>
!>    [ b11 b12 ] := [  csl  snl ] [ b11 b12 ] [  csr -snr ]
!>    [  0  b22 ]    [ -snl  csl ] [  0  b22 ] [  snr  csr ],
!>
!> 2) if the pencil (a,b) has a pair of complex conjugate eigenvalues,
!>    then
!>
!>    [ a11 a12 ] := [  csl  snl ] [ a11 a12 ] [  csr -snr ]
!>    [ a21 a22 ]    [ -snl  csl ] [ a21 a22 ] [  snr  csr ]
!>
!>    [ b11  0  ] := [  csl  snl ] [ b11 b12 ] [  csr -snr ]
!>    [  0  b22 ]    [ -snl  csl ] [  0  b22 ] [  snr  csr ]
!>
!>    where b11 >= b22 > 0.
!>
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in,out] a
!> \verbatim
!>          a is double precision array, dimension (lda, 2)
!>          on entry, the 2 x 2 matrix a.
!>          on exit, a is overwritten by the ``a-part'' of the
!>          generalized schur form.
!> \endverbatim
!>
!> \param[in] lda
!> \verbatim
!>          lda is integer
!>          the leading dimension of the array a.  lda >= 2.
!> \endverbatim
!>
!> \param[in,out] b
!> \verbatim
!>          b is double precision array, dimension (ldb, 2)
!>          on entry, the upper triangular 2 x 2 matrix b.
!>          on exit, b is overwritten by the ``b-part'' of the
!>          generalized schur form.
!> \endverbatim
!>
!> \param[in] ldb
!> \verbatim
!>          ldb is integer
!>          the leading dimension of the array b.  ldb >= 2.
!> \endverbatim
!>
!> \param[out] alphar
!> \verbatim
!>          alphar is double precision array, dimension (2)
!> \endverbatim
!>
!> \param[out] alphai
!> \verbatim
!>          alphai is double precision array, dimension (2)
!> \endverbatim
!>
!> \param[out] beta
!> \verbatim
!>          beta is double precision array, dimension (2)
!>          (alphar(k)+i*alphai(k))/beta(k) are the eigenvalues of the
!>          pencil (a,b), k=1,2, i = sqrt(-1).  note that beta(k) may
!>          be zero.
!> \endverbatim
!>
!> \param[out] csl
!> \verbatim
!>          csl is double precision
!>          the cosine of the left rotation matrix.
!> \endverbatim
!>
!> \param[out] snl
!> \verbatim
!>          snl is double precision
!>          the sine of the left rotation matrix.
!> \endverbatim
!>
!> \param[out] csr
!> \verbatim
!>          csr is double precision
!>          the cosine of the right rotation matrix.
!> \endverbatim
!>
!> \param[out] snr
!> \verbatim
!>          snr is double precision
!>          the sine of the right rotation matrix.
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup lagv2
!
!> \par contributors:
!  ==================
!>
!>     mark fahey, department of mathematics, univ. of kentucky, usa
!
!  =====================================================================
subroutine dlagv2( a, lda, b, ldb, alphar, alphai, beta, csl, snl, csr, snr )
implicit none
!
!  -- lapack auxiliary routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
integer, intent(in) :: lda !! the leading dimension of the array a.  lda >= 2.
integer, intent(in) :: ldb !! the leading dimension of the array b.  ldb >= 2.
real(dp), intent(out) :: csl !! the cosine of the left rotation matrix.
real(dp), intent(out) :: csr !! the cosine of the right rotation matrix.
real(dp), intent(out) :: snl !! the sine of the left rotation matrix.
real(dp), intent(out) :: snr !! the sine of the right rotation matrix.
!     ..
!     .. array arguments ..
real(dp), intent(inout) :: a( lda, * ) !! on entry, the 2 x 2 matrix a. on exit, a is overwritten
!! by the ``a-part'' of the generalized schur form.
real(dp), intent(out) :: alphai( 2 ) !! imaginary parts of the two homogeneous eigenvalue numerators.
real(dp), intent(out) :: alphar( 2 ) !! real parts of the two homogeneous eigenvalue numerators.
real(dp), intent(inout) :: b( ldb, * ) !! on entry, the upper triangular 2 x 2 matrix b. on exit,
!! b is overwritten by the ``b-part'' of the
real(dp), intent(out) :: beta( 2 ) !! (alphar(k)+i*alphai(k))/beta(k) are the eigenvalues of the
!! pencil (a,b), k=1,2, i = sqrt(-1).  note that beta(k) may
!     ..
!
!  =====================================================================
!
!     .. parameters ..
real(dp) :: zero, one
parameter          ( zero = 0.0d+0, one = 1.0d+0 )
!     ..
!     .. local scalars ..
real(dp) :: anorm, ascale, bnorm, bscale, h1, h2, h3, qq, r, rr, safmin, scale1, scale2, t, ulp, wi, wr1, wr2
!     ..
!     .. external subroutines ..
!     ..
!     .. external functions ..
!     ..
!     .. intrinsic functions ..
intrinsic          abs, max
!     ..
!     .. executable statements ..
!
safmin = dlamch( 's' )
ulp = dlamch( 'p' )
!
!     scale a
!
anorm = max( abs( a( 1, 1 ) )+abs( a( 2, 1 ) ), abs( a( 1, 2 ) )+abs( a( 2, 2 ) ), safmin )
ascale = one / anorm
a( 1, 1 ) = ascale*a( 1, 1 )
a( 1, 2 ) = ascale*a( 1, 2 )
a( 2, 1 ) = ascale*a( 2, 1 )
a( 2, 2 ) = ascale*a( 2, 2 )
!
!     scale b
!
bnorm = max( abs( b( 1, 1 ) ), abs( b( 1, 2 ) )+abs( b( 2, 2 ) ), safmin )
bscale = one / bnorm
b( 1, 1 ) = bscale*b( 1, 1 )
b( 1, 2 ) = bscale*b( 1, 2 )
b( 2, 2 ) = bscale*b( 2, 2 )
!
!     check if a can be deflated
!
if( abs( a( 2, 1 ) ).le.ulp ) then
csl = one
snl = zero
csr = one
snr = zero
a( 2, 1 ) = zero
b( 2, 1 ) = zero
wi = zero
!
!     check if b is singular
!
else if( abs( b( 1, 1 ) ).le.ulp ) then
call dlartg( a( 1, 1 ), a( 2, 1 ), csl, snl, r )
csr = one
snr = zero
call drot( 2, a( 1, 1 ), lda, a( 2, 1 ), lda, csl, snl )
call drot( 2, b( 1, 1 ), ldb, b( 2, 1 ), ldb, csl, snl )
a( 2, 1 ) = zero
b( 1, 1 ) = zero
b( 2, 1 ) = zero
wi = zero
!
else if( abs( b( 2, 2 ) ).le.ulp ) then
call dlartg( a( 2, 2 ), a( 2, 1 ), csr, snr, t )
snr = -snr
call drot( 2, a( 1, 1 ), 1, a( 1, 2 ), 1, csr, snr )
call drot( 2, b( 1, 1 ), 1, b( 1, 2 ), 1, csr, snr )
csl = one
snl = zero
a( 2, 1 ) = zero
b( 2, 1 ) = zero
b( 2, 2 ) = zero
wi = zero
!
else
!
!        b is nonsingular, first compute the eigenvalues of (a,b)
!
call dlag2( a, lda, b, ldb, safmin, scale1, scale2, wr1, wr2, wi )
!
if( wi.eq.zero ) then
!
!           two real eigenvalues, compute s*a-w*b
!
h1 = scale1*a( 1, 1 ) - wr1*b( 1, 1 )
h2 = scale1*a( 1, 2 ) - wr1*b( 1, 2 )
h3 = scale1*a( 2, 2 ) - wr1*b( 2, 2 )
!
rr = dlapy2( h1, h2 )
qq = dlapy2( scale1*a( 2, 1 ), h3 )
!
if( rr.gt.qq ) then
!
!              find right rotation matrix to zero 1,1 element of
!              (sa - wb)
!
call dlartg( h2, h1, csr, snr, t )
!
else
!
!              find right rotation matrix to zero 2,1 element of
!              (sa - wb)
!
call dlartg( h3, scale1*a( 2, 1 ), csr, snr, t )
!
end if
!
snr = -snr
call drot( 2, a( 1, 1 ), 1, a( 1, 2 ), 1, csr, snr )
call drot( 2, b( 1, 1 ), 1, b( 1, 2 ), 1, csr, snr )
!
!           compute inf norms of a and b
!
h1 = max( abs( a( 1, 1 ) )+abs( a( 1, 2 ) ), abs( a( 2, 1 ) )+abs( a( 2, 2 ) ) )
h2 = max( abs( b( 1, 1 ) )+abs( b( 1, 2 ) ), abs( b( 2, 1 ) )+abs( b( 2, 2 ) ) )
!
if( ( scale1*h1 ).ge.abs( wr1 )*h2 ) then
!
!              find left rotation matrix q to zero out b(2,1)
!
call dlartg( b( 1, 1 ), b( 2, 1 ), csl, snl, r )
!
else
!
!              find left rotation matrix q to zero out a(2,1)
!
call dlartg( a( 1, 1 ), a( 2, 1 ), csl, snl, r )
!
end if
!
call drot( 2, a( 1, 1 ), lda, a( 2, 1 ), lda, csl, snl )
call drot( 2, b( 1, 1 ), ldb, b( 2, 1 ), ldb, csl, snl )
!
a( 2, 1 ) = zero
b( 2, 1 ) = zero
!
else
!
!           a pair of complex conjugate eigenvalues
!           first compute the svd of the matrix b
!
call dlasv2( b( 1, 1 ), b( 1, 2 ), b( 2, 2 ), r, t, snr, csr, snl, csl )
!
!           form (a,b) := q(a,b)z**t where q is left rotation matrix and
!           z is right rotation matrix computed from dlasv2
!
call drot( 2, a( 1, 1 ), lda, a( 2, 1 ), lda, csl, snl )
call drot( 2, b( 1, 1 ), ldb, b( 2, 1 ), ldb, csl, snl )
call drot( 2, a( 1, 1 ), 1, a( 1, 2 ), 1, csr, snr )
call drot( 2, b( 1, 1 ), 1, b( 1, 2 ), 1, csr, snr )
!
b( 2, 1 ) = zero
b( 1, 2 ) = zero
!
end if
!
end if
!
!     unscaling
!
a( 1, 1 ) = anorm*a( 1, 1 )
a( 2, 1 ) = anorm*a( 2, 1 )
a( 1, 2 ) = anorm*a( 1, 2 )
a( 2, 2 ) = anorm*a( 2, 2 )
b( 1, 1 ) = bnorm*b( 1, 1 )
b( 2, 1 ) = bnorm*b( 2, 1 )
b( 1, 2 ) = bnorm*b( 1, 2 )
b( 2, 2 ) = bnorm*b( 2, 2 )
!
if( wi.eq.zero ) then
alphar( 1 ) = a( 1, 1 )
alphar( 2 ) = a( 2, 2 )
alphai( 1 ) = zero
alphai( 2 ) = zero
beta( 1 ) = b( 1, 1 )
beta( 2 ) = b( 2, 2 )
else
alphar( 1 ) = anorm*wr1 / scale1 / bnorm
alphai( 1 ) = anorm*wi / scale1 / bnorm
alphar( 2 ) = alphar( 1 )
alphai( 2 ) = -alphai( 1 )
beta( 1 ) = one
beta( 2 ) = one
end if
!
return
!
!     end of dlagv2
!
end subroutine dlagv2
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dorg2r.f
!> \brief \b dorg2r generates all or part of the orthogonal matrix q from a qr factorization
!determined by sgeqrf (unblocked algorithm).
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download dorg2r + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorg2r.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorg2r.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorg2r.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       subroutine dorg2r( m, n, k, a, lda, tau, work, info )
!
!       .. scalar arguments ..
!       integer            info, k, lda, m, n
!       ..
!       .. array arguments ..
!       double precision   a( lda, * ), tau( * ), work( * )
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dorg2r generates an m by n real matrix q with orthonormal columns,
!> which is defined as the first n columns of a product of k elementary
!> reflectors of order m
!>
!>       q  =  h(1) h(2) . . . h(k)
!>
!> as returned by dgeqrf.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] m
!> \verbatim
!>          m is integer
!>          the number of rows of the matrix q. m >= 0.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is integer
!>          the number of columns of the matrix q. m >= n >= 0.
!> \endverbatim
!>
!> \param[in] k
!> \verbatim
!>          k is integer
!>          the number of elementary reflectors whose product defines the
!>          matrix q. n >= k >= 0.
!> \endverbatim
!>
!> \param[in,out] a
!> \verbatim
!>          a is double precision array, dimension (lda,n)
!>          on entry, the i-th column must contain the vector which
!>          defines the elementary reflector h(i), for i = 1,2,...,k, as
!>          returned by dgeqrf in the first k columns of its array
!>          argument a.
!>          on exit, the m-by-n matrix q.
!> \endverbatim
!>
!> \param[in] lda
!> \verbatim
!>          lda is integer
!>          the first dimension of the array a. lda >= max(1,m).
!> \endverbatim
!>
!> \param[in] tau
!> \verbatim
!>          tau is double precision array, dimension (k)
!>          tau(i) must contain the scalar factor of the elementary
!>          reflector h(i), as returned by dgeqrf.
!> \endverbatim
!>
!> \param[out] work
!> \verbatim
!>          work is double precision array, dimension (n)
!> \endverbatim
!>
!> \param[out] info
!> \verbatim
!>          info is integer
!>          = 0: successful exit
!>          < 0: if info = -i, the i-th argument has an illegal value
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup ung2r
!
!  =====================================================================
subroutine dorg2r( m, n, k, a, lda, tau, work, info )
implicit none
!
!  -- lapack computational routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
integer, intent(out) :: info !! zero success; negative invalid argument or positive numerical failure.
integer, intent(in) :: k !! the number of elementary reflectors whose product defines the matrix
!! q. n >= k >= 0.
integer, intent(in) :: lda !! the first dimension of the array a. lda >= max(1,m).
integer, intent(in) :: m !! the number of rows of the matrix q. m >= 0.
integer, intent(in) :: n !! the number of columns of the matrix q. m >= n >= 0.
!     ..
!     .. array arguments ..
real(dp), intent(inout) :: a( lda, * ) !! on entry, the i-th column must contain the vector which
!! defines the elementary reflector h(i), for i =
!! 1,2,...,k, as
real(dp), intent(in) :: tau( * ) !! tau(i) must contain the scalar factor of the elementary
!! reflector h(i), as returned by dgeqrf.
real(dp), intent(out) :: work( * ) !! caller-owned scratch array, sized according to the retained contract.
!     ..
!
!  =====================================================================
!
!     .. parameters ..
real(dp) :: one, zero
parameter          ( one = 1.0d+0, zero = 0.0d+0 )
!     ..
!     .. local scalars ..
integer :: i, j, l
!     ..
!     .. external subroutines ..
!     ..
!     .. intrinsic functions ..
intrinsic          max
!     ..
!     .. executable statements ..
!
!     test the input arguments
!
info = 0
if( m.lt.0 ) then
info = -1
else if( n.lt.0 .or. n.gt.m ) then
info = -2
else if( k.lt.0 .or. k.gt.n ) then
info = -3
else if( lda.lt.max( 1, m ) ) then
info = -5
end if
if( info.ne.0 ) then
call xerbla( 'dorg2r', -info )
return
end if
!
!     quick return if possible
!
if( n.le.0 ) return
!
!     initialise columns k+1:n to columns of the unit matrix
!
do 20 j = k + 1, n
do 10 l = 1, m
a( l, j ) = zero
10 continue
a( j, j ) = one
20 continue
!
do 40 i = k, 1, -1
!
!        apply h(i) to a(i:m,i:n) from the left
!
if( i.lt.n ) then
call dlarf1f( 'left', m-i+1, n-i, a( i, i ), 1, tau( i ), a( i, i+1 ), lda, work )
end if
if( i.lt.m ) call dscal( m-i, -tau( i ), a( i+1, i ), 1 )
a( i, i ) = one - tau( i )
!
!        set a(1:i-1,i) to zero
!
do 30 l = 1, i - 1
a( l, i ) = zero
30 continue
40 continue
return
!
!     end of dorg2r
!
end subroutine dorg2r
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dorgr2.f
!> \brief \b dorgr2 generates all or part of the orthogonal matrix q from an rq factorization
!determined by sgerqf (unblocked algorithm).
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download dorgr2 + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorgr2.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorgr2.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorgr2.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       subroutine dorgr2( m, n, k, a, lda, tau, work, info )
!
!       .. scalar arguments ..
!       integer            info, k, lda, m, n
!       ..
!       .. array arguments ..
!       double precision   a( lda, * ), tau( * ), work( * )
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dorgr2 generates an m by n real matrix q with orthonormal rows,
!> which is defined as the last m rows of a product of k elementary
!> reflectors of order n
!>
!>       q  =  h(1) h(2) . . . h(k)
!>
!> as returned by dgerqf.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] m
!> \verbatim
!>          m is integer
!>          the number of rows of the matrix q. m >= 0.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is integer
!>          the number of columns of the matrix q. n >= m.
!> \endverbatim
!>
!> \param[in] k
!> \verbatim
!>          k is integer
!>          the number of elementary reflectors whose product defines the
!>          matrix q. m >= k >= 0.
!> \endverbatim
!>
!> \param[in,out] a
!> \verbatim
!>          a is double precision array, dimension (lda,n)
!>          on entry, the (m-k+i)-th row must contain the vector which
!>          defines the elementary reflector h(i), for i = 1,2,...,k, as
!>          returned by dgerqf in the last k rows of its array argument
!>          a.
!>          on exit, the m by n matrix q.
!> \endverbatim
!>
!> \param[in] lda
!> \verbatim
!>          lda is integer
!>          the first dimension of the array a. lda >= max(1,m).
!> \endverbatim
!>
!> \param[in] tau
!> \verbatim
!>          tau is double precision array, dimension (k)
!>          tau(i) must contain the scalar factor of the elementary
!>          reflector h(i), as returned by dgerqf.
!> \endverbatim
!>
!> \param[out] work
!> \verbatim
!>          work is double precision array, dimension (m)
!> \endverbatim
!>
!> \param[out] info
!> \verbatim
!>          info is integer
!>          = 0: successful exit
!>          < 0: if info = -i, the i-th argument has an illegal value
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup ungr2
!
!  =====================================================================
subroutine dorgr2( m, n, k, a, lda, tau, work, info )
implicit none
!
!  -- lapack computational routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
integer, intent(out) :: info !! zero success; negative invalid argument or positive numerical failure.
integer, intent(in) :: k !! the number of elementary reflectors whose product defines the matrix
!! q. m >= k >= 0.
integer, intent(in) :: lda !! the first dimension of the array a. lda >= max(1,m).
integer, intent(in) :: m !! the number of rows of the matrix q. m >= 0.
integer, intent(in) :: n !! the number of columns of the matrix q. n >= m.
!     ..
!     .. array arguments ..
real(dp), intent(inout) :: a( lda, * ) !! on entry, the (m-k+i)-th row must contain the vector
!! which defines the elementary reflector h(i), for i =
!! 1,2,...,k, as
real(dp), intent(in) :: tau( * ) !! tau(i) must contain the scalar factor of the elementary
!! reflector h(i), as returned by dgerqf.
real(dp), intent(out) :: work( * ) !! caller-owned scratch array, sized according to the retained contract.
!     ..
!
!  =====================================================================
!
!     .. parameters ..
real(dp) :: one, zero
parameter          ( one = 1.0d+0, zero = 0.0d+0 )
!     ..
!     .. local scalars ..
integer :: i, ii, j, l
!     ..
!     .. external subroutines ..
!     ..
!     .. intrinsic functions ..
intrinsic          max
!     ..
!     .. executable statements ..
!
!     test the input arguments
!
info = 0
if( m.lt.0 ) then
info = -1
else if( n.lt.m ) then
info = -2
else if( k.lt.0 .or. k.gt.m ) then
info = -3
else if( lda.lt.max( 1, m ) ) then
info = -5
end if
if( info.ne.0 ) then
call xerbla( 'dorgr2', -info )
return
end if
!
!     quick return if possible
!
if( m.le.0 ) return
!
if( k.lt.m ) then
!
!        initialise rows 1:m-k to rows of the unit matrix
!
do 20 j = 1, n
do 10 l = 1, m - k
a( l, j ) = zero
10 continue
if( j.gt.n-m .and. j.le.n-k ) a( m-n+j, j ) = one
20 continue
end if
!
do 40 i = 1, k
ii = m - k + i
!
!        apply h(i) to a(1:m-k+i,1:n-k+i) from the right
!
!a( ii, n-m+ii ) = one
call dlarf1l( 'right', ii-1, n-m+ii, a( ii, 1 ), lda, tau( i ), a, lda, work )
call dscal( n-m+ii-1, -tau( i ), a( ii, 1 ), lda )
a( ii, n-m+ii ) = one - tau( i )
!
!        set a(m-k+i,n-k+i+1:n) to zero
!
do 30 l = n - m + ii + 1, n
a( ii, l ) = zero
30 continue
40 continue
return
!
!     end of dorgr2
!
end subroutine dorgr2
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dorm2r.f
!> \brief \b dorm2r multiplies a general matrix by the orthogonal matrix from a qr factorization
!determined by sgeqrf (unblocked algorithm).
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download dorm2r + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorm2r.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorm2r.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorm2r.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       subroutine dorm2r( side, trans, m, n, k, a, lda, tau, c, ldc,
!                          work, info )
!
!       .. scalar arguments ..
!       character          side, trans
!       integer            info, k, lda, ldc, m, n
!       ..
!       .. array arguments ..
!       double precision   a( lda, * ), c( ldc, * ), tau( * ), work( * )
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dorm2r overwrites the general real m by n matrix c with
!>
!>       q * c  if side = 'l' and trans = 'n', or
!>
!>       q**t* c  if side = 'l' and trans = 't', or
!>
!>       c * q  if side = 'r' and trans = 'n', or
!>
!>       c * q**t if side = 'r' and trans = 't',
!>
!> where q is a real orthogonal matrix defined as the product of k
!> elementary reflectors
!>
!>       q = h(1) h(2) . . . h(k)
!>
!> as returned by dgeqrf. q is of order m if side = 'l' and of order n
!> if side = 'r'.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] side
!> \verbatim
!>          side is character*1
!>          = 'l': apply q or q**t from the left
!>          = 'r': apply q or q**t from the right
!> \endverbatim
!>
!> \param[in] trans
!> \verbatim
!>          trans is character*1
!>          = 'n': apply q  (no transpose)
!>          = 't': apply q**t (transpose)
!> \endverbatim
!>
!> \param[in] m
!> \verbatim
!>          m is integer
!>          the number of rows of the matrix c. m >= 0.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is integer
!>          the number of columns of the matrix c. n >= 0.
!> \endverbatim
!>
!> \param[in] k
!> \verbatim
!>          k is integer
!>          the number of elementary reflectors whose product defines
!>          the matrix q.
!>          if side = 'l', m >= k >= 0;
!>          if side = 'r', n >= k >= 0.
!> \endverbatim
!>
!> \param[in] a
!> \verbatim
!>          a is double precision array, dimension (lda,k)
!>          the i-th column must contain the vector which defines the
!>          elementary reflector h(i), for i = 1,2,...,k, as returned by
!>          dgeqrf in the first k columns of its array argument a.
!> \endverbatim
!>
!> \param[in] lda
!> \verbatim
!>          lda is integer
!>          the leading dimension of the array a.
!>          if side = 'l', lda >= max(1,m);
!>          if side = 'r', lda >= max(1,n).
!> \endverbatim
!>
!> \param[in] tau
!> \verbatim
!>          tau is double precision array, dimension (k)
!>          tau(i) must contain the scalar factor of the elementary
!>          reflector h(i), as returned by dgeqrf.
!> \endverbatim
!>
!> \param[in,out] c
!> \verbatim
!>          c is double precision array, dimension (ldc,n)
!>          on entry, the m by n matrix c.
!>          on exit, c is overwritten by q*c or q**t*c or c*q**t or c*q.
!> \endverbatim
!>
!> \param[in] ldc
!> \verbatim
!>          ldc is integer
!>          the leading dimension of the array c. ldc >= max(1,m).
!> \endverbatim
!>
!> \param[out] work
!> \verbatim
!>          work is double precision array, dimension
!>                                   (n) if side = 'l',
!>                                   (m) if side = 'r'
!> \endverbatim
!>
!> \param[out] info
!> \verbatim
!>          info is integer
!>          = 0: successful exit
!>          < 0: if info = -i, the i-th argument had an illegal value
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup unm2r
!
!  =====================================================================
subroutine dorm2r( side, trans, m, n, k, a, lda, tau, c, ldc, work, info )
implicit none
!
!  -- lapack computational routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
character, intent(in) :: side !! apply the transformation from the left (l) or right (r).
character, intent(in) :: trans !! apply the transformation directly (n) or transposed (t).
integer, intent(out) :: info !! zero success; negative invalid argument or positive numerical failure.
integer, intent(in) :: k !! the number of elementary reflectors whose product defines the matrix
!! q. if side = 'l', m >= k >= 0;
integer, intent(in) :: lda !! the leading dimension of the array a. if side = 'l', lda >=
!! max(1,m); if side = 'r', lda >= max(1,n).
integer, intent(in) :: ldc !! the leading dimension of the array c. ldc >= max(1,m).
integer, intent(in) :: m !! the number of rows of the matrix c. m >= 0.
integer, intent(in) :: n !! the number of columns of the matrix c. n >= 0.
!     ..
!     .. array arguments ..
real(dp), intent(in) :: a( lda, * ) !! the i-th column must contain the vector which defines the
!! elementary reflector h(i), for i = 1,2,...,k, as returned
!! by
real(dp), intent(inout) :: c( ldc, * ) !! on entry, the m by n matrix c. on exit, c is overwritten
!! by q*c or q**t*c or c*q**t or c*q.
real(dp), intent(in) :: tau( * ) !! tau(i) must contain the scalar factor of the elementary
!! reflector h(i), as returned by dgeqrf.
real(dp), intent(out) :: work( * ) !! (n) if side = 'l', (m) if side = 'r'
!     ..
!
!  =====================================================================
!
!     .. parameters ..
real(dp) :: one
parameter          ( one = 1.0d+0 )
!     ..
!     .. local scalars ..
logical :: left, notran
integer :: i, i1, i2, i3, ic, jc, mi, ni, nq
!     ..
!     .. external functions ..
!     ..
!     .. external subroutines ..
!     ..
!     .. intrinsic functions ..
intrinsic          max
!     ..
!     .. executable statements ..
!
!     test the input arguments
!
info = 0
left = lsame( side, 'l' )
notran = lsame( trans, 'n' )
!
!     nq is the order of q
!
if( left ) then
nq = m
else
nq = n
end if
if( .not.left .and. .not.lsame( side, 'r' ) ) then
info = -1
else if( .not.notran .and. .not.lsame( trans, 't' ) ) then
info = -2
else if( m.lt.0 ) then
info = -3
else if( n.lt.0 ) then
info = -4
else if( k.lt.0 .or. k.gt.nq ) then
info = -5
else if( lda.lt.max( 1, nq ) ) then
info = -7
else if( ldc.lt.max( 1, m ) ) then
info = -10
end if
if( info.ne.0 ) then
call xerbla( 'dorm2r', -info )
return
end if
!
!     quick return if possible
!
if( m.eq.0 .or. n.eq.0 .or. k.eq.0 ) return
!
if( ( left .and. .not.notran ) .or. ( .not.left .and. notran ) ) then
i1 = 1
i2 = k
i3 = 1
else
i1 = k
i2 = 1
i3 = -1
end if
!
if( left ) then
ni = n
jc = 1
else
mi = m
ic = 1
end if
!
do 10 i = i1, i2, i3
if( left ) then
!
!           h(i) is applied to c(i:m,1:n)
!
mi = m - i + 1
ic = i
else
!
!           h(i) is applied to c(1:m,i:n)
!
ni = n - i + 1
jc = i
end if
!
!        apply h(i)
!
call dlarf1f( side, mi, ni, a( i, i ), 1, tau( i ), c( ic, jc ), ldc, work )
10 continue
return
!
!     end of dorm2r
!
end subroutine dorm2r
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dormr2.f
!> \brief \b dormr2 multiplies a general matrix by the orthogonal matrix from a rq factorization
!determined by sgerqf (unblocked algorithm).
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download dormr2 + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dormr2.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dormr2.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dormr2.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       subroutine dormr2( side, trans, m, n, k, a, lda, tau, c, ldc,
!                          work, info )
!
!       .. scalar arguments ..
!       character          side, trans
!       integer            info, k, lda, ldc, m, n
!       ..
!       .. array arguments ..
!       double precision   a( lda, * ), c( ldc, * ), tau( * ), work( * )
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dormr2 overwrites the general real m by n matrix c with
!>
!>       q * c  if side = 'l' and trans = 'n', or
!>
!>       q**t* c  if side = 'l' and trans = 't', or
!>
!>       c * q  if side = 'r' and trans = 'n', or
!>
!>       c * q**t if side = 'r' and trans = 't',
!>
!> where q is a real orthogonal matrix defined as the product of k
!> elementary reflectors
!>
!>       q = h(1) h(2) . . . h(k)
!>
!> as returned by dgerqf. q is of order m if side = 'l' and of order n
!> if side = 'r'.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] side
!> \verbatim
!>          side is character*1
!>          = 'l': apply q or q**t from the left
!>          = 'r': apply q or q**t from the right
!> \endverbatim
!>
!> \param[in] trans
!> \verbatim
!>          trans is character*1
!>          = 'n': apply q  (no transpose)
!>          = 't': apply q' (transpose)
!> \endverbatim
!>
!> \param[in] m
!> \verbatim
!>          m is integer
!>          the number of rows of the matrix c. m >= 0.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is integer
!>          the number of columns of the matrix c. n >= 0.
!> \endverbatim
!>
!> \param[in] k
!> \verbatim
!>          k is integer
!>          the number of elementary reflectors whose product defines
!>          the matrix q.
!>          if side = 'l', m >= k >= 0;
!>          if side = 'r', n >= k >= 0.
!> \endverbatim
!>
!> \param[in] a
!> \verbatim
!>          a is double precision array, dimension
!>                               (lda,m) if side = 'l',
!>                               (lda,n) if side = 'r'
!>          the i-th row must contain the vector which defines the
!>          elementary reflector h(i), for i = 1,2,...,k, as returned by
!>          dgerqf in the last k rows of its array argument a.
!>          a is modified by the routine but restored on exit.
!> \endverbatim
!>
!> \param[in] lda
!> \verbatim
!>          lda is integer
!>          the leading dimension of the array a. lda >= max(1,k).
!> \endverbatim
!>
!> \param[in] tau
!> \verbatim
!>          tau is double precision array, dimension (k)
!>          tau(i) must contain the scalar factor of the elementary
!>          reflector h(i), as returned by dgerqf.
!> \endverbatim
!>
!> \param[in,out] c
!> \verbatim
!>          c is double precision array, dimension (ldc,n)
!>          on entry, the m by n matrix c.
!>          on exit, c is overwritten by q*c or q**t*c or c*q**t or c*q.
!> \endverbatim
!>
!> \param[in] ldc
!> \verbatim
!>          ldc is integer
!>          the leading dimension of the array c. ldc >= max(1,m).
!> \endverbatim
!>
!> \param[out] work
!> \verbatim
!>          work is double precision array, dimension
!>                                   (n) if side = 'l',
!>                                   (m) if side = 'r'
!> \endverbatim
!>
!> \param[out] info
!> \verbatim
!>          info is integer
!>          = 0: successful exit
!>          < 0: if info = -i, the i-th argument had an illegal value
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup unmr2
!
!  =====================================================================
subroutine dormr2( side, trans, m, n, k, a, lda, tau, c, ldc, work, info )
implicit none
!
!  -- lapack computational routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
character, intent(in) :: side !! apply the transformation from the left (l) or right (r).
character, intent(in) :: trans !! apply the transformation directly (n) or transposed (t).
integer, intent(out) :: info !! zero success; negative invalid argument or positive numerical failure.
integer, intent(in) :: k !! the number of elementary reflectors whose product defines the matrix
!! q. if side = 'l', m >= k >= 0;
integer, intent(in) :: lda !! the leading dimension of the array a. lda >= max(1,k).
integer, intent(in) :: ldc !! the leading dimension of the array c. ldc >= max(1,m).
integer, intent(in) :: m !! the number of rows of the matrix c. m >= 0.
integer, intent(in) :: n !! the number of columns of the matrix c. n >= 0.
!     ..
!     .. array arguments ..
real(dp), intent(inout) :: a( lda, * ) !! (lda,m) if side = 'l', (lda,n) if side = 'r' the i-th
!! row must contain the vector which defines the
real(dp), intent(inout) :: c( ldc, * ) !! on entry, the m by n matrix c. on exit, c is overwritten
!! by q*c or q**t*c or c*q**t or c*q.
real(dp), intent(in) :: tau( * ) !! tau(i) must contain the scalar factor of the elementary
!! reflector h(i), as returned by dgerqf.
real(dp), intent(out) :: work( * ) !! (n) if side = 'l', (m) if side = 'r'
!     ..
!
!  =====================================================================
!
!     .. parameters ..
real(dp) :: one
parameter          ( one = 1.0d+0 )
!     ..
!     .. local scalars ..
logical :: left, notran
integer :: i, i1, i2, i3, mi, ni, nq
real(dp) :: aii
!     ..
!     .. external functions ..
!     ..
!     .. external subroutines ..
!     ..
!     .. intrinsic functions ..
intrinsic          max
!     ..
!     .. executable statements ..
!
!     test the input arguments
!
info = 0
left = lsame( side, 'l' )
notran = lsame( trans, 'n' )
!
!     nq is the order of q
!
if( left ) then
nq = m
else
nq = n
end if
if( .not.left .and. .not.lsame( side, 'r' ) ) then
info = -1
else if( .not.notran .and. .not.lsame( trans, 't' ) ) then
info = -2
else if( m.lt.0 ) then
info = -3
else if( n.lt.0 ) then
info = -4
else if( k.lt.0 .or. k.gt.nq ) then
info = -5
else if( lda.lt.max( 1, k ) ) then
info = -7
else if( ldc.lt.max( 1, m ) ) then
info = -10
end if
if( info.ne.0 ) then
call xerbla( 'dormr2', -info )
return
end if
!
!     quick return if possible
!
if( m.eq.0 .or. n.eq.0 .or. k.eq.0 ) return
!
if( ( left .and. .not.notran ) .or. ( .not.left .and. notran ) ) then
i1 = 1
i2 = k
i3 = 1
else
i1 = k
i2 = 1
i3 = -1
end if
!
if( left ) then
ni = n
else
mi = m
end if
!
do 10 i = i1, i2, i3
if( left ) then
!
!           h(i) is applied to c(1:m-k+i,1:n)
!
mi = m - k + i
else
!
!           h(i) is applied to c(1:m,1:n-k+i)
!
ni = n - k + i
end if
!
!        apply h(i)
!
aii = a( i, nq-k+i )
a( i, nq-k+i ) = one
call dlarf( side, mi, ni, a( i, 1 ), lda, tau( i ), c, ldc, work )
a( i, nq-k+i ) = aii
10 continue
return
!
!     end of dormr2
!
end subroutine dormr2
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dtgsy2.f
!> \brief \b dtgsy2 solves the generalized sylvester equation (unblocked algorithm).
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download dtgsy2 + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtgsy2.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtgsy2.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtgsy2.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       subroutine dtgsy2( trans, ijob, m, n, a, lda, b, ldb, c, ldc, d,
!                          ldd, e, lde, f, ldf, scale, rdsum, rdscal,
!                          iwork, pq, info )
!
!       .. scalar arguments ..
!       character          trans
!       integer            ijob, info, lda, ldb, ldc, ldd, lde, ldf, m, n,
!      $                   pq
!       double precision   rdscal, rdsum, scale
!       ..
!       .. array arguments ..
!       integer            iwork( * )
!       double precision   a( lda, * ), b( ldb, * ), c( ldc, * ),
!      $                   d( ldd, * ), e( lde, * ), f( ldf, * )
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dtgsy2 solves the generalized sylvester equation:
!>
!>             a * r - l * b = scale * c                (1)
!>             d * r - l * e = scale * f,
!>
!> using level 1 and 2 blas. where r and l are unknown m-by-n matrices,
!> (a, d), (b, e) and (c, f) are given matrix pairs of size m-by-m,
!> n-by-n and m-by-n, respectively, with real entries. (a, d) and (b, e)
!> must be in generalized schur canonical form, i.e. a, b are upper
!> quasi triangular and d, e are upper triangular. the solution (r, l)
!> overwrites (c, f). 0 <= scale <= 1 is an output scaling factor
!> chosen to avoid overflow.
!>
!> in matrix notation solving equation (1) corresponds to solve
!> z*x = scale*b, where z is defined as
!>
!>        z = [ kron(in, a)  -kron(b**t, im) ]             (2)
!>            [ kron(in, d)  -kron(e**t, im) ],
!>
!> ik is the identity matrix of size k and x**t is the transpose of x.
!> kron(x, y) is the kronecker product between the matrices x and y.
!> in the process of solving (1), we solve a number of such systems
!> where dim(in), dim(in) = 1 or 2.
!>
!> if trans = 't', solve the transposed system z**t*y = scale*b for y,
!> which is equivalent to solve for r and l in
!>
!>             a**t * r  + d**t * l   = scale * c           (3)
!>             r  * b**t + l  * e**t  = scale * -f
!>
!> this case is used to compute an estimate of dif[(a, d), (b, e)] =
!> sigma_min(z) using reverse communication with dlacon.
!>
!> dtgsy2 also (ijob >= 1) contributes to the computation in dtgsyl
!> of an upper bound on the separation between to matrix pairs. then
!> the input (a, d), (b, e) are sub-pencils of the matrix pair in
!> dtgsyl. see dtgsyl for details.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] trans
!> \verbatim
!>          trans is character*1
!>          = 'n': solve the generalized sylvester equation (1).
!>          = 't': solve the 'transposed' system (3).
!> \endverbatim
!>
!> \param[in] ijob
!> \verbatim
!>          ijob is integer
!>          specifies what kind of functionality to be performed.
!>          = 0: solve (1) only.
!>          = 1: a contribution from this subsystem to a frobenius
!>               norm-based estimate of the separation between two matrix
!>               pairs is computed. (look ahead strategy is used).
!>          = 2: a contribution from this subsystem to a frobenius
!>               norm-based estimate of the separation between two matrix
!>               pairs is computed. (dgecon on sub-systems is used.)
!>          not referenced if trans = 't'.
!> \endverbatim
!>
!> \param[in] m
!> \verbatim
!>          m is integer
!>          on entry, m specifies the order of a and d, and the row
!>          dimension of c, f, r and l.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is integer
!>          on entry, n specifies the order of b and e, and the column
!>          dimension of c, f, r and l.
!> \endverbatim
!>
!> \param[in] a
!> \verbatim
!>          a is double precision array, dimension (lda, m)
!>          on entry, a contains an upper quasi triangular matrix.
!> \endverbatim
!>
!> \param[in] lda
!> \verbatim
!>          lda is integer
!>          the leading dimension of the matrix a. lda >= max(1, m).
!> \endverbatim
!>
!> \param[in] b
!> \verbatim
!>          b is double precision array, dimension (ldb, n)
!>          on entry, b contains an upper quasi triangular matrix.
!> \endverbatim
!>
!> \param[in] ldb
!> \verbatim
!>          ldb is integer
!>          the leading dimension of the matrix b. ldb >= max(1, n).
!> \endverbatim
!>
!> \param[in,out] c
!> \verbatim
!>          c is double precision array, dimension (ldc, n)
!>          on entry, c contains the right-hand-side of the first matrix
!>          equation in (1).
!>          on exit, if ijob = 0, c has been overwritten by the
!>          solution r.
!> \endverbatim
!>
!> \param[in] ldc
!> \verbatim
!>          ldc is integer
!>          the leading dimension of the matrix c. ldc >= max(1, m).
!> \endverbatim
!>
!> \param[in] d
!> \verbatim
!>          d is double precision array, dimension (ldd, m)
!>          on entry, d contains an upper triangular matrix.
!> \endverbatim
!>
!> \param[in] ldd
!> \verbatim
!>          ldd is integer
!>          the leading dimension of the matrix d. ldd >= max(1, m).
!> \endverbatim
!>
!> \param[in] e
!> \verbatim
!>          e is double precision array, dimension (lde, n)
!>          on entry, e contains an upper triangular matrix.
!> \endverbatim
!>
!> \param[in] lde
!> \verbatim
!>          lde is integer
!>          the leading dimension of the matrix e. lde >= max(1, n).
!> \endverbatim
!>
!> \param[in,out] f
!> \verbatim
!>          f is double precision array, dimension (ldf, n)
!>          on entry, f contains the right-hand-side of the second matrix
!>          equation in (1).
!>          on exit, if ijob = 0, f has been overwritten by the
!>          solution l.
!> \endverbatim
!>
!> \param[in] ldf
!> \verbatim
!>          ldf is integer
!>          the leading dimension of the matrix f. ldf >= max(1, m).
!> \endverbatim
!>
!> \param[out] scale
!> \verbatim
!>          scale is double precision
!>          on exit, 0 <= scale <= 1. if 0 < scale < 1, the solutions
!>          r and l (c and f on entry) will hold the solutions to a
!>          slightly perturbed system but the input matrices a, b, d and
!>          e have not been changed. if scale = 0, r and l will hold the
!>          solutions to the homogeneous system with c = f = 0. normally,
!>          scale = 1.
!> \endverbatim
!>
!> \param[in,out] rdsum
!> \verbatim
!>          rdsum is double precision
!>          on entry, the sum of squares of computed contributions to
!>          the dif-estimate under computation by dtgsyl, where the
!>          scaling factor rdscal (see below) has been factored out.
!>          on exit, the corresponding sum of squares updated with the
!>          contributions from the current sub-system.
!>          if trans = 't' rdsum is not touched.
!>          note: rdsum only makes sense when dtgsy2 is called by dtgsyl.
!> \endverbatim
!>
!> \param[in,out] rdscal
!> \verbatim
!>          rdscal is double precision
!>          on entry, scaling factor used to prevent overflow in rdsum.
!>          on exit, rdscal is updated w.r.t. the current contributions
!>          in rdsum.
!>          if trans = 't', rdscal is not touched.
!>          note: rdscal only makes sense when dtgsy2 is called by
!>                dtgsyl.
!> \endverbatim
!>
!> \param[out] iwork
!> \verbatim
!>          iwork is integer array, dimension (m+n+2)
!> \endverbatim
!>
!> \param[out] pq
!> \verbatim
!>          pq is integer
!>          on exit, the number of subsystems (of size 2-by-2, 4-by-4 and
!>          8-by-8) solved by this routine.
!> \endverbatim
!>
!> \param[out] info
!> \verbatim
!>          info is integer
!>          on exit, if info is set to
!>            =0: successful exit
!>            <0: if info = -i, the i-th argument had an illegal value.
!>            >0: the matrix pairs (a, d) and (b, e) have common or very
!>                close eigenvalues.
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup tgsy2
!
!> \par contributors:
!  ==================
!>
!>     bo kagstrom and peter poromaa, department of computing science,
!>     umea university, s-901 87 umea, sweden.
!
!  =====================================================================
subroutine dtgsy2( trans, ijob, m, n, a, lda, b, ldb, c, ldc, d, ldd, e, lde, f, ldf, scale, rdsum, &
    rdscal, iwork, pq, info )
implicit none
!
!  -- lapack auxiliary routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
character, intent(in) :: trans !! apply the transformation directly (n) or transposed (t).
integer, intent(in) :: ijob !! specifies what kind of functionality to be performed.
integer, intent(out) :: info !! on exit, if info is set to
integer, intent(in) :: lda !! the leading dimension of the matrix a. lda >= max(1, m).
integer, intent(in) :: ldb !! the leading dimension of the matrix b. ldb >= max(1, n).
integer, intent(in) :: ldc !! the leading dimension of the matrix c. ldc >= max(1, m).
integer, intent(in) :: ldd !! the leading dimension of the matrix d. ldd >= max(1, m).
integer, intent(in) :: lde !! the leading dimension of the matrix e. lde >= max(1, n).
integer, intent(in) :: ldf !! the leading dimension of the matrix f. ldf >= max(1, m).
integer, intent(in) :: m !! on entry, m specifies the order of a and d, and the row dimension of
!! c, f, r and l.
integer, intent(in) :: n !! on entry, n specifies the order of b and e, and the column dimension
!! of c, f, r and l.
integer, intent(out) :: pq !! on exit, the number of subsystems (of size 2-by-2, 4-by-4 and
!! 8-by-8) solved by this routine.
real(dp), intent(inout) :: rdscal !! on entry, scaling factor used to prevent overflow in rdsum.
!! on exit, rdscal is updated w.r.t. the current contributions
real(dp), intent(inout) :: rdsum !! on entry, the sum of squares of computed contributions to the
!! dif-estimate under computation by dtgsyl, where the
real(dp), intent(out) :: scale !! on exit, 0 <= scale <= 1. if 0 < scale < 1, the solutions r and
!! l (c and f on entry) will hold the solutions to a
!     ..
!     .. array arguments ..
integer, intent(out) :: iwork( * ) !! integer scratch and pivot indices for the active small systems.
real(dp), intent(in) :: a( lda, * ) !! on entry, a contains an upper quasi triangular matrix.
real(dp), intent(in) :: b( ldb, * ) !! on entry, b contains an upper quasi triangular matrix.
real(dp), intent(inout) :: c( ldc, * ) !! on entry, c contains the right-hand-side of the first
!! matrix equation in (1). on exit, if ijob = 0, c has been
!! overwritten by the
real(dp), intent(in) :: d( ldd, * ) !! on entry, d contains an upper triangular matrix.
real(dp), intent(in) :: e( lde, * ) !! on entry, e contains an upper triangular matrix.
real(dp), intent(inout) :: f( ldf, * ) !! on entry, f contains the right-hand-side of the second
!! matrix equation in (1). on exit, if ijob = 0, f has been
!! overwritten by the
!     ..
!
!  =====================================================================
!  replaced various illegal calls to dcopy by calls to dlaset.
!  sven hammarling, 27/5/02.
!
!     .. parameters ..
integer :: ldz
parameter          ( ldz = 8 )
real(dp) :: zero, one
parameter          ( zero = 0.0d+0, one = 1.0d+0 )
!     ..
!     .. local scalars ..
logical :: notran
integer :: i, ie, ierr, ii, is, isp1, j, je, jj, js, jsp1, k, mb, nb, p, q, zdim
real(dp) :: alpha, scaloc
!     ..
!     .. local arrays ..
integer :: ipiv( ldz ), jpiv( ldz )
real(dp) :: rhs( ldz ), z( ldz, ldz )
!     ..
!     .. external functions ..
!     ..
!     .. external subroutines ..
!     ..
!     .. intrinsic functions ..
intrinsic          max
!     ..
!     .. executable statements ..
!
!     decode and test input parameters
!
info = 0
ierr = 0
notran = lsame( trans, 'n' )
if( .not.notran .and. .not.lsame( trans, 't' ) ) then
info = -1
else if( notran ) then
if( ( ijob.lt.0 ) .or. ( ijob.gt.2 ) ) then
info = -2
end if
end if
if( info.eq.0 ) then
if( m.le.0 ) then
info = -3
else if( n.le.0 ) then
info = -4
else if( lda.lt.max( 1, m ) ) then
info = -6
else if( ldb.lt.max( 1, n ) ) then
info = -8
else if( ldc.lt.max( 1, m ) ) then
info = -10
else if( ldd.lt.max( 1, m ) ) then
info = -12
else if( lde.lt.max( 1, n ) ) then
info = -14
else if( ldf.lt.max( 1, m ) ) then
info = -16
end if
end if
if( info.ne.0 ) then
call xerbla( 'dtgsy2', -info )
return
end if
!
!     determine block structure of a
!
pq = 0
p = 0
i = 1
10 continue
if( i.gt.m ) go to 20
p = p + 1
iwork( p ) = i
if( i.eq.m ) go to 20
if( a( i+1, i ).ne.zero ) then
i = i + 2
else
i = i + 1
end if
go to 10
20 continue
iwork( p+1 ) = m + 1
!
!     determine block structure of b
!
q = p + 1
j = 1
30 continue
if( j.gt.n ) go to 40
q = q + 1
iwork( q ) = j
if( j.eq.n ) go to 40
if( b( j+1, j ).ne.zero ) then
j = j + 2
else
j = j + 1
end if
go to 30
40 continue
iwork( q+1 ) = n + 1
pq = p*( q-p-1 )
!
if( notran ) then
!
!        solve (i, j) - subsystem
!           a(i, i) * r(i, j) - l(i, j) * b(j, j) = c(i, j)
!           d(i, i) * r(i, j) - l(i, j) * e(j, j) = f(i, j)
!        for i = p, p - 1, ..., 1; j = 1, 2, ..., q
!
scale = one
scaloc = one
do 120 j = p + 2, q
js = iwork( j )
jsp1 = js + 1
je = iwork( j+1 ) - 1
nb = je - js + 1
do 110 i = p, 1, -1
!
is = iwork( i )
isp1 = is + 1
ie = iwork( i+1 ) - 1
mb = ie - is + 1
zdim = mb*nb*2
!
if( ( mb.eq.1 ) .and. ( nb.eq.1 ) ) then
!
!                 build a 2-by-2 system z * x = rhs
!
z( 1, 1 ) = a( is, is )
z( 2, 1 ) = d( is, is )
z( 1, 2 ) = -b( js, js )
z( 2, 2 ) = -e( js, js )
!
!                 set up right hand side(s)
!
rhs( 1 ) = c( is, js )
rhs( 2 ) = f( is, js )
!
!                 solve z * x = rhs
!
call dgetc2( zdim, z, ldz, ipiv, jpiv, ierr )
if( ierr.gt.0 ) info = ierr
!
if( ijob.eq.0 ) then
call dgesc2( zdim, z, ldz, rhs, ipiv, jpiv, scaloc )
if( scaloc.ne.one ) then
do 50 k = 1, n
call dscal( m, scaloc, c( 1, k ), 1 )
call dscal( m, scaloc, f( 1, k ), 1 )
50 continue
scale = scale*scaloc
end if
else
call dlatdf( ijob, zdim, z, ldz, rhs, rdsum, rdscal, ipiv, jpiv )
end if
!
!                 unpack solution vector(s)
!
c( is, js ) = rhs( 1 )
f( is, js ) = rhs( 2 )
!
!                 substitute r(i, j) and l(i, j) into remaining
!                 equation.
!
if( i.gt.1 ) then
alpha = -rhs( 1 )
call daxpy( is-1, alpha, a( 1, is ), 1, c( 1, js ), 1 )
call daxpy( is-1, alpha, d( 1, is ), 1, f( 1, js ), 1 )
end if
if( j.lt.q ) then
call daxpy( n-je, rhs( 2 ), b( js, je+1 ), ldb, c( is, je+1 ), ldc )
call daxpy( n-je, rhs( 2 ), e( js, je+1 ), lde, f( is, je+1 ), ldf )
end if
!
else if( ( mb.eq.1 ) .and. ( nb.eq.2 ) ) then
!
!                 build a 4-by-4 system z * x = rhs
!
z( 1, 1 ) = a( is, is )
z( 2, 1 ) = zero
z( 3, 1 ) = d( is, is )
z( 4, 1 ) = zero
!
z( 1, 2 ) = zero
z( 2, 2 ) = a( is, is )
z( 3, 2 ) = zero
z( 4, 2 ) = d( is, is )
!
z( 1, 3 ) = -b( js, js )
z( 2, 3 ) = -b( js, jsp1 )
z( 3, 3 ) = -e( js, js )
z( 4, 3 ) = -e( js, jsp1 )
!
z( 1, 4 ) = -b( jsp1, js )
z( 2, 4 ) = -b( jsp1, jsp1 )
z( 3, 4 ) = zero
z( 4, 4 ) = -e( jsp1, jsp1 )
!
!                 set up right hand side(s)
!
rhs( 1 ) = c( is, js )
rhs( 2 ) = c( is, jsp1 )
rhs( 3 ) = f( is, js )
rhs( 4 ) = f( is, jsp1 )
!
!                 solve z * x = rhs
!
call dgetc2( zdim, z, ldz, ipiv, jpiv, ierr )
if( ierr.gt.0 ) info = ierr
!
if( ijob.eq.0 ) then
call dgesc2( zdim, z, ldz, rhs, ipiv, jpiv, scaloc )
if( scaloc.ne.one ) then
do 60 k = 1, n
call dscal( m, scaloc, c( 1, k ), 1 )
call dscal( m, scaloc, f( 1, k ), 1 )
60 continue
scale = scale*scaloc
end if
else
call dlatdf( ijob, zdim, z, ldz, rhs, rdsum, rdscal, ipiv, jpiv )
end if
!
!                 unpack solution vector(s)
!
c( is, js ) = rhs( 1 )
c( is, jsp1 ) = rhs( 2 )
f( is, js ) = rhs( 3 )
f( is, jsp1 ) = rhs( 4 )
!
!                 substitute r(i, j) and l(i, j) into remaining
!                 equation.
!
if( i.gt.1 ) then
call dger( is-1, nb, -one, a( 1, is ), 1, rhs( 1 ), 1, c( 1, js ), ldc )
call dger( is-1, nb, -one, d( 1, is ), 1, rhs( 1 ), 1, f( 1, js ), ldf )
end if
if( j.lt.q ) then
call daxpy( n-je, rhs( 3 ), b( js, je+1 ), ldb, c( is, je+1 ), ldc )
call daxpy( n-je, rhs( 3 ), e( js, je+1 ), lde, f( is, je+1 ), ldf )
call daxpy( n-je, rhs( 4 ), b( jsp1, je+1 ), ldb, c( is, je+1 ), ldc )
call daxpy( n-je, rhs( 4 ), e( jsp1, je+1 ), lde, f( is, je+1 ), ldf )
end if
!
else if( ( mb.eq.2 ) .and. ( nb.eq.1 ) ) then
!
!                 build a 4-by-4 system z * x = rhs
!
z( 1, 1 ) = a( is, is )
z( 2, 1 ) = a( isp1, is )
z( 3, 1 ) = d( is, is )
z( 4, 1 ) = zero
!
z( 1, 2 ) = a( is, isp1 )
z( 2, 2 ) = a( isp1, isp1 )
z( 3, 2 ) = d( is, isp1 )
z( 4, 2 ) = d( isp1, isp1 )
!
z( 1, 3 ) = -b( js, js )
z( 2, 3 ) = zero
z( 3, 3 ) = -e( js, js )
z( 4, 3 ) = zero
!
z( 1, 4 ) = zero
z( 2, 4 ) = -b( js, js )
z( 3, 4 ) = zero
z( 4, 4 ) = -e( js, js )
!
!                 set up right hand side(s)
!
rhs( 1 ) = c( is, js )
rhs( 2 ) = c( isp1, js )
rhs( 3 ) = f( is, js )
rhs( 4 ) = f( isp1, js )
!
!                 solve z * x = rhs
!
call dgetc2( zdim, z, ldz, ipiv, jpiv, ierr )
if( ierr.gt.0 ) info = ierr
if( ijob.eq.0 ) then
call dgesc2( zdim, z, ldz, rhs, ipiv, jpiv, scaloc )
if( scaloc.ne.one ) then
do 70 k = 1, n
call dscal( m, scaloc, c( 1, k ), 1 )
call dscal( m, scaloc, f( 1, k ), 1 )
70 continue
scale = scale*scaloc
end if
else
call dlatdf( ijob, zdim, z, ldz, rhs, rdsum, rdscal, ipiv, jpiv )
end if
!
!                 unpack solution vector(s)
!
c( is, js ) = rhs( 1 )
c( isp1, js ) = rhs( 2 )
f( is, js ) = rhs( 3 )
f( isp1, js ) = rhs( 4 )
!
!                 substitute r(i, j) and l(i, j) into remaining
!                 equation.
!
if( i.gt.1 ) then
call dgemv( 'n', is-1, mb, -one, a( 1, is ), lda, rhs( 1 ), 1, one, c( 1, js ), 1 )
call dgemv( 'n', is-1, mb, -one, d( 1, is ), ldd, rhs( 1 ), 1, one, f( 1, js ), 1 )
end if
if( j.lt.q ) then
call dger( mb, n-je, one, rhs( 3 ), 1, b( js, je+1 ), ldb, c( is, je+1 ), ldc )
call dger( mb, n-je, one, rhs( 3 ), 1, e( js, je+1 ), lde, f( is, je+1 ), ldf )
end if
!
else if( ( mb.eq.2 ) .and. ( nb.eq.2 ) ) then
!
!                 build an 8-by-8 system z * x = rhs
!
call dlaset( 'f', ldz, ldz, zero, zero, z, ldz )
!
z( 1, 1 ) = a( is, is )
z( 2, 1 ) = a( isp1, is )
z( 5, 1 ) = d( is, is )
!
z( 1, 2 ) = a( is, isp1 )
z( 2, 2 ) = a( isp1, isp1 )
z( 5, 2 ) = d( is, isp1 )
z( 6, 2 ) = d( isp1, isp1 )
!
z( 3, 3 ) = a( is, is )
z( 4, 3 ) = a( isp1, is )
z( 7, 3 ) = d( is, is )
!
z( 3, 4 ) = a( is, isp1 )
z( 4, 4 ) = a( isp1, isp1 )
z( 7, 4 ) = d( is, isp1 )
z( 8, 4 ) = d( isp1, isp1 )
!
z( 1, 5 ) = -b( js, js )
z( 3, 5 ) = -b( js, jsp1 )
z( 5, 5 ) = -e( js, js )
z( 7, 5 ) = -e( js, jsp1 )
!
z( 2, 6 ) = -b( js, js )
z( 4, 6 ) = -b( js, jsp1 )
z( 6, 6 ) = -e( js, js )
z( 8, 6 ) = -e( js, jsp1 )
!
z( 1, 7 ) = -b( jsp1, js )
z( 3, 7 ) = -b( jsp1, jsp1 )
z( 7, 7 ) = -e( jsp1, jsp1 )
!
z( 2, 8 ) = -b( jsp1, js )
z( 4, 8 ) = -b( jsp1, jsp1 )
z( 8, 8 ) = -e( jsp1, jsp1 )
!
!                 set up right hand side(s)
!
k = 1
ii = mb*nb + 1
do 80 jj = 0, nb - 1
call dcopy( mb, c( is, js+jj ), 1, rhs( k ), 1 )
call dcopy( mb, f( is, js+jj ), 1, rhs( ii ), 1 )
k = k + mb
ii = ii + mb
80 continue
!
!                 solve z * x = rhs
!
call dgetc2( zdim, z, ldz, ipiv, jpiv, ierr )
if( ierr.gt.0 ) info = ierr
if( ijob.eq.0 ) then
call dgesc2( zdim, z, ldz, rhs, ipiv, jpiv, scaloc )
if( scaloc.ne.one ) then
do 90 k = 1, n
call dscal( m, scaloc, c( 1, k ), 1 )
call dscal( m, scaloc, f( 1, k ), 1 )
90 continue
scale = scale*scaloc
end if
else
call dlatdf( ijob, zdim, z, ldz, rhs, rdsum, rdscal, ipiv, jpiv )
end if
!
!                 unpack solution vector(s)
!
k = 1
ii = mb*nb + 1
do 100 jj = 0, nb - 1
call dcopy( mb, rhs( k ), 1, c( is, js+jj ), 1 )
call dcopy( mb, rhs( ii ), 1, f( is, js+jj ), 1 )
k = k + mb
ii = ii + mb
100 continue
!
!                 substitute r(i, j) and l(i, j) into remaining
!                 equation.
!
if( i.gt.1 ) then
call dgemm( 'n', 'n', is-1, nb, mb, -one, a( 1, is ), lda, rhs( 1 ), mb, one, c( 1, js ), ldc )
call dgemm( 'n', 'n', is-1, nb, mb, -one, d( 1, is ), ldd, rhs( 1 ), mb, one, f( 1, js ), ldf )
end if
if( j.lt.q ) then
k = mb*nb + 1
call dgemm( 'n', 'n', mb, n-je, nb, one, rhs( k ), mb, b( js, je+1 ), ldb, one, c( is, je+1 ), ldc )
call dgemm( 'n', 'n', mb, n-je, nb, one, rhs( k ), mb, e( js, je+1 ), lde, one, f( is, je+1 ), ldf )
end if
!
end if
!
110 continue
120 continue
else
!
!        solve (i, j) - subsystem
!             a(i, i)**t * r(i, j) + d(i, i)**t * l(j, j)  =  c(i, j)
!             r(i, i)  * b(j, j) + l(i, j)  * e(j, j)  = -f(i, j)
!        for i = 1, 2, ..., p, j = q, q - 1, ..., 1
!
scale = one
scaloc = one
do 200 i = 1, p
!
is = iwork( i )
isp1 = is + 1
ie = iwork ( i+1 ) - 1
mb = ie - is + 1
do 190 j = q, p + 2, -1
!
js = iwork( j )
jsp1 = js + 1
je = iwork( j+1 ) - 1
nb = je - js + 1
zdim = mb*nb*2
if( ( mb.eq.1 ) .and. ( nb.eq.1 ) ) then
!
!                 build a 2-by-2 system z**t * x = rhs
!
z( 1, 1 ) = a( is, is )
z( 2, 1 ) = -b( js, js )
z( 1, 2 ) = d( is, is )
z( 2, 2 ) = -e( js, js )
!
!                 set up right hand side(s)
!
rhs( 1 ) = c( is, js )
rhs( 2 ) = f( is, js )
!
!                 solve z**t * x = rhs
!
call dgetc2( zdim, z, ldz, ipiv, jpiv, ierr )
if( ierr.gt.0 ) info = ierr
!
call dgesc2( zdim, z, ldz, rhs, ipiv, jpiv, scaloc )
if( scaloc.ne.one ) then
do 130 k = 1, n
call dscal( m, scaloc, c( 1, k ), 1 )
call dscal( m, scaloc, f( 1, k ), 1 )
130 continue
scale = scale*scaloc
end if
!
!                 unpack solution vector(s)
!
c( is, js ) = rhs( 1 )
f( is, js ) = rhs( 2 )
!
!                 substitute r(i, j) and l(i, j) into remaining
!                 equation.
!
if( j.gt.p+2 ) then
alpha = rhs( 1 )
call daxpy( js-1, alpha, b( 1, js ), 1, f( is, 1 ), ldf )
alpha = rhs( 2 )
call daxpy( js-1, alpha, e( 1, js ), 1, f( is, 1 ), ldf )
end if
if( i.lt.p ) then
alpha = -rhs( 1 )
call daxpy( m-ie, alpha, a( is, ie+1 ), lda, c( ie+1, js ), 1 )
alpha = -rhs( 2 )
call daxpy( m-ie, alpha, d( is, ie+1 ), ldd, c( ie+1, js ), 1 )
end if
!
else if( ( mb.eq.1 ) .and. ( nb.eq.2 ) ) then
!
!                 build a 4-by-4 system z**t * x = rhs
!
z( 1, 1 ) = a( is, is )
z( 2, 1 ) = zero
z( 3, 1 ) = -b( js, js )
z( 4, 1 ) = -b( jsp1, js )
!
z( 1, 2 ) = zero
z( 2, 2 ) = a( is, is )
z( 3, 2 ) = -b( js, jsp1 )
z( 4, 2 ) = -b( jsp1, jsp1 )
!
z( 1, 3 ) = d( is, is )
z( 2, 3 ) = zero
z( 3, 3 ) = -e( js, js )
z( 4, 3 ) = zero
!
z( 1, 4 ) = zero
z( 2, 4 ) = d( is, is )
z( 3, 4 ) = -e( js, jsp1 )
z( 4, 4 ) = -e( jsp1, jsp1 )
!
!                 set up right hand side(s)
!
rhs( 1 ) = c( is, js )
rhs( 2 ) = c( is, jsp1 )
rhs( 3 ) = f( is, js )
rhs( 4 ) = f( is, jsp1 )
!
!                 solve z**t * x = rhs
!
call dgetc2( zdim, z, ldz, ipiv, jpiv, ierr )
if( ierr.gt.0 ) info = ierr
call dgesc2( zdim, z, ldz, rhs, ipiv, jpiv, scaloc )
if( scaloc.ne.one ) then
do 140 k = 1, n
call dscal( m, scaloc, c( 1, k ), 1 )
call dscal( m, scaloc, f( 1, k ), 1 )
140 continue
scale = scale*scaloc
end if
!
!                 unpack solution vector(s)
!
c( is, js ) = rhs( 1 )
c( is, jsp1 ) = rhs( 2 )
f( is, js ) = rhs( 3 )
f( is, jsp1 ) = rhs( 4 )
!
!                 substitute r(i, j) and l(i, j) into remaining
!                 equation.
!
if( j.gt.p+2 ) then
call daxpy( js-1, rhs( 1 ), b( 1, js ), 1, f( is, 1 ), ldf )
call daxpy( js-1, rhs( 2 ), b( 1, jsp1 ), 1, f( is, 1 ), ldf )
call daxpy( js-1, rhs( 3 ), e( 1, js ), 1, f( is, 1 ), ldf )
call daxpy( js-1, rhs( 4 ), e( 1, jsp1 ), 1, f( is, 1 ), ldf )
end if
if( i.lt.p ) then
call dger( m-ie, nb, -one, a( is, ie+1 ), lda, rhs( 1 ), 1, c( ie+1, js ), ldc )
call dger( m-ie, nb, -one, d( is, ie+1 ), ldd, rhs( 3 ), 1, c( ie+1, js ), ldc )
end if
!
else if( ( mb.eq.2 ) .and. ( nb.eq.1 ) ) then
!
!                 build a 4-by-4 system z**t * x = rhs
!
z( 1, 1 ) = a( is, is )
z( 2, 1 ) = a( is, isp1 )
z( 3, 1 ) = -b( js, js )
z( 4, 1 ) = zero
!
z( 1, 2 ) = a( isp1, is )
z( 2, 2 ) = a( isp1, isp1 )
z( 3, 2 ) = zero
z( 4, 2 ) = -b( js, js )
!
z( 1, 3 ) = d( is, is )
z( 2, 3 ) = d( is, isp1 )
z( 3, 3 ) = -e( js, js )
z( 4, 3 ) = zero
!
z( 1, 4 ) = zero
z( 2, 4 ) = d( isp1, isp1 )
z( 3, 4 ) = zero
z( 4, 4 ) = -e( js, js )
!
!                 set up right hand side(s)
!
rhs( 1 ) = c( is, js )
rhs( 2 ) = c( isp1, js )
rhs( 3 ) = f( is, js )
rhs( 4 ) = f( isp1, js )
!
!                 solve z**t * x = rhs
!
call dgetc2( zdim, z, ldz, ipiv, jpiv, ierr )
if( ierr.gt.0 ) info = ierr
!
call dgesc2( zdim, z, ldz, rhs, ipiv, jpiv, scaloc )
if( scaloc.ne.one ) then
do 150 k = 1, n
call dscal( m, scaloc, c( 1, k ), 1 )
call dscal( m, scaloc, f( 1, k ), 1 )
150 continue
scale = scale*scaloc
end if
!
!                 unpack solution vector(s)
!
c( is, js ) = rhs( 1 )
c( isp1, js ) = rhs( 2 )
f( is, js ) = rhs( 3 )
f( isp1, js ) = rhs( 4 )
!
!                 substitute r(i, j) and l(i, j) into remaining
!                 equation.
!
if( j.gt.p+2 ) then
call dger( mb, js-1, one, rhs( 1 ), 1, b( 1, js ), 1, f( is, 1 ), ldf )
call dger( mb, js-1, one, rhs( 3 ), 1, e( 1, js ), 1, f( is, 1 ), ldf )
end if
if( i.lt.p ) then
call dgemv( 't', mb, m-ie, -one, a( is, ie+1 ), lda, rhs( 1 ), 1, one, c( ie+1, js ), 1 )
call dgemv( 't', mb, m-ie, -one, d( is, ie+1 ), ldd, rhs( 3 ), 1, one, c( ie+1, js ), 1 )
end if
!
else if( ( mb.eq.2 ) .and. ( nb.eq.2 ) ) then
!
!                 build an 8-by-8 system z**t * x = rhs
!
call dlaset( 'f', ldz, ldz, zero, zero, z, ldz )
!
z( 1, 1 ) = a( is, is )
z( 2, 1 ) = a( is, isp1 )
z( 5, 1 ) = -b( js, js )
z( 7, 1 ) = -b( jsp1, js )
!
z( 1, 2 ) = a( isp1, is )
z( 2, 2 ) = a( isp1, isp1 )
z( 6, 2 ) = -b( js, js )
z( 8, 2 ) = -b( jsp1, js )
!
z( 3, 3 ) = a( is, is )
z( 4, 3 ) = a( is, isp1 )
z( 5, 3 ) = -b( js, jsp1 )
z( 7, 3 ) = -b( jsp1, jsp1 )
!
z( 3, 4 ) = a( isp1, is )
z( 4, 4 ) = a( isp1, isp1 )
z( 6, 4 ) = -b( js, jsp1 )
z( 8, 4 ) = -b( jsp1, jsp1 )
!
z( 1, 5 ) = d( is, is )
z( 2, 5 ) = d( is, isp1 )
z( 5, 5 ) = -e( js, js )
!
z( 2, 6 ) = d( isp1, isp1 )
z( 6, 6 ) = -e( js, js )
!
z( 3, 7 ) = d( is, is )
z( 4, 7 ) = d( is, isp1 )
z( 5, 7 ) = -e( js, jsp1 )
z( 7, 7 ) = -e( jsp1, jsp1 )
!
z( 4, 8 ) = d( isp1, isp1 )
z( 6, 8 ) = -e( js, jsp1 )
z( 8, 8 ) = -e( jsp1, jsp1 )
!
!                 set up right hand side(s)
!
k = 1
ii = mb*nb + 1
do 160 jj = 0, nb - 1
call dcopy( mb, c( is, js+jj ), 1, rhs( k ), 1 )
call dcopy( mb, f( is, js+jj ), 1, rhs( ii ), 1 )
k = k + mb
ii = ii + mb
160 continue
!
!
!                 solve z**t * x = rhs
!
call dgetc2( zdim, z, ldz, ipiv, jpiv, ierr )
if( ierr.gt.0 ) info = ierr
!
call dgesc2( zdim, z, ldz, rhs, ipiv, jpiv, scaloc )
if( scaloc.ne.one ) then
do 170 k = 1, n
call dscal( m, scaloc, c( 1, k ), 1 )
call dscal( m, scaloc, f( 1, k ), 1 )
170 continue
scale = scale*scaloc
end if
!
!                 unpack solution vector(s)
!
k = 1
ii = mb*nb + 1
do 180 jj = 0, nb - 1
call dcopy( mb, rhs( k ), 1, c( is, js+jj ), 1 )
call dcopy( mb, rhs( ii ), 1, f( is, js+jj ), 1 )
k = k + mb
ii = ii + mb
180 continue
!
!                 substitute r(i, j) and l(i, j) into remaining
!                 equation.
!
if( j.gt.p+2 ) then
call dgemm( 'n', 't', mb, js-1, nb, one, c( is, js ), ldc, b( 1, js ), ldb, one, f( is, 1 ), ldf )
call dgemm( 'n', 't', mb, js-1, nb, one, f( is, js ), ldf, e( 1, js ), lde, one, f( is, 1 ), ldf )
end if
if( i.lt.p ) then
call dgemm( 't', 'n', m-ie, nb, mb, -one, a( is, ie+1 ), lda, c( is, js ), ldc, one, c( ie+1, js ), ldc )
call dgemm( 't', 'n', m-ie, nb, mb, -one, d( is, ie+1 ), ldd, f( is, js ), ldf, one, c( ie+1, js ), ldc )
end if
!
end if
!
190 continue
200 continue
!
end if
return
!
!     end of dtgsy2
!
end subroutine dtgsy2
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlaisnan.f
!> \brief \b dlaisnan tests input for nan by comparing two arguments for inequality.
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download dlaisnan + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaisnan.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaisnan.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaisnan.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       logical function dlaisnan( din1, din2 )
!
!       .. scalar arguments ..
!       double precision, intent(in) :: din1, din2
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> this routine is not for general use.  it exists solely to avoid
!> over-optimization in disnan.
!>
!> dlaisnan checks for nans by comparing its two arguments for
!> inequality.  nan is the only floating-point value where nan != nan
!> returns .true.  to check for nans, pass the same variable as both
!> arguments.
!>
!> a compiler must assume that the two arguments are
!> not the same variable, and the test will not be optimized away.
!> interprocedural or whole-program optimization may delete this
!> test.  the isnan functions will be replaced by the correct
!> fortran 03 intrinsic once the intrinsic is widely available.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] din1
!> \verbatim
!>          din1 is double precision
!> \endverbatim
!>
!> \param[in] din2
!> \verbatim
!>          din2 is double precision
!>          two numbers to compare for inequality.
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup laisnan
!
!  =====================================================================
logical function dlaisnan( din1, din2 ) result(lapack_result)
implicit none
!
!  -- lapack auxiliary routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
real(dp), intent(in) :: din1, din2
!     ..
!
!  =====================================================================
!
!  .. executable statements ..
lapack_result = (din1.ne.din2)
return
end function dlaisnan
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlarf1f.f
!> \brief \b dlarf1f applies an elementary reflector to a general rectangular
!              matrix assuming v(1) = 1.
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download dlarf1f + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarf1f.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarf1f.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarf1f.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       subroutine dlarf1f( side, m, n, v, incv, tau, c, ldc, work )
!
!       .. scalar arguments ..
!       character          side
!       integer            incv, ldc, m, n
!       double precision   tau
!       ..
!       .. array arguments ..
!       double precision   c( ldc, * ), v( * ), work( * )
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dlarf1f applies a real elementary reflector h to a real m by n matrix
!> c, from either the left or the right. h is represented in the form
!>
!>       h = i - tau * v * v**t
!>
!> where tau is a real scalar and v is a real vector.
!>
!> if tau = 0, then h is taken to be the unit matrix.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] side
!> \verbatim
!>          side is character*1
!>          = 'l': form  h * c
!>          = 'r': form  c * h
!> \endverbatim
!>
!> \param[in] m
!> \verbatim
!>          m is integer
!>          the number of rows of the matrix c.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is integer
!>          the number of columns of the matrix c.
!> \endverbatim
!>
!> \param[in] v
!> \verbatim
!>          v is double precision array, dimension
!>                     (1 + (m-1)*abs(incv)) if side = 'l'
!>                  or (1 + (n-1)*abs(incv)) if side = 'r'
!>          the vector v in the representation of h. v is not used if
!>          tau = 0. v(1) is not referenced or modified.
!> \endverbatim
!>
!> \param[in] incv
!> \verbatim
!>          incv is integer
!>          the increment between elements of v. incv <> 0.
!> \endverbatim
!>
!> \param[in] tau
!> \verbatim
!>          tau is double precision
!>          the value tau in the representation of h.
!> \endverbatim
!>
!> \param[in,out] c
!> \verbatim
!>          c is double precision array, dimension (ldc,n)
!>          on entry, the m by n matrix c.
!>          on exit, c is overwritten by the matrix h * c if side = 'l',
!>          or c * h if side = 'r'.
!> \endverbatim
!>
!> \param[in] ldc
!> \verbatim
!>          ldc is integer
!>          the leading dimension of the array c. ldc >= max(1,m).
!> \endverbatim
!>
!> \param[out] work
!> \verbatim
!>          work is double precision array, dimension
!>                         (n) if side = 'l'
!>                      or (m) if side = 'r'
!> \endverbatim
!
!  to take advantage of the fact that v(1) = 1, we do the following
!     v = [ 1 v_2 ]**t
!     if side='l'
!           |-----|
!           | c_1 |
!        c =| c_2 |
!           |-----|
!        c_1\in\mathbb{r}^{1\times n}, c_2\in\mathbb{r}^{m-1\times n}
!        so we compute:
!        c = hc   = (i - \tau vv**t)c
!                 = c - \tau vv**t c
!        w = c**t v  = [ c_1**t c_2**t ] [ 1 v_2 ]**t
!                    = c_1**t + c_2**t v ( dgemm then daxpy )
!        c  = c - \tau vv**t c
!           = c - \tau vw**t
!        giving us   c_1 = c_1 - \tau w**t ( daxpy )
!                 and
!                    c_2 = c_2 - \tau v_2w**t ( dger )
!     if side='r'
!
!        c = [ c_1 c_2 ]
!        c_1\in\mathbb{r}^{m\times 1}, c_2\in\mathbb{r}^{m\times n-1}
!        so we compute:
!        c = ch   = c(i - \tau vv**t)
!                 = c - \tau cvv**t
!
!        w = cv   = [ c_1 c_2 ] [ 1 v_2 ]**t
!                 = c_1 + c_2v_2 ( dgemm then daxpy )
!        c  = c - \tau cvv**t
!           = c - \tau wv**t
!        giving us   c_1 = c_1 - \tau w ( daxpy )
!                 and
!                    c_2 = c_2 - \tau wv_2**t ( dger )
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup larf
!
!  =====================================================================
subroutine dlarf1f( side, m, n, v, incv, tau, c, ldc, work )
implicit none
!
!  -- lapack auxiliary routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
character, intent(in) :: side !! apply the transformation from the left (l) or right (r).
integer, intent(in) :: incv !! the increment between elements of v. incv <> 0.
integer, intent(in) :: ldc !! the leading dimension of the array c. ldc >= max(1,m).
integer, intent(in) :: m !! the number of rows of the matrix c.
integer, intent(in) :: n !! the number of columns of the matrix c.
real(dp), intent(in) :: tau !! the value tau in the representation of h.
!     ..
!     .. array arguments ..
real(dp), intent(inout) :: c( ldc, * ) !! on entry, the m by n matrix c. on exit, c is overwritten
!! by the matrix h * c if side = 'l',
real(dp), intent(in) :: v( * ) !! (1 + (m-1)*abs(incv)) if side = 'l' or (1 + (n-1)*abs(incv)) if
!! side = 'r' the vector v in the representation of h. v is not
!! used if
real(dp), intent(out) :: work( * ) !! (n) if side = 'l' or (m) if side = 'r' to take advantage of
!! the fact that v(1) = 1, we do the following
!     ..
!
!  =====================================================================
!
!     .. parameters ..
real(dp) :: one, zero
parameter          ( one = 1.0d+0, zero = 0.0d+0 )
!     ..
!     .. local scalars ..
logical :: applyleft
integer :: i, lastv, lastc
!     ..
!     .. external subroutines ..
!     ..
!     .. external functions ..
!     ..
!     .. executable statements ..
!
applyleft = lsame( side, 'l' )
lastv = 1
lastc = 0
if( tau.ne.zero ) then
!     set up variables for scanning v.  lastv begins pointing to the end
!     of v.
if( applyleft ) then
lastv = m
else
lastv = n
end if
if( incv.gt.0 ) then
i = 1 + (lastv-1) * incv
else
i = 1
end if
!     look for the last non-zero row in v.
!        since we are assuming that v(1) = 1, and it is not stored, so we
!        shouldn't access it.
do while( lastv.gt.1 .and. v( i ).eq.zero )
lastv = lastv - 1
i = i - incv
end do
if( applyleft ) then
!     scan for the last non-zero column in c(1:lastv,:).
lastc = iladlc(lastv, n, c, ldc)
else
!     scan for the last non-zero row in c(:,1:lastv).
lastc = iladlr(m, lastv, c, ldc)
end if
end if
if( lastc.eq.0 ) then
return
end if
if( applyleft ) then
!
!        form  h * c
!
! check if lastv = 1. this means v = 1, so we just need to comp
! c := hc = (1-\tau)c.
if( lastv.eq.1 ) then
!
!           c(1,1:lastc) := ( 1 - tau ) * c(1,1:lastc)
!
call dscal(lastc, one - tau, c, ldc)
else
!
!           w(1:lastc,1) := c(1:lastv,1:lastc)**t * v(1:lastv,1)
!
! w(1:lastc,1) := c(2:lastv,1:lastc)**t * v(2:lastv,1)
call dgemv( 'transpose', lastv-1, lastc, one, c(1+1,1), ldc, v(1+incv), incv, zero, work, 1)
! w(1:lastc,1) += c(1,1:lastc)**t * v(1,1) = c(1,1:lastc)**t
call daxpy(lastc, one, c, ldc, work, 1)
!
!           c(1:lastv,1:lastc) := c(...) - tau * v(1:lastv,1) * w(1:lastc,1)**t
!
! c(1, 1:lastc)   := c(...) - tau * v(1,1) * w(1:lastc,1)**t
!                  = c(...) - tau * w(1:lastc,1)**t
call daxpy(lastc, -tau, work, 1, c, ldc)
! c(2:lastv,1:lastc) := c(...) - tau * v(2:lastv,1)*w(1:last
call dger(lastv-1, lastc, -tau, v(1+incv), incv, work, 1, c(1+1,1), ldc)
end if
else
!
!        form  c * h
!
! check if n = 1. this means v = 1, so we just need to compute
! c := ch = c(1-\tau).
if( lastv.eq.1 ) then
!
!           c(1:lastc,1) := ( 1 - tau ) * c(1:lastc,1)
!
call dscal(lastc, one - tau, c, 1)
else
!
!           w(1:lastc,1) := c(1:lastc,1:lastv) * v(1:lastv,1)
!
! w(1:lastc,1) := c(1:lastc,2:lastv) * v(2:lastv,1)
call dgemv( 'no transpose', lastc, lastv-1, one, c(1,1+1), ldc, v(1+incv), incv, zero, work, 1 )
! w(1:lastc,1) += c(1:lastc,1) v(1,1) = c(1:lastc,1)
call daxpy(lastc, one, c, 1, work, 1)
!
!           c(1:lastc,1:lastv) := c(...) - tau * w(1:lastc,1) * v(1:lastv,1)**t
!
! c(1:lastc,1)     := c(...) - tau * w(1:lastc,1) * v(1,1)**
!                   = c(...) - tau * w(1:lastc,1)
call daxpy(lastc, -tau, work, 1, c, 1)
! c(1:lastc,2:lastv) := c(...) - tau * w(1:lastc,1) * v(2:la
call dger( lastc, lastv-1, -tau, work, 1, v(1+incv), incv, c(1,1+1), ldc )
end if
end if
return
!
!     end of dlarf1f
!
end subroutine dlarf1f
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlarf1l.f
!> \brief \b dlarf1l applies an elementary reflector to a general rectangular
!              matrix assuming v(lastv) = 1 where lastv is the last non-zero
!              element
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download dlarf1l + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarf1l.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarf1l.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarf1l.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       subroutine dlarf1l( side, m, n, v, incv, tau, c, ldc, work )
!
!       .. scalar arguments ..
!       character          side
!       integer            incv, ldc, m, n
!       double precision   tau
!       ..
!       .. array arguments ..
!       double precision   c( ldc, * ), v( * ), work( * )
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dlarf1l applies a real elementary reflector h to a real m by n matrix
!> c, from either the left or the right. h is represented in the form
!>
!>       h = i - tau * v * v**t
!>
!> where tau is a real scalar and v is a real vector.
!>
!> if tau = 0, then h is taken to be the unit matrix.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] side
!> \verbatim
!>          side is character*1
!>          = 'l': form  h * c
!>          = 'r': form  c * h
!> \endverbatim
!>
!> \param[in] m
!> \verbatim
!>          m is integer
!>          the number of rows of the matrix c.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is integer
!>          the number of columns of the matrix c.
!> \endverbatim
!>
!> \param[in] v
!> \verbatim
!>          v is double precision array, dimension
!>                     (1 + (m-1)*abs(incv)) if side = 'l'
!>                  or (1 + (n-1)*abs(incv)) if side = 'r'
!>          the vector v in the representation of h. v is not used if
!>          tau = 0.
!> \endverbatim
!>
!> \param[in] incv
!> \verbatim
!>          incv is integer
!>          the increment between elements of v. incv <> 0.
!> \endverbatim
!>
!> \param[in] tau
!> \verbatim
!>          tau is double precision
!>          the value tau in the representation of h.
!> \endverbatim
!>
!> \param[in,out] c
!> \verbatim
!>          c is double precision array, dimension (ldc,n)
!>          on entry, the m by n matrix c.
!>          on exit, c is overwritten by the matrix h * c if side = 'l',
!>          or c * h if side = 'r'.
!> \endverbatim
!>
!> \param[in] ldc
!> \verbatim
!>          ldc is integer
!>          the leading dimension of the array c. ldc >= max(1,m).
!> \endverbatim
!>
!> \param[out] work
!> \verbatim
!>          work is double precision array, dimension
!>                         (n) if side = 'l'
!>                      or (m) if side = 'r'
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup larf
!
!  =====================================================================
subroutine dlarf1l( side, m, n, v, incv, tau, c, ldc, work )
implicit none
!
!  -- lapack auxiliary routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
character, intent(in) :: side !! apply the transformation from the left (l) or right (r).
integer, intent(in) :: incv !! the increment between elements of v. incv <> 0.
integer, intent(in) :: ldc !! the leading dimension of the array c. ldc >= max(1,m).
integer, intent(in) :: m !! the number of rows of the matrix c.
integer, intent(in) :: n !! the number of columns of the matrix c.
real(dp), intent(in) :: tau !! the value tau in the representation of h.
!     ..
!     .. array arguments ..
real(dp), intent(inout) :: c( ldc, * ) !! on entry, the m by n matrix c. on exit, c is overwritten
!! by the matrix h * c if side = 'l',
real(dp), intent(in) :: v( * ) !! (1 + (m-1)*abs(incv)) if side = 'l' or (1 + (n-1)*abs(incv)) if
!! side = 'r' the vector v in the representation of h. v is not
!! used if
real(dp), intent(out) :: work( * ) !! (n) if side = 'l' or (m) if side = 'r'
!     ..
!
!  =====================================================================
!
!     .. parameters ..
real(dp) :: one, zero
parameter          ( one = 1.0d+0, zero = 0.0d+0 )
!     ..
!     .. local scalars ..
logical :: applyleft
integer :: i, firstv, lastv, lastc
!     ..
!     .. external subroutines ..
!     ..
!     .. external functions ..
!     ..
!     .. executable statements ..
!
applyleft = lsame( side, 'l' )
firstv = 1
lastc = 0
if( tau.ne.zero ) then
!     set up variables for scanning v.  lastv begins pointing to the end
!     of v.
if( applyleft ) then
lastv = m
else
lastv = n
end if
i = 1
!     look for the last non-zero row in v.
do while( lastv.gt.firstv .and. v( i ).eq.zero )
firstv = firstv + 1
i = i + incv
end do
if( applyleft ) then
!     scan for the last non-zero column in c(1:lastv,:).
lastc = iladlc(lastv, n, c, ldc)
else
!     scan for the last non-zero row in c(:,1:lastv).
lastc = iladlr(m, lastv, c, ldc)
end if
end if
if( lastc.eq.0 ) then
return
end if
if( applyleft ) then
!
!        form  h * c
!
if( lastv.gt.0 ) then
! check if m = 1. this means v = 1, so we just need to compu
! c := hc = (1-\tau)c.
if( lastv.eq.firstv ) then
call dscal(lastc, one - tau, c( firstv, 1), ldc)
else
!
!              w(1:lastc,1) := c(1:lastv,1:lastc)**t * v(1:lastv,1)
!
! w(1:lastc,1) := c(1:lastv-1,1:lastc)**t * v(1:lastv-1,1
call dgemv( 'transpose', lastv-firstv, lastc, one, c(firstv,1), ldc, v(i), incv, zero, work, 1)
! w(1:lastc,1) += c(lastv,1:lastc)**t * v(lastv,1) = c(la
call daxpy(lastc, one, c(lastv,1), ldc, work, 1)
!
!              c(1:lastv,1:lastc) := c(...) - tau * v(1:lastv,1) * w(1:lastc,1)**t
!
! c(lastv, 1:lastc)   := c(...) - tau * v(lastv,1) * w(1:
!                      = c(...) - tau * w(1:lastc,1)**t
call daxpy(lastc, -tau, work, 1, c(lastv,1), ldc)
! c(1:lastv-1,1:lastc) := c(...) - tau * v(1:lastv-1,1)*w
call dger(lastv-firstv, lastc, -tau, v(i), incv, work, 1, c(firstv,1), ldc)
end if
end if
else
!
!        form  c * h
!
if( lastv.gt.0 ) then
! check if n = 1. this means v = 1, so we just need to compu
! c := ch = c(1-\tau).
if( lastv.eq.firstv ) then
call dscal(lastc, one - tau, c, 1)
else
!
!              w(1:lastc,1) := c(1:lastc,1:lastv) * v(1:lastv,1)
!
! w(1:lastc,1) := c(1:lastc,1:lastv-1) * v(1:lastv-1,1)
call dgemv( 'no transpose', lastc, lastv-firstv, one, c(1,firstv), ldc, v(i), incv, zero, work, 1 )
! w(1:lastc,1) += c(1:lastc,lastv) * v(lastv,1) = c(1:las
call daxpy(lastc, one, c(1,lastv), 1, work, 1)
!
!              c(1:lastc,1:lastv) := c(...) - tau * w(1:lastc,1) * v(1:lastv,1)**t
!
! c(1:lastc,lastv)     := c(...) - tau * w(1:lastc,1) * v
!                       = c(...) - tau * w(1:lastc,1)
call daxpy(lastc, -tau, work, 1, c(1,lastv), 1)
! c(1:lastc,1:lastv-1) := c(...) - tau * w(1:lastc,1) * v
call dger( lastc, lastv-firstv, -tau, work, 1, v(i), incv, c(1,firstv), ldc )
end if
end if
end if
return
!
!     end of dlarf1l
!
end subroutine dlarf1l
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlarf.f
!> \brief \b dlarf applies an elementary reflector to a general rectangular matrix.
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download dlarf + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarf.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarf.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarf.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       subroutine dlarf( side, m, n, v, incv, tau, c, ldc, work )
!
!       .. scalar arguments ..
!       character          side
!       integer            incv, ldc, m, n
!       double precision   tau
!       ..
!       .. array arguments ..
!       double precision   c( ldc, * ), v( * ), work( * )
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dlarf applies a real elementary reflector h to a real m by n matrix
!> c, from either the left or the right. h is represented in the form
!>
!>       h = i - tau * v * v**t
!>
!> where tau is a real scalar and v is a real vector.
!>
!> if tau = 0, then h is taken to be the unit matrix.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] side
!> \verbatim
!>          side is character*1
!>          = 'l': form  h * c
!>          = 'r': form  c * h
!> \endverbatim
!>
!> \param[in] m
!> \verbatim
!>          m is integer
!>          the number of rows of the matrix c.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is integer
!>          the number of columns of the matrix c.
!> \endverbatim
!>
!> \param[in] v
!> \verbatim
!>          v is double precision array, dimension
!>                     (1 + (m-1)*abs(incv)) if side = 'l'
!>                  or (1 + (n-1)*abs(incv)) if side = 'r'
!>          the vector v in the representation of h. v is not used if
!>          tau = 0.
!> \endverbatim
!>
!> \param[in] incv
!> \verbatim
!>          incv is integer
!>          the increment between elements of v. incv <> 0.
!> \endverbatim
!>
!> \param[in] tau
!> \verbatim
!>          tau is double precision
!>          the value tau in the representation of h.
!> \endverbatim
!>
!> \param[in,out] c
!> \verbatim
!>          c is double precision array, dimension (ldc,n)
!>          on entry, the m by n matrix c.
!>          on exit, c is overwritten by the matrix h * c if side = 'l',
!>          or c * h if side = 'r'.
!> \endverbatim
!>
!> \param[in] ldc
!> \verbatim
!>          ldc is integer
!>          the leading dimension of the array c. ldc >= max(1,m).
!> \endverbatim
!>
!> \param[out] work
!> \verbatim
!>          work is double precision array, dimension
!>                         (n) if side = 'l'
!>                      or (m) if side = 'r'
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup larf
!
!  =====================================================================
subroutine dlarf( side, m, n, v, incv, tau, c, ldc, work )
implicit none
!
!  -- lapack auxiliary routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
character, intent(in) :: side !! apply the transformation from the left (l) or right (r).
integer, intent(in) :: incv !! the increment between elements of v. incv <> 0.
integer, intent(in) :: ldc !! the leading dimension of the array c. ldc >= max(1,m).
integer, intent(in) :: m !! the number of rows of the matrix c.
integer, intent(in) :: n !! the number of columns of the matrix c.
real(dp), intent(in) :: tau !! the value tau in the representation of h.
!     ..
!     .. array arguments ..
real(dp), intent(inout) :: c( ldc, * ) !! on entry, the m by n matrix c. on exit, c is overwritten
!! by the matrix h * c if side = 'l',
real(dp), intent(in) :: v( * ) !! (1 + (m-1)*abs(incv)) if side = 'l' or (1 + (n-1)*abs(incv)) if
!! side = 'r' the vector v in the representation of h. v is not
!! used if
real(dp), intent(out) :: work( * ) !! (n) if side = 'l' or (m) if side = 'r'
!     ..
!
!  =====================================================================
!
!     .. parameters ..
real(dp) :: one, zero
parameter          ( one = 1.0d+0, zero = 0.0d+0 )
!     ..
!     .. local scalars ..
logical :: applyleft
integer :: i, lastv, lastc
!     ..
!     .. external subroutines ..
!     ..
!     .. external functions ..
!     ..
!     .. executable statements ..
!
applyleft = lsame( side, 'l' )
lastv = 0
lastc = 0
if( tau.ne.zero ) then
!     set up variables for scanning v.  lastv begins pointing to the end
!     of v.
if( applyleft ) then
lastv = m
else
lastv = n
end if
if( incv.gt.0 ) then
i = 1 + (lastv-1) * incv
else
i = 1
end if
!     look for the last non-zero row in v.
do while( lastv.gt.0 .and. v( i ).eq.zero )
lastv = lastv - 1
i = i - incv
end do
if( applyleft ) then
!     scan for the last non-zero column in c(1:lastv,:).
lastc = iladlc(lastv, n, c, ldc)
else
!     scan for the last non-zero row in c(:,1:lastv).
lastc = iladlr(m, lastv, c, ldc)
end if
end if
!     note that lastc.eq.0 renders the blas operations null; no special
!     case is needed at this level.
if( applyleft ) then
!
!        form  h * c
!
if( lastv.gt.0 ) then
!
!           w(1:lastc,1) := c(1:lastv,1:lastc)**t * v(1:lastv,1)
!
call dgemv( 'transpose', lastv, lastc, one, c, ldc, v, incv, zero, work, 1 )
!
!           c(1:lastv,1:lastc) := c(...) - v(1:lastv,1) * w(1:lastc,1)**t
!
call dger( lastv, lastc, -tau, v, incv, work, 1, c, ldc )
end if
else
!
!        form  c * h
!
if( lastv.gt.0 ) then
!
!           w(1:lastc,1) := c(1:lastc,1:lastv) * v(1:lastv,1)
!
call dgemv( 'no transpose', lastc, lastv, one, c, ldc, v, incv, zero, work, 1 )
!
!           c(1:lastc,1:lastv) := c(...) - w(1:lastc,1) * v(1:lastv,1)**t
!
call dger( lastc, lastv, -tau, work, 1, v, incv, c, ldc )
end if
end if
return
!
!     end of dlarf
!
end subroutine dlarf
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/BLAS/SRC/daxpy.f
!> \brief \b daxpy
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  definition:
!  ===========
!
!       subroutine daxpy(n,da,dx,incx,dy,incy)
!
!       .. scalar arguments ..
!       double precision da
!       integer incx,incy,n
!       ..
!       .. array arguments ..
!       double precision dx(*),dy(*)
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!>    daxpy constant times a vector plus a vector.
!>    uses unrolled loops for increments equal to one.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] n
!> \verbatim
!>          n is integer
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in] da
!> \verbatim
!>          da is double precision
!>           on entry, da specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] dx
!> \verbatim
!>          dx is double precision array, dimension ( 1 + ( n - 1 )*abs( incx ) )
!> \endverbatim
!>
!> \param[in] incx
!> \verbatim
!>          incx is integer
!>         storage spacing between elements of dx
!> \endverbatim
!>
!> \param[in,out] dy
!> \verbatim
!>          dy is double precision array, dimension ( 1 + ( n - 1 )*abs( incy ) )
!> \endverbatim
!>
!> \param[in] incy
!> \verbatim
!>          incy is integer
!>         storage spacing between elements of dy
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup axpy
!
!> \par further details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
subroutine daxpy(n,da,dx,incx,dy,incy)
implicit none
!
!  -- reference blas level1 routine --
!  -- reference blas is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
real(dp), intent(in) :: da !! on entry, da specifies the scalar alpha.
integer, intent(in) :: incx !! storage spacing between elements of dx
integer, intent(in) :: incy !! storage spacing between elements of dy
integer, intent(in) :: n !! number of elements in input vector(s)
!     ..
!     .. array arguments ..
real(dp), intent(in) :: dx(*) !! first strided vector, updated in place when requested by the operation.
real(dp), intent(inout) :: dy(*) !! second strided vector, updated in place when requested by the operation.
!     ..
!
!  =====================================================================
!
!     .. local scalars ..
integer :: i, ix, iy, m, mp1
!     ..
!     .. intrinsic functions ..
intrinsic mod
!     ..
if (n.le.0) return
if (da.eq.0.0d0) return
if (incx.eq.1 .and. incy.eq.1) then
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
m = mod(n,4)
if (m.ne.0) then
do i = 1,m
dy(i) = dy(i) + da*dx(i)
end do
end if
if (n.lt.4) return
mp1 = m + 1
do i = mp1,n,4
dy(i) = dy(i) + da*dx(i)
dy(i+1) = dy(i+1) + da*dx(i+1)
dy(i+2) = dy(i+2) + da*dx(i+2)
dy(i+3) = dy(i+3) + da*dx(i+3)
end do
else
!
!        code for unequal increments or equal increments
!          not equal to 1
!
ix = 1
iy = 1
if (incx.lt.0) ix = (-n+1)*incx + 1
if (incy.lt.0) iy = (-n+1)*incy + 1
do i = 1,n
dy(iy) = dy(iy) + da*dx(ix)
ix = ix + incx
iy = iy + incy
end do
end if
return
!
!     end of daxpy
!
end subroutine daxpy
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/BLAS/SRC/dcopy.f
!> \brief \b dcopy
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  definition:
!  ===========
!
!       subroutine dcopy(n,dx,incx,dy,incy)
!
!       .. scalar arguments ..
!       integer incx,incy,n
!       ..
!       .. array arguments ..
!       double precision dx(*),dy(*)
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!>    dcopy copies a vector, x, to a vector, y.
!>    uses unrolled loops for increments equal to 1.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] n
!> \verbatim
!>          n is integer
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in] dx
!> \verbatim
!>          dx is double precision array, dimension ( 1 + ( n - 1 )*abs( incx ) )
!> \endverbatim
!>
!> \param[in] incx
!> \verbatim
!>          incx is integer
!>         storage spacing between elements of dx
!> \endverbatim
!>
!> \param[out] dy
!> \verbatim
!>          dy is double precision array, dimension ( 1 + ( n - 1 )*abs( incy ) )
!> \endverbatim
!>
!> \param[in] incy
!> \verbatim
!>          incy is integer
!>         storage spacing between elements of dy
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup copy
!
!> \par further details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
subroutine dcopy(n,dx,incx,dy,incy)
implicit none
!
!  -- reference blas level1 routine --
!  -- reference blas is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
integer, intent(in) :: incx !! storage spacing between elements of dx
integer, intent(in) :: incy !! storage spacing between elements of dy
integer, intent(in) :: n !! number of elements in input vector(s)
!     ..
!     .. array arguments ..
real(dp), intent(in) :: dx(*) !! first strided vector, updated in place when requested by the operation.
real(dp), intent(out) :: dy(*) !! second strided vector, updated in place when requested by the operation.
!     ..
!
!  =====================================================================
!
!     .. local scalars ..
integer :: i, ix, iy, m, mp1
!     ..
!     .. intrinsic functions ..
intrinsic mod
!     ..
if (n.le.0) return
if (incx.eq.1 .and. incy.eq.1) then
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
m = mod(n,7)
if (m.ne.0) then
do i = 1,m
dy(i) = dx(i)
end do
if (n.lt.7) return
end if
mp1 = m + 1
do i = mp1,n,7
dy(i) = dx(i)
dy(i+1) = dx(i+1)
dy(i+2) = dx(i+2)
dy(i+3) = dx(i+3)
dy(i+4) = dx(i+4)
dy(i+5) = dx(i+5)
dy(i+6) = dx(i+6)
end do
else
!
!        code for unequal increments or equal increments
!          not equal to 1
!
ix = 1
iy = 1
if (incx.lt.0) ix = (-n+1)*incx + 1
if (incy.lt.0) iy = (-n+1)*incy + 1
do i = 1,n
dy(iy) = dx(ix)
ix = ix + incx
iy = iy + incy
end do
end if
return
!
!     end of dcopy
!
end subroutine dcopy
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/BLAS/SRC/dgemv.f
!> \brief \b dgemv
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  definition:
!  ===========
!
!       subroutine dgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
!
!       .. scalar arguments ..
!       double precision alpha,beta
!       integer incx,incy,lda,m,n
!       character trans
!       ..
!       .. array arguments ..
!       double precision a(lda,*),x(*),y(*)
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dgemv  performs one of the matrix-vector operations
!>
!>    y := alpha*a*x + beta*y,   or   y := alpha*a**t*x + beta*y,
!>
!> where alpha and beta are scalars, x and y are vectors and a is an
!> m by n matrix.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] trans
!> \verbatim
!>          trans is character*1
!>           on entry, trans specifies the operation to be performed as
!>           follows:
!>
!>              trans = 'n' or 'n'   y := alpha*a*x + beta*y.
!>
!>              trans = 't' or 't'   y := alpha*a**t*x + beta*y.
!>
!>              trans = 'c' or 'c'   y := alpha*a**t*x + beta*y.
!> \endverbatim
!>
!> \param[in] m
!> \verbatim
!>          m is integer
!>           on entry, m specifies the number of rows of the matrix a.
!>           m must be at least zero.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is integer
!>           on entry, n specifies the number of columns of the matrix a.
!>           n must be at least zero.
!> \endverbatim
!>
!> \param[in] alpha
!> \verbatim
!>          alpha is double precision.
!>           on entry, alpha specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] a
!> \verbatim
!>          a is double precision array, dimension ( lda, n )
!>           before entry, the leading m by n part of the array a must
!>           contain the matrix of coefficients.
!> \endverbatim
!>
!> \param[in] lda
!> \verbatim
!>          lda is integer
!>           on entry, lda specifies the first dimension of a as declared
!>           in the calling (sub) program. lda must be at least
!>           max( 1, m ).
!> \endverbatim
!>
!> \param[in] x
!> \verbatim
!>          x is double precision array, dimension at least
!>           ( 1 + ( n - 1 )*abs( incx ) ) when trans = 'n' or 'n'
!>           and at least
!>           ( 1 + ( m - 1 )*abs( incx ) ) otherwise.
!>           before entry, the incremented array x must contain the
!>           vector x.
!> \endverbatim
!>
!> \param[in] incx
!> \verbatim
!>          incx is integer
!>           on entry, incx specifies the increment for the elements of
!>           x. incx must not be zero.
!> \endverbatim
!>
!> \param[in] beta
!> \verbatim
!>          beta is double precision.
!>           on entry, beta specifies the scalar beta. when beta is
!>           supplied as zero then y need not be set on input.
!> \endverbatim
!>
!> \param[in,out] y
!> \verbatim
!>          y is double precision array, dimension at least
!>           ( 1 + ( m - 1 )*abs( incy ) ) when trans = 'n' or 'n'
!>           and at least
!>           ( 1 + ( n - 1 )*abs( incy ) ) otherwise.
!>           before entry with beta non-zero, the incremented array y
!>           must contain the vector y. on exit, y is overwritten by the
!>           updated vector y.
!>           if either m or n is zero, then y not referenced and the function
!>           performs a quick return.
!> \endverbatim
!>
!> \param[in] incy
!> \verbatim
!>          incy is integer
!>           on entry, incy specifies the increment for the elements of
!>           y. incy must not be zero.
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup gemv
!
!> \par further details:
!  =====================
!>
!> \verbatim
!>
!>  level 2 blas routine.
!>  the vector and matrix arguments are not referenced when n = 0, or m = 0
!>
!>  -- written on 22-october-1986.
!>     jack dongarra, argonne national lab.
!>     jeremy du croz, nag central office.
!>     sven hammarling, nag central office.
!>     richard hanson, sandia national labs.
!> \endverbatim
!>
!  =====================================================================
subroutine dgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
implicit none
!
!  -- reference blas level2 routine --
!  -- reference blas is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
real(dp), intent(in) :: alpha !! on entry, alpha specifies the scalar alpha.
real(dp), intent(in) :: beta !! on entry, beta specifies the scalar beta. when beta is supplied as
!! zero then y need not be set on input.
integer, intent(in) :: incx !! on entry, incx specifies the increment for the elements of x. incx
!! must not be zero.
integer, intent(in) :: incy !! on entry, incy specifies the increment for the elements of y. incy
!! must not be zero.
integer, intent(in) :: lda !! on entry, lda specifies the first dimension of a as declared in the
!! calling (sub) program. lda must be at least
integer, intent(in) :: m !! on entry, m specifies the number of rows of the matrix a. m must be at
!! least zero.
integer, intent(in) :: n !! on entry, n specifies the number of columns of the matrix a. n must be
!! at least zero.
character, intent(in) :: trans !! on entry, trans specifies the operation to be performed as
!! follows: trans = 'n' or 'n'   y := alpha*a*x + beta*y.
!     ..
!     .. array arguments ..
real(dp), intent(in) :: a(lda,*) !! before entry, the leading m by n part of the array a must
!! contain the matrix of coefficients.
real(dp), intent(in) :: x(*) !! ( 1 + ( n - 1 )*abs( incx ) ) when trans = 'n' or 'n' and at least
!! ( 1 + ( m - 1 )*abs( incx ) ) otherwise.
real(dp), intent(inout) :: y(*) !! ( 1 + ( m - 1 )*abs( incy ) ) when trans = 'n' or 'n' and at
!! least ( 1 + ( n - 1 )*abs( incy ) ) otherwise.
!     ..
!
!  =====================================================================
!
!     .. parameters ..
real(dp) :: one, zero
parameter (one=1.0d+0,zero=0.0d+0)
!     ..
!     .. local scalars ..
real(dp) :: temp
integer :: i, info, ix, iy, j, jx, jy, kx, ky, lenx, leny
!     ..
!     .. external functions ..
!     ..
!     .. external subroutines ..
!     ..
!     .. intrinsic functions ..
intrinsic max
!     ..
!
!     test the input parameters.
!
info = 0
if (.not.lsame(trans,'n') .and. .not.lsame(trans,'t') .and. .not.lsame(trans,'c')) then
info = 1
else if (m.lt.0) then
info = 2
else if (n.lt.0) then
info = 3
else if (lda.lt.max(1,m)) then
info = 6
else if (incx.eq.0) then
info = 8
else if (incy.eq.0) then
info = 11
end if
if (info.ne.0) then
call xerbla('dgemv ',info)
return
end if
!
!     quick return if possible.
!
if ((m.eq.0) .or. (n.eq.0) .or. ((alpha.eq.zero).and. (beta.eq.one))) return
!
!     set  lenx  and  leny, the lengths of the vectors x and y, and set
!     up the start points in  x  and  y.
!
if (lsame(trans,'n')) then
lenx = n
leny = m
else
lenx = m
leny = n
end if
if (incx.gt.0) then
kx = 1
else
kx = 1 - (lenx-1)*incx
end if
if (incy.gt.0) then
ky = 1
else
ky = 1 - (leny-1)*incy
end if
!
!     start the operations. in this version the elements of a are
!     accessed sequentially with one pass through a.
!
!     first form  y := beta*y.
!
if (beta.ne.one) then
if (incy.eq.1) then
if (beta.eq.zero) then
do 10 i = 1,leny
y(i) = zero
10 continue
else
do 20 i = 1,leny
y(i) = beta*y(i)
20 continue
end if
else
iy = ky
if (beta.eq.zero) then
do 30 i = 1,leny
y(iy) = zero
iy = iy + incy
30 continue
else
do 40 i = 1,leny
y(iy) = beta*y(iy)
iy = iy + incy
40 continue
end if
end if
end if
if (alpha.eq.zero) return
if (lsame(trans,'n')) then
!
!        form  y := alpha*a*x + y.
!
jx = kx
if (incy.eq.1) then
do 60 j = 1,n
temp = alpha*x(jx)
do 50 i = 1,m
y(i) = y(i) + temp*a(i,j)
50 continue
jx = jx + incx
60 continue
else
do 80 j = 1,n
temp = alpha*x(jx)
iy = ky
do 70 i = 1,m
y(iy) = y(iy) + temp*a(i,j)
iy = iy + incy
70 continue
jx = jx + incx
80 continue
end if
else
!
!        form  y := alpha*a**t*x + y.
!
jy = ky
if (incx.eq.1) then
do 100 j = 1,n
temp = zero
do 90 i = 1,m
temp = temp + a(i,j)*x(i)
90 continue
y(jy) = y(jy) + alpha*temp
jy = jy + incy
100 continue
else
do 120 j = 1,n
temp = zero
ix = kx
do 110 i = 1,m
temp = temp + a(i,j)*x(ix)
ix = ix + incx
110 continue
y(jy) = y(jy) + alpha*temp
jy = jy + incy
120 continue
end if
end if
!
return
!
!     end of dgemv
!
end subroutine dgemv
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/BLAS/SRC/dger.f
!> \brief \b dger
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  definition:
!  ===========
!
!       subroutine dger(m,n,alpha,x,incx,y,incy,a,lda)
!
!       .. scalar arguments ..
!       double precision alpha
!       integer incx,incy,lda,m,n
!       ..
!       .. array arguments ..
!       double precision a(lda,*),x(*),y(*)
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dger   performs the rank 1 operation
!>
!>    a := alpha*x*y**t + a,
!>
!> where alpha is a scalar, x is an m element vector, y is an n element
!> vector and a is an m by n matrix.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] m
!> \verbatim
!>          m is integer
!>           on entry, m specifies the number of rows of the matrix a.
!>           m must be at least zero.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is integer
!>           on entry, n specifies the number of columns of the matrix a.
!>           n must be at least zero.
!> \endverbatim
!>
!> \param[in] alpha
!> \verbatim
!>          alpha is double precision.
!>           on entry, alpha specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] x
!> \verbatim
!>          x is double precision array, dimension at least
!>           ( 1 + ( m - 1 )*abs( incx ) ).
!>           before entry, the incremented array x must contain the m
!>           element vector x.
!> \endverbatim
!>
!> \param[in] incx
!> \verbatim
!>          incx is integer
!>           on entry, incx specifies the increment for the elements of
!>           x. incx must not be zero.
!> \endverbatim
!>
!> \param[in] y
!> \verbatim
!>          y is double precision array, dimension at least
!>           ( 1 + ( n - 1 )*abs( incy ) ).
!>           before entry, the incremented array y must contain the n
!>           element vector y.
!> \endverbatim
!>
!> \param[in] incy
!> \verbatim
!>          incy is integer
!>           on entry, incy specifies the increment for the elements of
!>           y. incy must not be zero.
!> \endverbatim
!>
!> \param[in,out] a
!> \verbatim
!>          a is double precision array, dimension ( lda, n )
!>           before entry, the leading m by n part of the array a must
!>           contain the matrix of coefficients. on exit, a is
!>           overwritten by the updated matrix.
!> \endverbatim
!>
!> \param[in] lda
!> \verbatim
!>          lda is integer
!>           on entry, lda specifies the first dimension of a as declared
!>           in the calling (sub) program. lda must be at least
!>           max( 1, m ).
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup ger
!
!> \par further details:
!  =====================
!>
!> \verbatim
!>
!>  level 2 blas routine.
!>
!>  -- written on 22-october-1986.
!>     jack dongarra, argonne national lab.
!>     jeremy du croz, nag central office.
!>     sven hammarling, nag central office.
!>     richard hanson, sandia national labs.
!> \endverbatim
!>
!  =====================================================================
subroutine dger(m,n,alpha,x,incx,y,incy,a,lda)
implicit none
!
!  -- reference blas level2 routine --
!  -- reference blas is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
real(dp), intent(in) :: alpha !! on entry, alpha specifies the scalar alpha.
integer, intent(in) :: incx !! on entry, incx specifies the increment for the elements of x. incx
!! must not be zero.
integer, intent(in) :: incy !! on entry, incy specifies the increment for the elements of y. incy
!! must not be zero.
integer, intent(in) :: lda !! on entry, lda specifies the first dimension of a as declared in the
!! calling (sub) program. lda must be at least
integer, intent(in) :: m !! on entry, m specifies the number of rows of the matrix a. m must be at
!! least zero.
integer, intent(in) :: n !! on entry, n specifies the number of columns of the matrix a. n must be
!! at least zero.
!     ..
!     .. array arguments ..
real(dp), intent(inout) :: a(lda,*) !! before entry, the leading m by n part of the array a must
!! contain the matrix of coefficients. on exit, a is
real(dp), intent(in) :: x(*) !! ( 1 + ( m - 1 )*abs( incx ) ). before entry, the incremented array
!! x must contain the m element vector x.
real(dp), intent(in) :: y(*) !! ( 1 + ( n - 1 )*abs( incy ) ). before entry, the incremented array
!! y must contain the n element vector y.
!     ..
!
!  =====================================================================
!
!     .. parameters ..
real(dp) :: zero
parameter (zero=0.0d+0)
!     ..
!     .. local scalars ..
real(dp) :: temp
integer :: i, info, ix, j, jy, kx
!     ..
!     .. external subroutines ..
!     ..
!     .. intrinsic functions ..
intrinsic max
!     ..
!
!     test the input parameters.
!
info = 0
if (m.lt.0) then
info = 1
else if (n.lt.0) then
info = 2
else if (incx.eq.0) then
info = 5
else if (incy.eq.0) then
info = 7
else if (lda.lt.max(1,m)) then
info = 9
end if
if (info.ne.0) then
call xerbla('dger  ',info)
return
end if
!
!     quick return if possible.
!
if ((m.eq.0) .or. (n.eq.0) .or. (alpha.eq.zero)) return
!
!     start the operations. in this version the elements of a are
!     accessed sequentially with one pass through a.
!
if (incy.gt.0) then
jy = 1
else
jy = 1 - (n-1)*incy
end if
if (incx.eq.1) then
do 20 j = 1,n
if (y(jy).ne.zero) then
temp = alpha*y(jy)
do 10 i = 1,m
a(i,j) = a(i,j) + x(i)*temp
10 continue
end if
jy = jy + incy
20 continue
else
if (incx.gt.0) then
kx = 1
else
kx = 1 - (m-1)*incx
end if
do 40 j = 1,n
if (y(jy).ne.zero) then
temp = alpha*y(jy)
ix = kx
do 30 i = 1,m
a(i,j) = a(i,j) + x(ix)*temp
ix = ix + incx
30 continue
end if
jy = jy + incy
40 continue
end if
!
return
!
!     end of dger
!
end subroutine dger
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dgesc2.f
!> \brief \b dgesc2 solves a system of linear equations using the lu factorization with complete
!pivoting computed by sgetc2.
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download dgesc2 + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgesc2.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgesc2.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgesc2.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       subroutine dgesc2( n, a, lda, rhs, ipiv, jpiv, scale )
!
!       .. scalar arguments ..
!       integer            lda, n
!       double precision   scale
!       ..
!       .. array arguments ..
!       integer            ipiv( * ), jpiv( * )
!       double precision   a( lda, * ), rhs( * )
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dgesc2 solves a system of linear equations
!>
!>           a * x = scale* rhs
!>
!> with a general n-by-n matrix a using the lu factorization with
!> complete pivoting computed by dgetc2.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] n
!> \verbatim
!>          n is integer
!>          the order of the matrix a.
!> \endverbatim
!>
!> \param[in] a
!> \verbatim
!>          a is double precision array, dimension (lda,n)
!>          on entry, the  lu part of the factorization of the n-by-n
!>          matrix a computed by dgetc2:  a = p * l * u * q
!> \endverbatim
!>
!> \param[in] lda
!> \verbatim
!>          lda is integer
!>          the leading dimension of the array a.  lda >= max(1, n).
!> \endverbatim
!>
!> \param[in,out] rhs
!> \verbatim
!>          rhs is double precision array, dimension (n).
!>          on entry, the right hand side vector b.
!>          on exit, the solution vector x.
!> \endverbatim
!>
!> \param[in] ipiv
!> \verbatim
!>          ipiv is integer array, dimension (n).
!>          the pivot indices; for 1 <= i <= n, row i of the
!>          matrix has been interchanged with row ipiv(i).
!> \endverbatim
!>
!> \param[in] jpiv
!> \verbatim
!>          jpiv is integer array, dimension (n).
!>          the pivot indices; for 1 <= j <= n, column j of the
!>          matrix has been interchanged with column jpiv(j).
!> \endverbatim
!>
!> \param[out] scale
!> \verbatim
!>          scale is double precision
!>          on exit, scale contains the scale factor. scale is chosen
!>          0 <= scale <= 1 to prevent overflow in the solution.
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup gesc2
!
!> \par contributors:
!  ==================
!>
!>     bo kagstrom and peter poromaa, department of computing science,
!>     umea university, s-901 87 umea, sweden.
!
!  =====================================================================
subroutine dgesc2( n, a, lda, rhs, ipiv, jpiv, scale )
implicit none
!
!  -- lapack auxiliary routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
integer, intent(in) :: lda !! the leading dimension of the array a.  lda >= max(1, n).
integer, intent(in) :: n !! the order of the matrix a.
real(dp), intent(out) :: scale !! on exit, scale contains the scale factor. scale is chosen 0 <=
!! scale <= 1 to prevent overflow in the solution.
!     ..
!     .. array arguments ..
integer, intent(in) :: ipiv( * ) !! the pivot indices; for 1 <= i <= n, row i of the matrix has
!! been interchanged with row ipiv(i).
integer, intent(in) :: jpiv( * ) !! the pivot indices; for 1 <= j <= n, column j of the matrix has
!! been interchanged with column jpiv(j).
real(dp), intent(in) :: a( lda, * ) !! on entry, the  lu part of the factorization of the n-by-n
!! matrix a computed by dgetc2:  a = p * l * u * q
real(dp), intent(inout) :: rhs( * ) !! on entry, the right hand side vector b. on exit, the
!! solution vector x.
!     ..
!
!  =====================================================================
!
!     .. parameters ..
real(dp) :: one, two
parameter          ( one = 1.0d+0, two = 2.0d+0 )
!     ..
!     .. local scalars ..
integer :: i, j
real(dp) :: bignum, eps, smlnum, temp
!     ..
!     .. external subroutines ..
!     ..
!     .. external functions ..
!     ..
!     .. intrinsic functions ..
intrinsic          abs
!     ..
!     .. executable statements ..
!
!      set constant to control overflow
!
eps = dlamch( 'p' )
smlnum = dlamch( 's' ) / eps
bignum = one / smlnum
!
!     apply permutations ipiv to rhs
!
call dlaswp( 1, rhs, lda, 1, n-1, ipiv, 1 )
!
!     solve for l part
!
do 20 i = 1, n - 1
do 10 j = i + 1, n
rhs( j ) = rhs( j ) - a( j, i )*rhs( i )
10 continue
20 continue
!
!     solve for u part
!
scale = one
!
!     check for scaling
!
i = idamax( n, rhs, 1 )
if( two*smlnum*abs( rhs( i ) ).gt.abs( a( n, n ) ) ) then
temp = ( one / two ) / abs( rhs( i ) )
call dscal( n, temp, rhs( 1 ), 1 )
scale = scale*temp
end if
!
do 40 i = n, 1, -1
temp = one / a( i, i )
rhs( i ) = rhs( i )*temp
do 30 j = i + 1, n
rhs( i ) = rhs( i ) - rhs( j )*( a( i, j )*temp )
30 continue
40 continue
!
!     apply permutations jpiv to the solution (rhs)
!
call dlaswp( 1, rhs, lda, 1, n-1, jpiv, -1 )
return
!
!     end of dgesc2
!
end subroutine dgesc2
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dgetc2.f
!> \brief \b dgetc2 computes the lu factorization with complete pivoting of the general n-by-n
!matrix.
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download dgetc2 + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgetc2.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgetc2.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgetc2.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       subroutine dgetc2( n, a, lda, ipiv, jpiv, info )
!
!       .. scalar arguments ..
!       integer            info, lda, n
!       ..
!       .. array arguments ..
!       integer            ipiv( * ), jpiv( * )
!       double precision   a( lda, * )
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dgetc2 computes an lu factorization with complete pivoting of the
!> n-by-n matrix a. the factorization has the form a = p * l * u * q,
!> where p and q are permutation matrices, l is lower triangular with
!> unit diagonal elements and u is upper triangular.
!>
!> this is the level 2 blas algorithm.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] n
!> \verbatim
!>          n is integer
!>          the order of the matrix a. n >= 0.
!> \endverbatim
!>
!> \param[in,out] a
!> \verbatim
!>          a is double precision array, dimension (lda, n)
!>          on entry, the n-by-n matrix a to be factored.
!>          on exit, the factors l and u from the factorization
!>          a = p*l*u*q; the unit diagonal elements of l are not stored.
!>          if u(k, k) appears to be less than smin, u(k, k) is given the
!>          value of smin, i.e., giving a nonsingular perturbed system.
!> \endverbatim
!>
!> \param[in] lda
!> \verbatim
!>          lda is integer
!>          the leading dimension of the array a.  lda >= max(1,n).
!> \endverbatim
!>
!> \param[out] ipiv
!> \verbatim
!>          ipiv is integer array, dimension(n).
!>          the pivot indices; for 1 <= i <= n, row i of the
!>          matrix has been interchanged with row ipiv(i).
!> \endverbatim
!>
!> \param[out] jpiv
!> \verbatim
!>          jpiv is integer array, dimension(n).
!>          the pivot indices; for 1 <= j <= n, column j of the
!>          matrix has been interchanged with column jpiv(j).
!> \endverbatim
!>
!> \param[out] info
!> \verbatim
!>          info is integer
!>           = 0: successful exit
!>           > 0: if info = k, u(k, k) is likely to produce overflow if
!>                we try to solve for x in ax = b. so u is perturbed to
!>                avoid the overflow.
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup getc2
!
!> \par contributors:
!  ==================
!>
!>     bo kagstrom and peter poromaa, department of computing science,
!>     umea university, s-901 87 umea, sweden.
!
!  =====================================================================
subroutine dgetc2( n, a, lda, ipiv, jpiv, info )
implicit none
!
!  -- lapack auxiliary routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
integer, intent(out) :: info !! zero success; negative invalid argument or positive numerical failure.
integer, intent(in) :: lda !! the leading dimension of the array a.  lda >= max(1,n).
integer, intent(in) :: n !! the order of the matrix a. n >= 0.
!     ..
!     .. array arguments ..
integer, intent(out) :: ipiv( * ) !! the pivot indices; for 1 <= i <= n, row i of the matrix has
!! been interchanged with row ipiv(i).
integer, intent(out) :: jpiv( * ) !! the pivot indices; for 1 <= j <= n, column j of the matrix
!! has been interchanged with column jpiv(j).
real(dp), intent(inout) :: a( lda, * ) !! on entry, the n-by-n matrix a to be factored. on exit,
!! the factors l and u from the factorization
!     ..
!
!  =====================================================================
!
!     .. parameters ..
real(dp) :: zero, one
parameter          ( zero = 0.0d+0, one = 1.0d+0 )
!     ..
!     .. local scalars ..
integer :: i, ip, ipv, j, jp, jpv
real(dp) :: bignum, eps, smin, smlnum, xmax
!     ..
!     .. external subroutines ..
!     ..
!     .. external functions ..
!     ..
!     .. intrinsic functions ..
intrinsic          abs, max
!     ..
!     .. executable statements ..
!
info = 0
!
!     quick return if possible
!
if( n.eq.0 ) return
!
!     set constants to control overflow
!
eps = dlamch( 'p' )
smlnum = dlamch( 's' ) / eps
bignum = one / smlnum
!
!     handle the case n=1 by itself
!
if( n.eq.1 ) then
ipiv( 1 ) = 1
jpiv( 1 ) = 1
if( abs( a( 1, 1 ) ).lt.smlnum ) then
info = 1
a( 1, 1 ) = smlnum
end if
return
end if
!
!     factorize a using complete pivoting.
!     set pivots less than smin to smin.
!
do 40 i = 1, n - 1
!
!        find max element in matrix a
!
xmax = zero
do 20 jp = i, n
do 10 ip = i, n
if( abs( a( ip, jp ) ).ge.xmax ) then
xmax = abs( a( ip, jp ) )
ipv = ip
jpv = jp
end if
10 continue
20 continue
if( i.eq.1 ) smin = max( eps*xmax, smlnum )
!
!        swap rows
!
if( ipv.ne.i ) call dswap( n, a( ipv, 1 ), lda, a( i, 1 ), lda )
ipiv( i ) = ipv
!
!        swap columns
!
if( jpv.ne.i ) call dswap( n, a( 1, jpv ), 1, a( 1, i ), 1 )
jpiv( i ) = jpv
!
!        check for singularity
!
if( abs( a( i, i ) ).lt.smin ) then
info = i
a( i, i ) = smin
end if
do 30 j = i + 1, n
a( j, i ) = a( j, i ) / a( i, i )
30 continue
call dger( n-i, n-i, -one, a( i+1, i ), 1, a( i, i+1 ), lda, a( i+1, i+1 ), lda )
40 continue
!
if( abs( a( n, n ) ).lt.smin ) then
info = n
a( n, n ) = smin
end if
!
!     set last pivots to n
!
ipiv( n ) = n
jpiv( n ) = n
!
return
!
!     end of dgetc2
!
end subroutine dgetc2
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlatdf.f
!> \brief \b dlatdf uses the lu factorization of the n-by-n matrix computed by sgetc2 and
!computes a contribution to the reciprocal dif-estimate.
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download dlatdf + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlatdf.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlatdf.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlatdf.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       subroutine dlatdf( ijob, n, z, ldz, rhs, rdsum, rdscal, ipiv,
!                          jpiv )
!
!       .. scalar arguments ..
!       integer            ijob, ldz, n
!       double precision   rdscal, rdsum
!       ..
!       .. array arguments ..
!       integer            ipiv( * ), jpiv( * )
!       double precision   rhs( * ), z( ldz, * )
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dlatdf uses the lu factorization of the n-by-n matrix z computed by
!> dgetc2 and computes a contribution to the reciprocal dif-estimate
!> by solving z * x = b for x, and choosing the r.h.s. b such that
!> the norm of x is as large as possible. on entry rhs = b holds the
!> contribution from earlier solved sub-systems, and on return rhs = x.
!>
!> the factorization of z returned by dgetc2 has the form z = p*l*u*q,
!> where p and q are permutation matrices. l is lower triangular with
!> unit diagonal elements and u is upper triangular.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] ijob
!> \verbatim
!>          ijob is integer
!>          ijob = 2: first compute an approximative null-vector e
!>              of z using dgecon, e is normalized and solve for
!>              zx = +-e - f with the sign giving the greater value
!>              of 2-norm(x). about 5 times as expensive as default.
!>          ijob .ne. 2: local look ahead strategy where all entries of
!>              the r.h.s. b is chosen as either +1 or -1 (default).
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is integer
!>          the number of columns of the matrix z.
!> \endverbatim
!>
!> \param[in] z
!> \verbatim
!>          z is double precision array, dimension (ldz, n)
!>          on entry, the lu part of the factorization of the n-by-n
!>          matrix z computed by dgetc2:  z = p * l * u * q
!> \endverbatim
!>
!> \param[in] ldz
!> \verbatim
!>          ldz is integer
!>          the leading dimension of the array z.  lda >= max(1, n).
!> \endverbatim
!>
!> \param[in,out] rhs
!> \verbatim
!>          rhs is double precision array, dimension (n)
!>          on entry, rhs contains contributions from other subsystems.
!>          on exit, rhs contains the solution of the subsystem with
!>          entries according to the value of ijob (see above).
!> \endverbatim
!>
!> \param[in,out] rdsum
!> \verbatim
!>          rdsum is double precision
!>          on entry, the sum of squares of computed contributions to
!>          the dif-estimate under computation by dtgsyl, where the
!>          scaling factor rdscal (see below) has been factored out.
!>          on exit, the corresponding sum of squares updated with the
!>          contributions from the current sub-system.
!>          if trans = 't' rdsum is not touched.
!>          note: rdsum only makes sense when dtgsy2 is called by stgsyl.
!> \endverbatim
!>
!> \param[in,out] rdscal
!> \verbatim
!>          rdscal is double precision
!>          on entry, scaling factor used to prevent overflow in rdsum.
!>          on exit, rdscal is updated w.r.t. the current contributions
!>          in rdsum.
!>          if trans = 't', rdscal is not touched.
!>          note: rdscal only makes sense when dtgsy2 is called by
!>                dtgsyl.
!> \endverbatim
!>
!> \param[in] ipiv
!> \verbatim
!>          ipiv is integer array, dimension (n).
!>          the pivot indices; for 1 <= i <= n, row i of the
!>          matrix has been interchanged with row ipiv(i).
!> \endverbatim
!>
!> \param[in] jpiv
!> \verbatim
!>          jpiv is integer array, dimension (n).
!>          the pivot indices; for 1 <= j <= n, column j of the
!>          matrix has been interchanged with column jpiv(j).
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup latdf
!
!> \par further details:
!  =====================
!>
!>  this routine is a further developed implementation of algorithm
!>  bsolve in [1] using complete pivoting in the lu factorization.
!
!> \par contributors:
!  ==================
!>
!>     bo kagstrom and peter poromaa, department of computing science,
!>     umea university, s-901 87 umea, sweden.
!
!> \par references:
!  ================
!>
!> \verbatim
!>
!>
!>  [1] bo kagstrom and lars westin,
!>      generalized schur methods with condition estimators for
!>      solving the generalized sylvester equation, ieee transactions
!>      on automatic control, vol. 34, no. 7, july 1989, pp 745-751.
!>
!>  [2] peter poromaa,
!>      on efficient and robust estimators for the separation
!>      between two regular matrix pairs with applications in
!>      condition estimation. report iminf-95.05, departement of
!>      computing science, umea university, s-901 87 umea, sweden, 1995.
!> \endverbatim
!>
!  =====================================================================
subroutine dlatdf( ijob, n, z, ldz, rhs, rdsum, rdscal, ipiv, jpiv )
implicit none
!
!  -- lapack auxiliary routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
integer, intent(in) :: ijob !! ijob = 2: first compute an approximative null-vector e of z using
!! dgecon, e is normalized and solve for
integer, intent(in) :: ldz !! the leading dimension of the array z.  lda >= max(1, n).
integer, intent(in) :: n !! the number of columns of the matrix z.
real(dp), intent(inout) :: rdscal !! on entry, scaling factor used to prevent overflow in rdsum.
!! on exit, rdscal is updated w.r.t. the current contributions
real(dp), intent(inout) :: rdsum !! on entry, the sum of squares of computed contributions to the
!! dif-estimate under computation by dtgsyl, where the
!     ..
!     .. array arguments ..
integer, intent(in) :: ipiv( * ) !! the pivot indices; for 1 <= i <= n, row i of the matrix has
!! been interchanged with row ipiv(i).
integer, intent(in) :: jpiv( * ) !! the pivot indices; for 1 <= j <= n, column j of the matrix has
!! been interchanged with column jpiv(j).
real(dp), intent(inout) :: rhs( * ) !! on entry, rhs contains contributions from other subsystems.
!! on exit, rhs contains the solution of the subsystem with
real(dp), intent(in) :: z( ldz, * ) !! on entry, the lu part of the factorization of the n-by-n
!! matrix z computed by dgetc2:  z = p * l * u * q
!     ..
!
!  =====================================================================
!
!     .. parameters ..
integer :: maxdim
parameter          ( maxdim = 8 )
real(dp) :: zero, one
parameter          ( zero = 0.0d+0, one = 1.0d+0 )
!     ..
!     .. local scalars ..
integer :: i, info, j, k
real(dp) :: bm, bp, pmone, sminu, splus, temp
!     ..
!     .. local arrays ..
integer :: iwork( maxdim )
real(dp) :: work( 4*maxdim ), xm( maxdim ), xp( maxdim )
!     ..
!     .. external subroutines ..
!     ..
!     .. external functions ..
!     ..
!     .. intrinsic functions ..
intrinsic          abs, sqrt
!     ..
!     .. executable statements ..
!
if( ijob.ne.2 ) then
!
!        apply permutations ipiv to rhs
!
call dlaswp( 1, rhs, ldz, 1, n-1, ipiv, 1 )
!
!        solve for l-part choosing rhs either to +1 or -1.
!
pmone = -one
!
do 10 j = 1, n - 1
bp = rhs( j ) + one
bm = rhs( j ) - one
splus = one
!
!           look-ahead for l-part rhs(1:n-1) = + or -1, splus and
!           smin computed more efficiently than in bsolve [1].
!
splus = splus + ddot( n-j, z( j+1, j ), 1, z( j+1, j ), 1 )
sminu = ddot( n-j, z( j+1, j ), 1, rhs( j+1 ), 1 )
splus = splus*rhs( j )
if( splus.gt.sminu ) then
rhs( j ) = bp
else if( sminu.gt.splus ) then
rhs( j ) = bm
else
!
!              in this case the updating sums are equal and we can
!              choose rhs(j) +1 or -1. the first time this happens
!              we choose -1, thereafter +1. this is a simple way to
!              get good estimates of matrices like byers well-known
!              example (see [1]). (not done in bsolve.)
!
rhs( j ) = rhs( j ) + pmone
pmone = one
end if
!
!           compute the remaining r.h.s.
!
temp = -rhs( j )
call daxpy( n-j, temp, z( j+1, j ), 1, rhs( j+1 ), 1 )
!
10 continue
!
!        solve for u-part, look-ahead for rhs(n) = +-1. this is not done
!        in bsolve and will hopefully give us a better estimate because
!        any ill-conditioning of the original matrix is transferred to u
!        and not to l. u(n, n) is an approximation to sigma_min(lu).
!
call dcopy( n-1, rhs, 1, xp, 1 )
xp( n ) = rhs( n ) + one
rhs( n ) = rhs( n ) - one
splus = zero
sminu = zero
do 30 i = n, 1, -1
temp = one / z( i, i )
xp( i ) = xp( i )*temp
rhs( i ) = rhs( i )*temp
do 20 k = i + 1, n
xp( i ) = xp( i ) - xp( k )*( z( i, k )*temp )
rhs( i ) = rhs( i ) - rhs( k )*( z( i, k )*temp )
20 continue
splus = splus + abs( xp( i ) )
sminu = sminu + abs( rhs( i ) )
30 continue
if( splus.gt.sminu ) call dcopy( n, xp, 1, rhs, 1 )
!
!        apply the permutations jpiv to the computed solution (rhs)
!
call dlaswp( 1, rhs, ldz, 1, n-1, jpiv, -1 )
!
!        compute the sum of squares
!
call dlassq( n, rhs, 1, rdscal, rdsum )
!
else
!
!        ijob = 2, compute approximate nullvector xm of z
!
call dgecon( 'i', n, z, ldz, one, temp, work, iwork, info )
call dcopy( n, work( n+1 ), 1, xm, 1 )
!
!        compute rhs
!
call dlaswp( 1, xm, ldz, 1, n-1, ipiv, -1 )
temp = one / sqrt( ddot( n, xm, 1, xm, 1 ) )
call dscal( n, temp, xm, 1 )
call dcopy( n, xm, 1, xp, 1 )
call daxpy( n, one, rhs, 1, xp, 1 )
call daxpy( n, -one, xm, 1, rhs, 1 )
call dgesc2( n, z, ldz, rhs, ipiv, jpiv, temp )
call dgesc2( n, z, ldz, xp, ipiv, jpiv, temp )
if( dasum( n, xp, 1 ).gt.dasum( n, rhs, 1 ) ) call dcopy( n, xp, 1, rhs, 1 )
!
!        compute the sum of squares
!
call dlassq( n, rhs, 1, rdscal, rdsum )
!
end if
!
return
!
!     end of dlatdf
!
end subroutine dlatdf
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/iladlr.f
!> \brief \b iladlr scans a matrix for its last non-zero row.
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download iladlr + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/iladlr.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/iladlr.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/iladlr.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       integer function iladlr( m, n, a, lda )
!
!       .. scalar arguments ..
!       integer            m, n, lda
!       ..
!       .. array arguments ..
!       double precision   a( lda, * )
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> iladlr scans a for its last non-zero row.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] m
!> \verbatim
!>          m is integer
!>          the number of rows of the matrix a.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is integer
!>          the number of columns of the matrix a.
!> \endverbatim
!>
!> \param[in] a
!> \verbatim
!>          a is double precision array, dimension (lda,n)
!>          the m by n matrix a.
!> \endverbatim
!>
!> \param[in] lda
!> \verbatim
!>          lda is integer
!>          the leading dimension of the array a. lda >= max(1,m).
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup ilalr
!
!  =====================================================================
integer function iladlr( m, n, a, lda ) result(lapack_result)
implicit none
!
!  -- lapack auxiliary routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
integer, intent(in) :: m !! the number of rows of the matrix a.
integer, intent(in) :: n !! the number of columns of the matrix a.
integer, intent(in) :: lda !! the leading dimension of the array a. lda >= max(1,m).
!     ..
!     .. array arguments ..
real(dp), intent(in) :: a( lda, * ) !! the m by n matrix a.
!     ..
!
!  =====================================================================
!
!     .. parameters ..
real(dp) :: zero
parameter ( zero = 0.0d+0 )
!     ..
!     .. local scalars ..
integer :: i, j
!     ..
!     .. executable statements ..
!
!     quick test for the common case where one corner is non-zero.
if( m.eq.0 ) then
lapack_result = m
else if( a(m, 1).ne.zero .or. a(m, n).ne.zero ) then
lapack_result = m
else
!     scan up each column tracking the last zero row seen.
lapack_result = 0
do j = 1, n
i=m
do while((a(max(i,1),j).eq.zero).and.(i.ge.1))
i=i-1
enddo
lapack_result = max( lapack_result, i )
end do
end if
return
end function iladlr
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/iladlc.f
!> \brief \b iladlc scans a matrix for its last non-zero column.
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download iladlc + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/iladlc.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/iladlc.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/iladlc.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       integer function iladlc( m, n, a, lda )
!
!       .. scalar arguments ..
!       integer            m, n, lda
!       ..
!       .. array arguments ..
!       double precision   a( lda, * )
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> iladlc scans a for its last non-zero column.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] m
!> \verbatim
!>          m is integer
!>          the number of rows of the matrix a.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is integer
!>          the number of columns of the matrix a.
!> \endverbatim
!>
!> \param[in] a
!> \verbatim
!>          a is double precision array, dimension (lda,n)
!>          the m by n matrix a.
!> \endverbatim
!>
!> \param[in] lda
!> \verbatim
!>          lda is integer
!>          the leading dimension of the array a. lda >= max(1,m).
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup ilalc
!
!  =====================================================================
integer function iladlc( m, n, a, lda ) result(lapack_result)
implicit none
!
!  -- lapack auxiliary routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
integer, intent(in) :: m !! the number of rows of the matrix a.
integer, intent(in) :: n !! the number of columns of the matrix a.
integer, intent(in) :: lda !! the leading dimension of the array a. lda >= max(1,m).
!     ..
!     .. array arguments ..
real(dp), intent(in) :: a( lda, * ) !! the m by n matrix a.
!     ..
!
!  =====================================================================
!
!     .. parameters ..
real(dp) :: zero
parameter ( zero = 0.0d+0 )
!     ..
!     .. local scalars ..
integer :: i
!     ..
!     .. executable statements ..
!
!     quick test for the common case where one corner is non-zero.
if( n.eq.0 ) then
lapack_result = n
else if( a(1, n).ne.zero .or. a(m, n).ne.zero ) then
lapack_result = n
else
!     now scan each column from the end, returning with the first non-zero.
do lapack_result = n, 1, -1
do i = 1, m
if( a(i, lapack_result).ne.zero ) return
end do
end do
end if
return
end function iladlc
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlaswp.f
!> \brief \b dlaswp performs a series of row interchanges on a general rectangular matrix.
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download dlaswp + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaswp.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaswp.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaswp.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       subroutine dlaswp( n, a, lda, k1, k2, ipiv, incx )
!
!       .. scalar arguments ..
!       integer            incx, k1, k2, lda, n
!       ..
!       .. array arguments ..
!       integer            ipiv( * )
!       double precision   a( lda, * )
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dlaswp performs a series of row interchanges on the matrix a.
!> one row interchange is initiated for each of rows k1 through k2 of a.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] n
!> \verbatim
!>          n is integer
!>          the number of columns of the matrix a.
!> \endverbatim
!>
!> \param[in,out] a
!> \verbatim
!>          a is double precision array, dimension (lda,n)
!>          on entry, the matrix of column dimension n to which the row
!>          interchanges will be applied.
!>          on exit, the permuted matrix.
!> \endverbatim
!>
!> \param[in] lda
!> \verbatim
!>          lda is integer
!>          the leading dimension of the array a.
!> \endverbatim
!>
!> \param[in] k1
!> \verbatim
!>          k1 is integer
!>          the first element of ipiv for which a row interchange will
!>          be done.
!> \endverbatim
!>
!> \param[in] k2
!> \verbatim
!>          k2 is integer
!>          (k2-k1+1) is the number of elements of ipiv for which a row
!>          interchange will be done.
!> \endverbatim
!>
!> \param[in] ipiv
!> \verbatim
!>          ipiv is integer array, dimension (k1+(k2-k1)*abs(incx))
!>          the vector of pivot indices. only the elements in positions
!>          k1 through k1+(k2-k1)*abs(incx) of ipiv are accessed.
!>          ipiv(k1+(k-k1)*abs(incx)) = l implies rows k and l are to be
!>          interchanged.
!> \endverbatim
!>
!> \param[in] incx
!> \verbatim
!>          incx is integer
!>          the increment between successive values of ipiv. if incx
!>          is negative, the pivots are applied in reverse order.
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup laswp
!
!> \par further details:
!  =====================
!>
!> \verbatim
!>
!>  modified by
!>   r. c. whaley, computer science dept., univ. of tenn., knoxville, usa
!> \endverbatim
!>
!  =====================================================================
subroutine dlaswp( n, a, lda, k1, k2, ipiv, incx )
implicit none
!
!  -- lapack auxiliary routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
integer, intent(in) :: incx !! the increment between successive values of ipiv. if incx is
!! negative, the pivots are applied in reverse order.
integer, intent(in) :: k1 !! the first element of ipiv for which a row interchange will be done.
integer, intent(in) :: k2 !! (k2-k1+1) is the number of elements of ipiv for which a row
!! interchange will be done.
integer, intent(in) :: lda !! the leading dimension of the array a.
integer, intent(in) :: n !! the number of columns of the matrix a.
!     ..
!     .. array arguments ..
integer, intent(in) :: ipiv( * ) !! the vector of pivot indices. only the elements in positions k1
!! through k1+(k2-k1)*abs(incx) of ipiv are accessed.
real(dp), intent(inout) :: a( lda, * ) !! on entry, the matrix of column dimension n to which the
!! row interchanges will be applied. on exit, the permuted
!! matrix.
!     ..
!
! =====================================================================
!
!     .. local scalars ..
integer :: i, i1, i2, inc, ip, ix, ix0, j, k, n32
real(dp) :: temp
!     ..
!     .. executable statements ..
!
!     interchange row i with row ipiv(k1+(i-k1)*abs(incx)) for each of rows
!     k1 through k2.
!
if( incx.gt.0 ) then
ix0 = k1
i1 = k1
i2 = k2
inc = 1
else if( incx.lt.0 ) then
ix0 = k1 + ( k1-k2 )*incx
i1 = k2
i2 = k1
inc = -1
else
return
end if
!
n32 = ( n / 32 )*32
if( n32.ne.0 ) then
do 30 j = 1, n32, 32
ix = ix0
do 20 i = i1, i2, inc
ip = ipiv( ix )
if( ip.ne.i ) then
do 10 k = j, j + 31
temp = a( i, k )
a( i, k ) = a( ip, k )
a( ip, k ) = temp
10 continue
end if
ix = ix + incx
20 continue
30 continue
end if
if( n32.ne.n ) then
n32 = n32 + 1
ix = ix0
do 50 i = i1, i2, inc
ip = ipiv( ix )
if( ip.ne.i ) then
do 40 k = n32, n
temp = a( i, k )
a( i, k ) = a( ip, k )
a( ip, k ) = temp
40 continue
end if
ix = ix + incx
50 continue
end if
!
return
!
!     end of dlaswp
!
end subroutine dlaswp
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/BLAS/SRC/idamax.f
!> \brief \b idamax
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  definition:
!  ===========
!
!       integer function idamax(n,dx,incx)
!
!       .. scalar arguments ..
!       integer incx,n
!       ..
!       .. array arguments ..
!       double precision dx(*)
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!>    idamax finds the index of the first element having maximum absolute value.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] n
!> \verbatim
!>          n is integer
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in] dx
!> \verbatim
!>          dx is double precision array, dimension ( 1 + ( n - 1 )*abs( incx ) )
!> \endverbatim
!>
!> \param[in] incx
!> \verbatim
!>          incx is integer
!>         storage spacing between elements of dx
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup iamax
!
!> \par further details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 3/93 to return if incx .le. 0.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
integer function idamax(n,dx,incx) result(lapack_result)
implicit none
!
!  -- reference blas level1 routine --
!  -- reference blas is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
integer, intent(in) :: incx !! storage spacing between elements of dx
integer, intent(in) :: n !! number of elements in input vector(s)
!     ..
!     .. array arguments ..
real(dp), intent(in) :: dx(*) !! first strided vector, updated in place when requested by the operation.
!     ..
!
!  =====================================================================
!
!     .. local scalars ..
real(dp) :: dmax
integer :: i, ix
!     ..
!     .. intrinsic functions ..
intrinsic dabs
!     ..
lapack_result = 0
if (n.lt.1 .or. incx.le.0) return
lapack_result = 1
if (n.eq.1) return
if (incx.eq.1) then
!
!        code for increment equal to 1
!
dmax = dabs(dx(1))
do i = 2,n
if (dabs(dx(i)).gt.dmax) then
lapack_result = i
dmax = dabs(dx(i))
end if
end do
else
!
!        code for increment not equal to 1
!
ix = 1
dmax = dabs(dx(1))
ix = ix + incx
do i = 2,n
if (dabs(dx(ix)).gt.dmax) then
lapack_result = i
dmax = dabs(dx(ix))
end if
ix = ix + incx
end do
end if
return
!
!     end of idamax
!
end function idamax
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/BLAS/SRC/dswap.f
!> \brief \b dswap
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  definition:
!  ===========
!
!       subroutine dswap(n,dx,incx,dy,incy)
!
!       .. scalar arguments ..
!       integer incx,incy,n
!       ..
!       .. array arguments ..
!       double precision dx(*),dy(*)
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!>    dswap interchanges two vectors.
!>    uses unrolled loops for increments equal to 1.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] n
!> \verbatim
!>          n is integer
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in,out] dx
!> \verbatim
!>          dx is double precision array, dimension ( 1 + ( n - 1 )*abs( incx ) )
!> \endverbatim
!>
!> \param[in] incx
!> \verbatim
!>          incx is integer
!>         storage spacing between elements of dx
!> \endverbatim
!>
!> \param[in,out] dy
!> \verbatim
!>          dy is double precision array, dimension ( 1 + ( n - 1 )*abs( incy ) )
!> \endverbatim
!>
!> \param[in] incy
!> \verbatim
!>          incy is integer
!>         storage spacing between elements of dy
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup swap
!
!> \par further details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
subroutine dswap(n,dx,incx,dy,incy)
implicit none
!
!  -- reference blas level1 routine --
!  -- reference blas is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
integer, intent(in) :: incx !! storage spacing between elements of dx
integer, intent(in) :: incy !! storage spacing between elements of dy
integer, intent(in) :: n !! number of elements in input vector(s)
!     ..
!     .. array arguments ..
real(dp), intent(inout) :: dx(*) !! first strided vector, updated in place when requested by the operation.
real(dp), intent(inout) :: dy(*) !! second strided vector, updated in place when requested by the operation.
!     ..
!
!  =====================================================================
!
!     .. local scalars ..
real(dp) :: dtemp
integer :: i, ix, iy, m, mp1
!     ..
!     .. intrinsic functions ..
intrinsic mod
!     ..
if (n.le.0) return
if (incx.eq.1 .and. incy.eq.1) then
!
!       code for both increments equal to 1
!
!
!       clean-up loop
!
m = mod(n,3)
if (m.ne.0) then
do i = 1,m
dtemp = dx(i)
dx(i) = dy(i)
dy(i) = dtemp
end do
if (n.lt.3) return
end if
mp1 = m + 1
do i = mp1,n,3
dtemp = dx(i)
dx(i) = dy(i)
dy(i) = dtemp
dtemp = dx(i+1)
dx(i+1) = dy(i+1)
dy(i+1) = dtemp
dtemp = dx(i+2)
dx(i+2) = dy(i+2)
dy(i+2) = dtemp
end do
else
!
!       code for unequal increments or equal increments not equal
!         to 1
!
ix = 1
iy = 1
if (incx.lt.0) ix = (-n+1)*incx + 1
if (incy.lt.0) iy = (-n+1)*incy + 1
do i = 1,n
dtemp = dx(ix)
dx(ix) = dy(iy)
dy(iy) = dtemp
ix = ix + incx
iy = iy + incy
end do
end if
return
!
!     end of dswap
!
end subroutine dswap
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dgecon.f
!> \brief \b dgecon
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download dgecon + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgecon.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgecon.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgecon.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       subroutine dgecon( norm, n, a, lda, anorm, rcond, work, iwork,
!                          info )
!
!       .. scalar arguments ..
!       character          norm
!       integer            info, lda, n
!       double precision   anorm, rcond
!       ..
!       .. array arguments ..
!       integer            iwork( * )
!       double precision   a( lda, * ), work( * )
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dgecon estimates the reciprocal of the condition number of a general
!> real matrix a, in either the 1-norm or the infinity-norm, using
!> the lu factorization computed by dgetrf.
!>
!> an estimate is obtained for norm(inv(a)), and the reciprocal of the
!> condition number is computed as
!>    rcond = 1 / ( norm(a) * norm(inv(a)) ).
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] norm
!> \verbatim
!>          norm is character*1
!>          specifies whether the 1-norm condition number or the
!>          infinity-norm condition number is required:
!>          = '1' or 'o':  1-norm;
!>          = 'i':         infinity-norm.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is integer
!>          the order of the matrix a.  n >= 0.
!> \endverbatim
!>
!> \param[in] a
!> \verbatim
!>          a is double precision array, dimension (lda,n)
!>          the factors l and u from the factorization a = p*l*u
!>          as computed by dgetrf.
!> \endverbatim
!>
!> \param[in] lda
!> \verbatim
!>          lda is integer
!>          the leading dimension of the array a.  lda >= max(1,n).
!> \endverbatim
!>
!> \param[in] anorm
!> \verbatim
!>          anorm is double precision
!>          if norm = '1' or 'o', the 1-norm of the original matrix a.
!>          if norm = 'i', the infinity-norm of the original matrix a.
!> \endverbatim
!>
!> \param[out] rcond
!> \verbatim
!>          rcond is double precision
!>          the reciprocal of the condition number of the matrix a,
!>          computed as rcond = 1/(norm(a) * norm(inv(a))).
!> \endverbatim
!>
!> \param[out] work
!> \verbatim
!>          work is double precision array, dimension (4*n)
!> \endverbatim
!>
!> \param[out] iwork
!> \verbatim
!>          iwork is integer array, dimension (n)
!> \endverbatim
!>
!> \param[out] info
!> \verbatim
!>          info is integer
!>          = 0:  successful exit
!>          < 0:  if info = -i, the i-th argument had an illegal value.
!>                nans are illegal values for anorm, and they propagate to
!>                the output parameter rcond.
!>                infinity is illegal for anorm, and it propagates to the output
!>                parameter rcond as 0.
!>          = 1:  if rcond = nan, or
!>                   rcond = inf, or
!>                   the computed norm of the inverse of a is 0.
!>                in the latter, rcond = 0 is returned.
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup gecon
!
!  =====================================================================
subroutine dgecon( norm, n, a, lda, anorm, rcond, work, iwork, info )
implicit none
!
!  -- lapack computational routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
character, intent(in) :: norm !! specifies whether the 1-norm condition number or the infinity-
!! norm condition number is required:
integer, intent(out) :: info !! zero success; negative invalid argument or positive numerical failure.
integer, intent(in) :: lda !! the leading dimension of the array a.  lda >= max(1,n).
integer, intent(in) :: n !! the order of the matrix a.  n >= 0.
real(dp), intent(in) :: anorm !! if norm = '1' or 'o', the 1-norm of the original matrix a. if
!! norm = 'i', the infinity-norm of the original matrix a.
real(dp), intent(out) :: rcond !! the reciprocal of the condition number of the matrix a, computed
!! as rcond = 1/(norm(a) * norm(inv(a))).
!     ..
!     .. array arguments ..
integer, intent(out) :: iwork( * ) !! integer scratch and pivot indices for the active small systems.
real(dp), intent(in) :: a( lda, * ) !! the factors l and u from the factorization a = p*l*u as
!! computed by dgetrf.
real(dp), intent(out) :: work( * ) !! caller-owned scratch array, sized according to the retained contract.
!     ..
!
!  =====================================================================
!
!     .. parameters ..
real(dp) :: one, zero
parameter          ( one = 1.0d+0, zero = 0.0d+0 )
!     ..
!     .. local scalars ..
logical :: onenrm
character :: normin
integer :: ix, kase, kase1
real(dp) :: ainvnm, scale, sl, smlnum, su, hugeval
!     ..
!     .. local arrays ..
integer :: isave( 3 )
!     ..
!     .. external functions ..
!     ..
!     .. external subroutines ..
!     ..
!     .. intrinsic functions ..
intrinsic          abs, max
!     ..
!     .. executable statements ..
!
hugeval = dlamch( 'overflow' )
!
!     test the input parameters.
!
info = 0
onenrm = norm.eq.'1' .or. lsame( norm, 'o' )
if( .not.onenrm .and. .not.lsame( norm, 'i' ) ) then
info = -1
else if( n.lt.0 ) then
info = -2
else if( lda.lt.max( 1, n ) ) then
info = -4
else if( anorm.lt.zero ) then
info = -5
end if
if( info.ne.0 ) then
call xerbla( 'dgecon', -info )
return
end if
!
!     quick return if possible
!
rcond = zero
if( n.eq.0 ) then
rcond = one
return
else if( anorm.eq.zero ) then
return
else if( disnan( anorm ) ) then
rcond = anorm
info = -5
return
else if( anorm.gt.hugeval ) then
info = -5
return
end if
!
smlnum = dlamch( 'safe minimum' )
!
!     estimate the norm of inv(a).
!
ainvnm = zero
normin = 'n'
if( onenrm ) then
kase1 = 1
else
kase1 = 2
end if
kase = 0
10 continue
call dlacn2( n, work( n+1 ), work, iwork, ainvnm, kase, isave )
if( kase.ne.0 ) then
if( kase.eq.kase1 ) then
!
!           multiply by inv(l).
!
call dlatrs( 'lower', 'no transpose', 'unit', normin, n, a, lda, work, sl, work( 2*n+1 ), info )
!
!           multiply by inv(u).
!
call dlatrs( 'upper', 'no transpose', 'non-unit', normin, n, a, lda, work, su, work( 3*n+1 ), info )
else
!
!           multiply by inv(u**t).
!
call dlatrs( 'upper', 'transpose', 'non-unit', normin, n, a, lda, work, su, work( 3*n+1 ), info )
!
!           multiply by inv(l**t).
!
call dlatrs( 'lower', 'transpose', 'unit', normin, n, a, lda, work, sl, work( 2*n+1 ), info )
end if
!
!        divide x by 1/(sl*su) if doing so will not cause overflow.
!
scale = sl*su
normin = 'y'
if( scale.ne.one ) then
ix = idamax( n, work, 1 )
if( scale.lt.abs( work( ix ) )*smlnum .or. scale.eq.zero ) go to 20
call drscl( n, scale, work, 1 )
end if
go to 10
end if
!
!     compute the estimate of the reciprocal condition number.
!
if( ainvnm.ne.zero ) then
rcond = ( one / ainvnm ) / anorm
else
info = 1
return
end if
!
!     check for nans and infs
!
if( disnan( rcond ) .or. rcond.gt.hugeval ) info = 1
!
20 continue
return
!
!     end of dgecon
!
end subroutine dgecon
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/BLAS/SRC/dasum.f
!> \brief \b dasum
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  definition:
!  ===========
!
!       double precision function dasum(n,dx,incx)
!
!       .. scalar arguments ..
!       integer incx,n
!       ..
!       .. array arguments ..
!       double precision dx(*)
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!>    dasum takes the sum of the absolute values.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] n
!> \verbatim
!>          n is integer
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in] dx
!> \verbatim
!>          dx is double precision array, dimension ( 1 + ( n - 1 )*abs( incx ) )
!> \endverbatim
!>
!> \param[in] incx
!> \verbatim
!>          incx is integer
!>         storage spacing between elements of dx
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup asum
!
!> \par further details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 3/93 to return if incx .le. 0.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
real(dp) function dasum(n,dx,incx) result(lapack_result)
implicit none
!
!  -- reference blas level1 routine --
!  -- reference blas is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
integer, intent(in) :: incx !! storage spacing between elements of dx
integer, intent(in) :: n !! number of elements in input vector(s)
!     ..
!     .. array arguments ..
real(dp), intent(in) :: dx(*) !! first strided vector, updated in place when requested by the operation.
!     ..
!
!  =====================================================================
!
!     .. local scalars ..
real(dp) :: dtemp
integer :: i, m, mp1, nincx
!     ..
!     .. intrinsic functions ..
intrinsic dabs,mod
!     ..
lapack_result = 0.0d0
dtemp = 0.0d0
if (n.le.0 .or. incx.le.0) return
if (incx.eq.1) then
!        code for increment equal to 1
!
!
!        clean-up loop
!
m = mod(n,6)
if (m.ne.0) then
do i = 1,m
dtemp = dtemp + dabs(dx(i))
end do
if (n.lt.6) then
lapack_result = dtemp
return
end if
end if
mp1 = m + 1
do i = mp1,n,6
dtemp = dtemp + dabs(dx(i)) + dabs(dx(i+1)) + dabs(dx(i+2)) + dabs(dx(i+3)) + dabs(dx(i+4)) + dabs(dx(i+5))
end do
else
!
!        code for increment not equal to 1
!
nincx = n*incx
do i = 1,nincx,incx
dtemp = dtemp + dabs(dx(i))
end do
end if
lapack_result = dtemp
return
!
!     end of dasum
!
end function dasum
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/BLAS/SRC/ddot.f
!> \brief \b ddot
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  definition:
!  ===========
!
!       double precision function ddot(n,dx,incx,dy,incy)
!
!       .. scalar arguments ..
!       integer incx,incy,n
!       ..
!       .. array arguments ..
!       double precision dx(*),dy(*)
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!>    ddot forms the dot product of two vectors.
!>    uses unrolled loops for increments equal to one.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] n
!> \verbatim
!>          n is integer
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in] dx
!> \verbatim
!>          dx is double precision array, dimension ( 1 + ( n - 1 )*abs( incx ) )
!> \endverbatim
!>
!> \param[in] incx
!> \verbatim
!>          incx is integer
!>         storage spacing between elements of dx
!> \endverbatim
!>
!> \param[in] dy
!> \verbatim
!>          dy is double precision array, dimension ( 1 + ( n - 1 )*abs( incy ) )
!> \endverbatim
!>
!> \param[in] incy
!> \verbatim
!>          incy is integer
!>         storage spacing between elements of dy
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup dot
!
!> \par further details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
real(dp) function ddot(n,dx,incx,dy,incy) result(lapack_result)
implicit none
!
!  -- reference blas level1 routine --
!  -- reference blas is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
integer, intent(in) :: incx !! storage spacing between elements of dx
integer, intent(in) :: incy !! storage spacing between elements of dy
integer, intent(in) :: n !! number of elements in input vector(s)
!     ..
!     .. array arguments ..
real(dp), intent(in) :: dx(*) !! first strided vector, updated in place when requested by the operation.
real(dp), intent(in) :: dy(*) !! second strided vector, updated in place when requested by the operation.
!     ..
!
!  =====================================================================
!
!     .. local scalars ..
real(dp) :: dtemp
integer :: i, ix, iy, m, mp1
!     ..
!     .. intrinsic functions ..
intrinsic mod
!     ..
lapack_result = 0.0d0
dtemp = 0.0d0
if (n.le.0) return
if (incx.eq.1 .and. incy.eq.1) then
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
m = mod(n,5)
if (m.ne.0) then
do i = 1,m
dtemp = dtemp + dx(i)*dy(i)
end do
if (n.lt.5) then
lapack_result=dtemp
return
end if
end if
mp1 = m + 1
do i = mp1,n,5
dtemp = dtemp + dx(i)*dy(i) + dx(i+1)*dy(i+1) + dx(i+2)*dy(i+2) + dx(i+3)*dy(i+3) + dx(i+4)*dy(i+4)
end do
else
!
!        code for unequal increments or equal increments
!          not equal to 1
!
ix = 1
iy = 1
if (incx.lt.0) ix = (-n+1)*incx + 1
if (incy.lt.0) iy = (-n+1)*incy + 1
do i = 1,n
dtemp = dtemp + dx(ix)*dy(iy)
ix = ix + incx
iy = iy + incy
end do
end if
lapack_result = dtemp
return
!
!     end of ddot
!
end function ddot
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlacn2.f
!> \brief \b dlacn2 estimates the 1-norm of a square matrix, using reverse communication for
!evaluating matrix-vector products.
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download dlacn2 + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlacn2.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlacn2.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlacn2.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       subroutine dlacn2( n, v, x, isgn, est, kase, isave )
!
!       .. scalar arguments ..
!       integer            kase, n
!       double precision   est
!       ..
!       .. array arguments ..
!       integer            isgn( * ), isave( 3 )
!       double precision   v( * ), x( * )
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dlacn2 estimates the 1-norm of a square, real matrix a.
!> reverse communication is used for evaluating matrix-vector products.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] n
!> \verbatim
!>          n is integer
!>         the order of the matrix.  n >= 1.
!> \endverbatim
!>
!> \param[out] v
!> \verbatim
!>          v is double precision array, dimension (n)
!>         on the final return, v = a*w,  where  est = norm(v)/norm(w)
!>         (w is not returned).
!> \endverbatim
!>
!> \param[in,out] x
!> \verbatim
!>          x is double precision array, dimension (n)
!>         on an intermediate return, x should be overwritten by
!>               a * x,   if kase=1,
!>               a**t * x,  if kase=2,
!>         and dlacn2 must be re-called with all the other parameters
!>         unchanged.
!> \endverbatim
!>
!> \param[out] isgn
!> \verbatim
!>          isgn is integer array, dimension (n)
!> \endverbatim
!>
!> \param[in,out] est
!> \verbatim
!>          est is double precision
!>         on entry with kase = 1 or 2 and isave(1) = 3, est should be
!>         unchanged from the previous call to dlacn2.
!>         on exit, est is an estimate (a lower bound) for norm(a).
!> \endverbatim
!>
!> \param[in,out] kase
!> \verbatim
!>          kase is integer
!>         on the initial call to dlacn2, kase should be 0.
!>         on an intermediate return, kase will be 1 or 2, indicating
!>         whether x should be overwritten by a * x  or a**t * x.
!>         on the final return from dlacn2, kase will again be 0.
!> \endverbatim
!>
!> \param[in,out] isave
!> \verbatim
!>          isave is integer array, dimension (3)
!>         isave is used to save variables between calls to dlacn2
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup lacn2
!
!> \par further details:
!  =====================
!>
!> \verbatim
!>
!>  originally named sonest, dated march 16, 1988.
!>
!>  this is a thread safe version of dlacon, which uses the array isave
!>  in place of a save statement, as follows:
!>
!>     dlacon     dlacn2
!>      jump     isave(1)
!>      j        isave(2)
!>      iter     isave(3)
!> \endverbatim
!
!> \par contributors:
!  ==================
!>
!>     nick higham, university of manchester
!
!> \par references:
!  ================
!>
!>  n.j. higham, "fortran codes for estimating the one-norm of
!>  a real or complex matrix, with applications to condition estimation",
!>  acm trans. math. soft., vol. 14, no. 4, pp. 381-396, december 1988.
!>
!  =====================================================================
subroutine dlacn2( n, v, x, isgn, est, kase, isave )
implicit none
!
!  -- lapack auxiliary routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
integer, intent(inout) :: kase !! on the initial call to dlacn2, kase should be 0. on an
!! intermediate return, kase will be 1 or 2, indicating
integer, intent(in) :: n !! the order of the matrix.  n >= 1.
real(dp), intent(inout) :: est !! on entry with kase = 1 or 2 and isave(1) = 3, est should be
!! unchanged from the previous call to dlacn2.
!     ..
!     .. array arguments ..
integer, intent(out) :: isgn( * ) !! saved component signs for reverse-communication norm estimation.
integer, intent(inout) :: isave( 3 ) !! caller-owned reverse-communication iteration history, length three.
real(dp), intent(out) :: v( * ) !! on the final return, v = a*w,  where  est = norm(v)/norm(w) (w
!! is not returned).
real(dp), intent(inout) :: x( * ) !! on an intermediate return, x should be overwritten by a * x,
!! if kase=1, a**t * x,  if kase=2,
!     ..
!
!  =====================================================================
!
!     .. parameters ..
integer :: itmax
parameter          ( itmax = 5 )
real(dp) :: zero, one, two
parameter          ( zero = 0.0d+0, one = 1.0d+0, two = 2.0d+0 )
!     ..
!     .. local scalars ..
integer :: i, jlast
real(dp) :: altsgn, estold, temp, xs
!     ..
!     .. external functions ..
!     ..
!     .. external subroutines ..
!     ..
!     .. intrinsic functions ..
intrinsic          abs, dble, nint
!     ..
!     .. executable statements ..
!
if( kase.eq.0 ) then
do 10 i = 1, n
x( i ) = one / dble( n )
10 continue
kase = 1
isave( 1 ) = 1
return
end if
!
go to ( 20, 40, 70, 110, 140 )isave( 1 )
!
!     ................ entry   (isave( 1 ) = 1)
!     first iteration.  x has been overwritten by a*x.
!
20 continue
if( n.eq.1 ) then
v( 1 ) = x( 1 )
est = abs( v( 1 ) )
!        ... quit
go to 150
end if
est = dasum( n, x, 1 )
!
do 30 i = 1, n
if( x(i).ge.zero ) then
x(i) = one
else
x(i) = -one
end if
isgn( i ) = nint( x( i ) )
30 continue
kase = 2
isave( 1 ) = 2
return
!
!     ................ entry   (isave( 1 ) = 2)
!     first iteration.  x has been overwritten by transpose(a)*x.
!
40 continue
isave( 2 ) = idamax( n, x, 1 )
isave( 3 ) = 2
!
!     main loop - iterations 2,3,...,itmax.
!
50 continue
do 60 i = 1, n
x( i ) = zero
60 continue
x( isave( 2 ) ) = one
kase = 1
isave( 1 ) = 3
return
!
!     ................ entry   (isave( 1 ) = 3)
!     x has been overwritten by a*x.
!
70 continue
call dcopy( n, x, 1, v, 1 )
estold = est
est = dasum( n, v, 1 )
do 80 i = 1, n
if( x(i).ge.zero ) then
xs = one
else
xs = -one
end if
if( nint( xs ).ne.isgn( i ) ) go to 90
80 continue
!     repeated sign vector detected, hence algorithm has converged.
go to 120
!
90 continue
!     test for cycling.
if( est.le.estold ) go to 120
!
do 100 i = 1, n
if( x(i).ge.zero ) then
x(i) = one
else
x(i) = -one
end if
isgn( i ) = nint( x( i ) )
100 continue
kase = 2
isave( 1 ) = 4
return
!
!     ................ entry   (isave( 1 ) = 4)
!     x has been overwritten by transpose(a)*x.
!
110 continue
jlast = isave( 2 )
isave( 2 ) = idamax( n, x, 1 )
if( ( x( jlast ).ne.abs( x( isave( 2 ) ) ) ) .and. ( isave( 3 ).lt.itmax ) ) then
isave( 3 ) = isave( 3 ) + 1
go to 50
end if
!
!     iteration complete.  final stage.
!
120 continue
altsgn = one
do 130 i = 1, n
x( i ) = altsgn*( one+dble( i-1 ) / dble( n-1 ) )
altsgn = -altsgn
130 continue
kase = 1
isave( 1 ) = 5
return
!
!     ................ entry   (isave( 1 ) = 5)
!     x has been overwritten by a*x.
!
140 continue
temp = two*( dasum( n, x, 1 ) / dble( 3*n ) )
if( temp.gt.est ) then
call dcopy( n, x, 1, v, 1 )
est = temp
end if
!
150 continue
kase = 0
return
!
!     end of dlacn2
!
end subroutine dlacn2
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlatrs.f
!> \brief \b dlatrs solves a triangular system of equations with the scale factor set to prevent
!overflow.
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download dlatrs + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlatrs.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlatrs.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlatrs.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       subroutine dlatrs( uplo, trans, diag, normin, n, a, lda, x, scale,
!                          cnorm, info )
!
!       .. scalar arguments ..
!       character          diag, normin, trans, uplo
!       integer            info, lda, n
!       double precision   scale
!       ..
!       .. array arguments ..
!       double precision   a( lda, * ), cnorm( * ), x( * )
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dlatrs solves one of the triangular systems
!>
!>    a *x = s*b  or  a**t *x = s*b
!>
!> with scaling to prevent overflow.  here a is an upper or lower
!> triangular matrix, a**t denotes the transpose of a, x and b are
!> n-element vectors, and s is a scaling factor, usually less than
!> or equal to 1, chosen so that the components of x will be less than
!> the overflow threshold.  if the unscaled problem will not cause
!> overflow, the level 2 blas routine dtrsv is called.  if the matrix a
!> is singular (a(j,j) = 0 for some j), then s is set to 0 and a
!> non-trivial solution to a*x = 0 is returned.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] uplo
!> \verbatim
!>          uplo is character*1
!>          specifies whether the matrix a is upper or lower triangular.
!>          = 'u':  upper triangular
!>          = 'l':  lower triangular
!> \endverbatim
!>
!> \param[in] trans
!> \verbatim
!>          trans is character*1
!>          specifies the operation applied to a.
!>          = 'n':  solve a * x = s*b  (no transpose)
!>          = 't':  solve a**t* x = s*b  (transpose)
!>          = 'c':  solve a**t* x = s*b  (conjugate transpose = transpose)
!> \endverbatim
!>
!> \param[in] diag
!> \verbatim
!>          diag is character*1
!>          specifies whether or not the matrix a is unit triangular.
!>          = 'n':  non-unit triangular
!>          = 'u':  unit triangular
!> \endverbatim
!>
!> \param[in] normin
!> \verbatim
!>          normin is character*1
!>          specifies whether cnorm has been set or not.
!>          = 'y':  cnorm contains the column norms on entry
!>          = 'n':  cnorm is not set on entry.  on exit, the norms will
!>                  be computed and stored in cnorm.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is integer
!>          the order of the matrix a.  n >= 0.
!> \endverbatim
!>
!> \param[in] a
!> \verbatim
!>          a is double precision array, dimension (lda,n)
!>          the triangular matrix a.  if uplo = 'u', the leading n by n
!>          upper triangular part of the array a contains the upper
!>          triangular matrix, and the strictly lower triangular part of
!>          a is not referenced.  if uplo = 'l', the leading n by n lower
!>          triangular part of the array a contains the lower triangular
!>          matrix, and the strictly upper triangular part of a is not
!>          referenced.  if diag = 'u', the diagonal elements of a are
!>          also not referenced and are assumed to be 1.
!> \endverbatim
!>
!> \param[in] lda
!> \verbatim
!>          lda is integer
!>          the leading dimension of the array a.  lda >= max (1,n).
!> \endverbatim
!>
!> \param[in,out] x
!> \verbatim
!>          x is double precision array, dimension (n)
!>          on entry, the right hand side b of the triangular system.
!>          on exit, x is overwritten by the solution vector x.
!> \endverbatim
!>
!> \param[out] scale
!> \verbatim
!>          scale is double precision
!>          the scaling factor s for the triangular system
!>             a * x = s*b  or  a**t* x = s*b.
!>          if scale = 0, the matrix a is singular or badly scaled, and
!>          the vector x is an exact or approximate solution to a*x = 0.
!> \endverbatim
!>
!> \param[in,out] cnorm
!> \verbatim
!>          cnorm is double precision array, dimension (n)
!>
!>          if normin = 'y', cnorm is an input argument and cnorm(j)
!>          contains the norm of the off-diagonal part of the j-th column
!>          of a.  if trans = 'n', cnorm(j) must be greater than or equal
!>          to the infinity-norm, and if trans = 't' or 'c', cnorm(j)
!>          must be greater than or equal to the 1-norm.
!>
!>          if normin = 'n', cnorm is an output argument and cnorm(j)
!>          returns the 1-norm of the offdiagonal part of the j-th column
!>          of a.
!> \endverbatim
!>
!> \param[out] info
!> \verbatim
!>          info is integer
!>          = 0:  successful exit
!>          < 0:  if info = -k, the k-th argument had an illegal value
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup latrs
!
!> \par further details:
!  =====================
!>
!> \verbatim
!>
!>  a rough bound on x is computed; if that is less than overflow, dtrsv
!>  is called, otherwise, specific code is used which checks for possible
!>  overflow or divide-by-zero at every operation.
!>
!>  a columnwise scheme is used for solving a*x = b.  the basic algorithm
!>  if a is lower triangular is
!>
!>       x[1:n] := b[1:n]
!>       for j = 1, ..., n
!>            x(j) := x(j) / a(j,j)
!>            x[j+1:n] := x[j+1:n] - x(j) * a[j+1:n,j]
!>       end
!>
!>  define bounds on the components of x after j iterations of the loop:
!>     m(j) = bound on x[1:j]
!>     g(j) = bound on x[j+1:n]
!>  initially, let m(0) = 0 and g(0) = max{x(i), i=1,...,n}.
!>
!>  then for iteration j+1 we have
!>     m(j+1) <= g(j) / | a(j+1,j+1) |
!>     g(j+1) <= g(j) + m(j+1) * | a[j+2:n,j+1] |
!>            <= g(j) ( 1 + cnorm(j+1) / | a(j+1,j+1) | )
!>
!>  where cnorm(j+1) is greater than or equal to the infinity-norm of
!>  column j+1 of a, not counting the diagonal.  hence
!>
!>     g(j) <= g(0) product ( 1 + cnorm(i) / | a(i,i) | )
!>                  1<=i<=j
!>  and
!>
!>     |x(j)| <= ( g(0) / |a(j,j)| ) product ( 1 + cnorm(i) / |a(i,i)| )
!>                                   1<=i< j
!>
!>  since |x(j)| <= m(j), we use the level 2 blas routine dtrsv if the
!>  reciprocal of the largest m(j), j=1,..,n, is larger than
!>  max(underflow, 1/overflow).
!>
!>  the bound on x(j) is also used to determine when a step in the
!>  columnwise method can be performed without fear of overflow.  if
!>  the computed bound is greater than a large constant, x is scaled to
!>  prevent overflow, but if the bound overflows, x is set to 0, x(j) to
!>  1, and scale to 0, and a non-trivial solution to a*x = 0 is found.
!>
!>  similarly, a row-wise scheme is used to solve a**t*x = b.  the basic
!>  algorithm for a upper triangular is
!>
!>       for j = 1, ..., n
!>            x(j) := ( b(j) - a[1:j-1,j]**t * x[1:j-1] ) / a(j,j)
!>       end
!>
!>  we simultaneously compute two bounds
!>       g(j) = bound on ( b(i) - a[1:i-1,i]**t * x[1:i-1] ), 1<=i<=j
!>       m(j) = bound on x(i), 1<=i<=j
!>
!>  the initial values are g(0) = 0, m(0) = max{b(i), i=1,..,n}, and we
!>  add the constraint g(j) >= g(j-1) and m(j) >= m(j-1) for j >= 1.
!>  then the bound on x(j) is
!>
!>       m(j) <= m(j-1) * ( 1 + cnorm(j) ) / | a(j,j) |
!>
!>            <= m(0) * product ( ( 1 + cnorm(i) ) / |a(i,i)| )
!>                      1<=i<=j
!>
!>  and we can safely call dtrsv if 1/m(n) and 1/g(n) are both greater
!>  than max(underflow, 1/overflow).
!> \endverbatim
!>
!  =====================================================================
subroutine dlatrs( uplo, trans, diag, normin, n, a, lda, x, scale, cnorm, info )
implicit none
!
!  -- lapack auxiliary routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
character, intent(in) :: diag !! specifies whether or not the matrix a is unit triangular.
character, intent(in) :: normin !! specifies whether cnorm has been set or not.
character, intent(in) :: trans !! specifies the operation applied to a.
character, intent(in) :: uplo !! specifies whether the matrix a is upper or lower triangular.
integer, intent(out) :: info !! zero success; negative invalid argument or positive numerical failure.
integer, intent(in) :: lda !! the leading dimension of the array a.  lda >= max (1,n).
integer, intent(in) :: n !! the order of the matrix a.  n >= 0.
real(dp), intent(out) :: scale !! the scaling factor s for the triangular system a * x = s*b  or
!! a**t* x = s*b. if scale = 0, the matrix a is singular or badly
!! scaled, and
!     ..
!     .. array arguments ..
real(dp), intent(in) :: a( lda, * ) !! the triangular matrix a.  if uplo = 'u', the leading n by n
!! upper triangular part of the array a contains the upper
real(dp), intent(inout) :: cnorm( * ) !! if normin = 'y', cnorm is an input argument and cnorm(j)
!! contains the norm of the off-diagonal part of the j-th
!! column
real(dp), intent(inout) :: x( * ) !! on entry, the right hand side b of the triangular system. on
!! exit, x is overwritten by the solution vector x.
!     ..
!
!  =====================================================================
!
!     .. parameters ..
real(dp) :: zero, half, one
parameter          ( zero = 0.0d+0, half = 0.5d+0, one = 1.0d+0 )
!     ..
!     .. local scalars ..
logical :: notran, nounit, upper
integer :: i, imax, j, jfirst, jinc, jlast
real(dp) :: bignum, grow, rec, smlnum, sumj, tjj, tjjs, tmax, tscal, uscal, xbnd, xj, xmax
!     ..
!     .. local arrays ..
real(dp) :: work(1)
!     ..
!     .. external functions ..
!     ..
!     .. external subroutines ..
!     ..
!     .. intrinsic functions ..
intrinsic          abs, max, min
!     ..
!     .. executable statements ..
!
info = 0
upper = lsame( uplo, 'u' )
notran = lsame( trans, 'n' )
nounit = lsame( diag, 'n' )
!
!     test the input parameters.
!
if( .not.upper .and. .not.lsame( uplo, 'l' ) ) then
info = -1
else if( .not.notran .and. .not.lsame( trans, 't' ) .and. .not. lsame( trans, 'c' ) ) then
info = -2
else if( .not.nounit .and. .not.lsame( diag, 'u' ) ) then
info = -3
else if( .not.lsame( normin, 'y' ) .and. .not. lsame( normin, 'n' ) ) then
info = -4
else if( n.lt.0 ) then
info = -5
else if( lda.lt.max( 1, n ) ) then
info = -7
end if
if( info.ne.0 ) then
call xerbla( 'dlatrs', -info )
return
end if
!
!     quick return if possible
!
scale = one
if( n.eq.0 ) return
!
!     determine machine dependent parameters to control overflow.
!
smlnum = dlamch( 'safe minimum' ) / dlamch( 'precision' )
bignum = one / smlnum
!
if( lsame( normin, 'n' ) ) then
!
!        compute the 1-norm of each column, not including the diagonal.
!
if( upper ) then
!
!           a is upper triangular.
!
do 10 j = 1, n
cnorm( j ) = dasum( j-1, a( 1, j ), 1 )
10 continue
else
!
!           a is lower triangular.
!
do 20 j = 1, n - 1
cnorm( j ) = dasum( n-j, a( j+1, j ), 1 )
20 continue
cnorm( n ) = zero
end if
end if
!
!     scale the column norms by tscal if the maximum element in cnorm is
!     greater than bignum.
!
imax = idamax( n, cnorm, 1 )
tmax = cnorm( imax )
if( tmax.le.bignum ) then
tscal = one
else
!
!        avoid nan generation if entries in cnorm exceed the
!        overflow threshold
!
if( tmax.le.dlamch('overflow') ) then
!           case 1: all entries in cnorm are valid floating-point numbers
tscal = one / ( smlnum*tmax )
call dscal( n, tscal, cnorm, 1 )
else
!           case 2: at least one column norm of a cannot be represented
!           as floating-point number. find the offdiagonal entry a( i, j )
!           with the largest absolute value. if this entry is not +/- infinity,
!           use this value as tscal.
tmax = zero
if( upper ) then
!
!              a is upper triangular.
!
do j = 2, n
tmax = max( dlange( 'm', j-1, 1, a( 1, j ), 1, work ), tmax )
end do
else
!
!              a is lower triangular.
!
do j = 1, n - 1
tmax = max( dlange( 'm', n-j, 1, a( j+1, j ), 1, work ), tmax )
end do
end if
!
if( tmax.le.dlamch('overflow') ) then
tscal = one / ( smlnum*tmax )
do j = 1, n
if( cnorm( j ).le.dlamch('overflow') ) then
cnorm( j ) = cnorm( j )*tscal
else
!                    recompute the 1-norm without introducing infinity
!                    in the summation
cnorm( j ) = zero
if( upper ) then
do i = 1, j - 1
cnorm( j ) = cnorm( j ) + tscal * abs( a( i, j ) )
end do
else
do i = j + 1, n
cnorm( j ) = cnorm( j ) + tscal * abs( a( i, j ) )
end do
end if
end if
end do
else
!              at least one entry of a is not a valid floating-point entry.
!              rely on trsv to propagate inf and nan.
call dtrsv( uplo, trans, diag, n, a, lda, x, 1 )
return
end if
end if
end if
!
!     compute a bound on the computed solution vector to see if the
!     level 2 blas routine dtrsv can be used.
!
j = idamax( n, x, 1 )
xmax = abs( x( j ) )
xbnd = xmax
if( notran ) then
!
!        compute the growth in a * x = b.
!
if( upper ) then
jfirst = n
jlast = 1
jinc = -1
else
jfirst = 1
jlast = n
jinc = 1
end if
!
if( tscal.ne.one ) then
grow = zero
go to 50
end if
!
if( nounit ) then
!
!           a is non-unit triangular.
!
!           compute grow = 1/g(j) and xbnd = 1/m(j).
!           initially, g(0) = max{x(i), i=1,...,n}.
!
grow = one / max( xbnd, smlnum )
xbnd = grow
do 30 j = jfirst, jlast, jinc
!
!              exit the loop if the growth factor is too small.
!
if( grow.le.smlnum ) go to 50
!
!              m(j) = g(j-1) / abs(a(j,j))
!
tjj = abs( a( j, j ) )
xbnd = min( xbnd, min( one, tjj )*grow )
if( tjj+cnorm( j ).ge.smlnum ) then
!
!                 g(j) = g(j-1)*( 1 + cnorm(j) / abs(a(j,j)) )
!
grow = grow*( tjj / ( tjj+cnorm( j ) ) )
else
!
!                 g(j) could overflow, set grow to 0.
!
grow = zero
end if
30 continue
grow = xbnd
else
!
!           a is unit triangular.
!
!           compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
!
grow = min( one, one / max( xbnd, smlnum ) )
do 40 j = jfirst, jlast, jinc
!
!              exit the loop if the growth factor is too small.
!
if( grow.le.smlnum ) go to 50
!
!              g(j) = g(j-1)*( 1 + cnorm(j) )
!
grow = grow*( one / ( one+cnorm( j ) ) )
40 continue
end if
50 continue
!
else
!
!        compute the growth in a**t * x = b.
!
if( upper ) then
jfirst = 1
jlast = n
jinc = 1
else
jfirst = n
jlast = 1
jinc = -1
end if
!
if( tscal.ne.one ) then
grow = zero
go to 80
end if
!
if( nounit ) then
!
!           a is non-unit triangular.
!
!           compute grow = 1/g(j) and xbnd = 1/m(j).
!           initially, m(0) = max{x(i), i=1,...,n}.
!
grow = one / max( xbnd, smlnum )
xbnd = grow
do 60 j = jfirst, jlast, jinc
!
!              exit the loop if the growth factor is too small.
!
if( grow.le.smlnum ) go to 80
!
!              g(j) = max( g(j-1), m(j-1)*( 1 + cnorm(j) ) )
!
xj = one + cnorm( j )
grow = min( grow, xbnd / xj )
!
!              m(j) = m(j-1)*( 1 + cnorm(j) ) / abs(a(j,j))
!
tjj = abs( a( j, j ) )
if( xj.gt.tjj ) xbnd = xbnd*( tjj / xj )
60 continue
grow = min( grow, xbnd )
else
!
!           a is unit triangular.
!
!           compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
!
grow = min( one, one / max( xbnd, smlnum ) )
do 70 j = jfirst, jlast, jinc
!
!              exit the loop if the growth factor is too small.
!
if( grow.le.smlnum ) go to 80
!
!              g(j) = ( 1 + cnorm(j) )*g(j-1)
!
xj = one + cnorm( j )
grow = grow / xj
70 continue
end if
80 continue
end if
!
if( ( grow*tscal ).gt.smlnum ) then
!
!        use the level 2 blas solve if the reciprocal of the bound on
!        elements of x is not too small.
!
call dtrsv( uplo, trans, diag, n, a, lda, x, 1 )
else
!
!        use a level 1 blas solve, scaling intermediate results.
!
if( xmax.gt.bignum ) then
!
!           scale x so that its components are less than or equal to
!           bignum in absolute value.
!
scale = bignum / xmax
call dscal( n, scale, x, 1 )
xmax = bignum
end if
!
if( notran ) then
!
!           solve a * x = b
!
do 110 j = jfirst, jlast, jinc
!
!              compute x(j) = b(j) / a(j,j), scaling x if necessary.
!
xj = abs( x( j ) )
if( nounit ) then
tjjs = a( j, j )*tscal
else
tjjs = tscal
if( tscal.eq.one ) go to 100
end if
tjj = abs( tjjs )
if( tjj.gt.smlnum ) then
!
!                    abs(a(j,j)) > smlnum:
!
if( tjj.lt.one ) then
if( xj.gt.tjj*bignum ) then
!
!                          scale x by 1/b(j).
!
rec = one / xj
call dscal( n, rec, x, 1 )
scale = scale*rec
xmax = xmax*rec
end if
end if
x( j ) = x( j ) / tjjs
xj = abs( x( j ) )
else if( tjj.gt.zero ) then
!
!                    0 < abs(a(j,j)) <= smlnum:
!
if( xj.gt.tjj*bignum ) then
!
!                       scale x by (1/abs(x(j)))*abs(a(j,j))*bignum
!                       to avoid overflow when dividing by a(j,j).
!
rec = ( tjj*bignum ) / xj
if( cnorm( j ).gt.one ) then
!
!                          scale by 1/cnorm(j) to avoid overflow when
!                          multiplying x(j) times column j.
!
rec = rec / cnorm( j )
end if
call dscal( n, rec, x, 1 )
scale = scale*rec
xmax = xmax*rec
end if
x( j ) = x( j ) / tjjs
xj = abs( x( j ) )
else
!
!                    a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
!                    scale = 0, and compute a solution to a*x = 0.
!
do 90 i = 1, n
x( i ) = zero
90 continue
x( j ) = one
xj = one
scale = zero
xmax = zero
end if
100 continue
!
!              scale x if necessary to avoid overflow when adding a
!              multiple of column j of a.
!
if( xj.gt.one ) then
rec = one / xj
if( cnorm( j ).gt.( bignum-xmax )*rec ) then
!
!                    scale x by 1/(2*abs(x(j))).
!
rec = rec*half
call dscal( n, rec, x, 1 )
scale = scale*rec
end if
else if( xj*cnorm( j ).gt.( bignum-xmax ) ) then
!
!                 scale x by 1/2.
!
call dscal( n, half, x, 1 )
scale = scale*half
end if
!
if( upper ) then
if( j.gt.1 ) then
!
!                    compute the update
!                       x(1:j-1) := x(1:j-1) - x(j) * a(1:j-1,j)
!
call daxpy( j-1, -x( j )*tscal, a( 1, j ), 1, x, 1 )
i = idamax( j-1, x, 1 )
xmax = abs( x( i ) )
end if
else
if( j.lt.n ) then
!
!                    compute the update
!                       x(j+1:n) := x(j+1:n) - x(j) * a(j+1:n,j)
!
call daxpy( n-j, -x( j )*tscal, a( j+1, j ), 1, x( j+1 ), 1 )
i = j + idamax( n-j, x( j+1 ), 1 )
xmax = abs( x( i ) )
end if
end if
110 continue
!
else
!
!           solve a**t * x = b
!
do 160 j = jfirst, jlast, jinc
!
!              compute x(j) = b(j) - sum a(k,j)*x(k).
!                                    k<>j
!
xj = abs( x( j ) )
uscal = tscal
rec = one / max( xmax, one )
if( cnorm( j ).gt.( bignum-xj )*rec ) then
!
!                 if x(j) could overflow, scale x by 1/(2*xmax).
!
rec = rec*half
if( nounit ) then
tjjs = a( j, j )*tscal
else
tjjs = tscal
end if
tjj = abs( tjjs )
if( tjj.gt.one ) then
!
!                       divide by a(j,j) when scaling x if a(j,j) > 1.
!
rec = min( one, rec*tjj )
uscal = uscal / tjjs
end if
if( rec.lt.one ) then
call dscal( n, rec, x, 1 )
scale = scale*rec
xmax = xmax*rec
end if
end if
!
sumj = zero
if( uscal.eq.one ) then
!
!                 if the scaling needed for a in the dot product is 1,
!                 call ddot to perform the dot product.
!
if( upper ) then
sumj = ddot( j-1, a( 1, j ), 1, x, 1 )
else if( j.lt.n ) then
sumj = ddot( n-j, a( j+1, j ), 1, x( j+1 ), 1 )
end if
else
!
!                 otherwise, use in-line code for the dot product.
!
if( upper ) then
do 120 i = 1, j - 1
sumj = sumj + ( a( i, j )*uscal )*x( i )
120 continue
else if( j.lt.n ) then
do 130 i = j + 1, n
sumj = sumj + ( a( i, j )*uscal )*x( i )
130 continue
end if
end if
!
if( uscal.eq.tscal ) then
!
!                 compute x(j) := ( x(j) - sumj ) / a(j,j) if 1/a(j,j)
!                 was not used to scale the dotproduct.
!
x( j ) = x( j ) - sumj
xj = abs( x( j ) )
if( nounit ) then
tjjs = a( j, j )*tscal
else
tjjs = tscal
if( tscal.eq.one ) go to 150
end if
!
!                    compute x(j) = x(j) / a(j,j), scaling if necessary.
!
tjj = abs( tjjs )
if( tjj.gt.smlnum ) then
!
!                       abs(a(j,j)) > smlnum:
!
if( tjj.lt.one ) then
if( xj.gt.tjj*bignum ) then
!
!                             scale x by 1/abs(x(j)).
!
rec = one / xj
call dscal( n, rec, x, 1 )
scale = scale*rec
xmax = xmax*rec
end if
end if
x( j ) = x( j ) / tjjs
else if( tjj.gt.zero ) then
!
!                       0 < abs(a(j,j)) <= smlnum:
!
if( xj.gt.tjj*bignum ) then
!
!                          scale x by (1/abs(x(j)))*abs(a(j,j))*bignum.
!
rec = ( tjj*bignum ) / xj
call dscal( n, rec, x, 1 )
scale = scale*rec
xmax = xmax*rec
end if
x( j ) = x( j ) / tjjs
else
!
!                       a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
!                       scale = 0, and compute a solution to a**t*x = 0.
!
do 140 i = 1, n
x( i ) = zero
140 continue
x( j ) = one
scale = zero
xmax = zero
end if
150 continue
else
!
!                 compute x(j) := x(j) / a(j,j)  - sumj if the dot
!                 product has already been divided by 1/a(j,j).
!
x( j ) = x( j ) / tjjs - sumj
end if
xmax = max( xmax, abs( x( j ) ) )
160 continue
end if
scale = scale / tscal
end if
!
!     scale the column norms by 1/tscal for return.
!
if( tscal.ne.one ) then
call dscal( n, one / tscal, cnorm, 1 )
end if
!
return
!
!     end of dlatrs
!
end subroutine dlatrs
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/drscl.f
!> \brief \b drscl multiplies a vector by the reciprocal of a real scalar.
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download drscl + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/drscl.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/drscl.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/drscl.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       subroutine drscl( n, sa, sx, incx )
!
!       .. scalar arguments ..
!       integer            incx, n
!       double precision   sa
!       ..
!       .. array arguments ..
!       double precision   sx( * )
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> drscl multiplies an n-element real vector x by the real scalar 1/a.
!> this is done without overflow or underflow as long as
!> the final result x/a does not overflow or underflow.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] n
!> \verbatim
!>          n is integer
!>          the number of components of the vector x.
!> \endverbatim
!>
!> \param[in] sa
!> \verbatim
!>          sa is double precision
!>          the scalar a which is used to divide each component of x.
!>          sa must be >= 0, or the subroutine will divide by zero.
!> \endverbatim
!>
!> \param[in,out] sx
!> \verbatim
!>          sx is double precision array, dimension
!>                         (1+(n-1)*abs(incx))
!>          the n-element vector x.
!> \endverbatim
!>
!> \param[in] incx
!> \verbatim
!>          incx is integer
!>          the increment between successive values of the vector sx.
!>          > 0:  sx(1) = x(1) and sx(1+(i-1)*incx) = x(i),     1< i<= n
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup rscl
!
!  =====================================================================
subroutine drscl( n, sa, sx, incx )
implicit none
!
!  -- lapack auxiliary routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
integer, intent(in) :: incx !! the increment between successive values of the vector sx. > 0:
!! sx(1) = x(1) and sx(1+(i-1)*incx) = x(i),     1< i<= n
integer, intent(in) :: n !! the number of components of the vector x.
real(dp), intent(in) :: sa !! the scalar a which is used to divide each component of x. sa must be
!! >= 0, or the subroutine will divide by zero.
!     ..
!     .. array arguments ..
real(dp), intent(inout) :: sx( * ) !! (1+(n-1)*abs(incx)) the n-element vector x.
!     ..
!
! =====================================================================
!
!     .. parameters ..
real(dp) :: one, zero
parameter          ( one = 1.0d+0, zero = 0.0d+0 )
!     ..
!     .. local scalars ..
logical :: done
real(dp) :: bignum, cden, cden1, cnum, cnum1, mul, smlnum
!     ..
!     .. external functions ..
!     ..
!     .. external subroutines ..
!     ..
!     .. intrinsic functions ..
intrinsic          abs
!     ..
!     .. executable statements ..
!
!     quick return if possible
!
if( n.le.0 ) return
!
!     get machine parameters
!
smlnum = dlamch( 's' )
bignum = one / smlnum
!
!     initialize the denominator to sa and the numerator to 1.
!
cden = sa
cnum = one
!
10 continue
cden1 = cden*smlnum
cnum1 = cnum / bignum
if( abs( cden1 ).gt.abs( cnum ) .and. cnum.ne.zero ) then
!
!        pre-multiply x by smlnum if cden is large compared to cnum.
!
mul = smlnum
done = .false.
cden = cden1
else if( abs( cnum1 ).gt.abs( cden ) ) then
!
!        pre-multiply x by bignum if cden is small compared to cnum.
!
mul = bignum
done = .false.
cnum = cnum1
else
!
!        multiply x by cnum / cden and return.
!
mul = cnum / cden
done = .true.
end if
!
!     scale the vector x by mul
!
call dscal( n, mul, sx, incx )
!
if( .not.done ) go to 10
!
return
!
!     end of drscl
!
end subroutine drscl
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlange.f
!> \brief \b dlange returns the value of the 1-norm, frobenius norm, infinity-norm, or the
!largest absolute value of any element of a general rectangular matrix.
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> download dlange + dependencies
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlange.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlange.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-
!bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlange.f">
!> [txt]</a>
!> \endhtmlonly
!
!  definition:
!  ===========
!
!       double precision function dlange( norm, m, n, a, lda, work )
!
!       .. scalar arguments ..
!       character          norm
!       integer            lda, m, n
!       ..
!       .. array arguments ..
!       double precision   a( lda, * ), work( * )
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dlange  returns the value of the one norm,  or the frobenius norm, or
!> the  infinity norm,  or the  element of  largest absolute value  of a
!> real matrix a.
!> \endverbatim
!>
!> \return dlange
!> \verbatim
!>
!>    dlange = ( max(abs(a(i,j))), norm = 'm' or 'm'
!>             (
!>             ( norm1(a),         norm = '1', 'o' or 'o'
!>             (
!>             ( normi(a),         norm = 'i' or 'i'
!>             (
!>             ( normf(a),         norm = 'f', 'f', 'e' or 'e'
!>
!> where  norm1  denotes the  one norm of a matrix (maximum column sum),
!> normi  denotes the  infinity norm  of a matrix  (maximum row sum) and
!> normf  denotes the  frobenius norm of a matrix (square root of sum of
!> squares).  note that  max(abs(a(i,j)))  is not a consistent matrix norm.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] norm
!> \verbatim
!>          norm is character*1
!>          specifies the value to be returned in dlange as described
!>          above.
!> \endverbatim
!>
!> \param[in] m
!> \verbatim
!>          m is integer
!>          the number of rows of the matrix a.  m >= 0.  when m = 0,
!>          dlange is set to zero.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is integer
!>          the number of columns of the matrix a.  n >= 0.  when n = 0,
!>          dlange is set to zero.
!> \endverbatim
!>
!> \param[in] a
!> \verbatim
!>          a is double precision array, dimension (lda,n)
!>          the m by n matrix a.
!> \endverbatim
!>
!> \param[in] lda
!> \verbatim
!>          lda is integer
!>          the leading dimension of the array a.  lda >= max(m,1).
!> \endverbatim
!>
!> \param[out] work
!> \verbatim
!>          work is double precision array, dimension (max(1,lwork)),
!>          where lwork >= m when norm = 'i'; otherwise, work is not
!>          referenced.
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup lange
!
!  =====================================================================
real(dp) function dlange( norm, m, n, a, lda, work ) result(lapack_result)
implicit none
!
!  -- lapack auxiliary routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
character, intent(in) :: norm !! specifies the value to be returned in dlange as described above.
integer, intent(in) :: lda !! the leading dimension of the array a.  lda >= max(m,1).
integer, intent(in) :: m !! the number of rows of the matrix a.  m >= 0.  when m = 0, dlange is
!! set to zero.
integer, intent(in) :: n !! the number of columns of the matrix a.  n >= 0.  when n = 0, dlange is
!! set to zero.
!     ..
!     .. array arguments ..
real(dp), intent(in) :: a( lda, * ) !! the m by n matrix a.
real(dp), intent(out) :: work( * ) !! where lwork >= m when norm = 'i'; otherwise, work is not
!! referenced.
!     ..
!
! =====================================================================
!
!     .. parameters ..
real(dp) :: one, zero
parameter          ( one = 1.0d+0, zero = 0.0d+0 )
!     ..
!     .. local scalars ..
integer :: i, j
real(dp) :: scale, sum, value, temp
!     ..
!     .. external subroutines ..
!     ..
!     .. external functions ..
!     ..
!     .. intrinsic functions ..
intrinsic          abs, min, sqrt
!     ..
!     .. executable statements ..
!
if( min( m, n ).eq.0 ) then
value = zero
else if( lsame( norm, 'm' ) ) then
!
!        find max(abs(a(i,j))).
!
value = zero
do 20 j = 1, n
do 10 i = 1, m
temp = abs( a( i, j ) )
if( value.lt.temp .or. disnan( temp ) ) value = temp
10 continue
20 continue
else if( ( lsame( norm, 'o' ) ) .or. ( norm.eq.'1' ) ) then
!
!        find norm1(a).
!
value = zero
do 40 j = 1, n
sum = zero
do 30 i = 1, m
sum = sum + abs( a( i, j ) )
30 continue
if( value.lt.sum .or. disnan( sum ) ) value = sum
40 continue
else if( lsame( norm, 'i' ) ) then
!
!        find normi(a).
!
do 50 i = 1, m
work( i ) = zero
50 continue
do 70 j = 1, n
do 60 i = 1, m
work( i ) = work( i ) + abs( a( i, j ) )
60 continue
70 continue
value = zero
do 80 i = 1, m
temp = work( i )
if( value.lt.temp .or. disnan( temp ) ) value = temp
80 continue
else if( ( lsame( norm, 'f' ) ) .or. ( lsame( norm, 'e' ) ) ) then
!
!        find normf(a).
!
scale = zero
sum = one
do 90 j = 1, n
call dlassq( m, a( 1, j ), 1, scale, sum )
90 continue
value = scale*sqrt( sum )
end if
!
lapack_result = value
return
!
!     end of dlange
!
end function dlange
! Upstream: https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/BLAS/SRC/dtrsv.f
!> \brief \b dtrsv
!
!  =========== documentation ===========
!
! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  definition:
!  ===========
!
!       subroutine dtrsv(uplo,trans,diag,n,a,lda,x,incx)
!
!       .. scalar arguments ..
!       integer incx,lda,n
!       character diag,trans,uplo
!       ..
!       .. array arguments ..
!       double precision a(lda,*),x(*)
!       ..
!
!
!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dtrsv  solves one of the systems of equations
!>
!>    a*x = b,   or   a**t*x = b,
!>
!> where b and x are n element vectors and a is an n by n unit, or
!> non-unit, upper or lower triangular matrix.
!>
!> no test for singularity or near-singularity is included in this
!> routine. such tests must be performed before calling this routine.
!> \endverbatim
!
!  arguments:
!  ==========
!
!> \param[in] uplo
!> \verbatim
!>          uplo is character*1
!>           on entry, uplo specifies whether the matrix is an upper or
!>           lower triangular matrix as follows:
!>
!>              uplo = 'u' or 'u'   a is an upper triangular matrix.
!>
!>              uplo = 'l' or 'l'   a is a lower triangular matrix.
!> \endverbatim
!>
!> \param[in] trans
!> \verbatim
!>          trans is character*1
!>           on entry, trans specifies the equations to be solved as
!>           follows:
!>
!>              trans = 'n' or 'n'   a*x = b.
!>
!>              trans = 't' or 't'   a**t*x = b.
!>
!>              trans = 'c' or 'c'   a**t*x = b.
!> \endverbatim
!>
!> \param[in] diag
!> \verbatim
!>          diag is character*1
!>           on entry, diag specifies whether or not a is unit
!>           triangular as follows:
!>
!>              diag = 'u' or 'u'   a is assumed to be unit triangular.
!>
!>              diag = 'n' or 'n'   a is not assumed to be unit
!>                                  triangular.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is integer
!>           on entry, n specifies the order of the matrix a.
!>           n must be at least zero.
!> \endverbatim
!>
!> \param[in] a
!> \verbatim
!>          a is double precision array, dimension ( lda, n )
!>           before entry with  uplo = 'u' or 'u', the leading n by n
!>           upper triangular part of the array a must contain the upper
!>           triangular matrix and the strictly lower triangular part of
!>           a is not referenced.
!>           before entry with uplo = 'l' or 'l', the leading n by n
!>           lower triangular part of the array a must contain the lower
!>           triangular matrix and the strictly upper triangular part of
!>           a is not referenced.
!>           note that when  diag = 'u' or 'u', the diagonal elements of
!>           a are not referenced either, but are assumed to be unity.
!> \endverbatim
!>
!> \param[in] lda
!> \verbatim
!>          lda is integer
!>           on entry, lda specifies the first dimension of a as declared
!>           in the calling (sub) program. lda must be at least
!>           max( 1, n ).
!> \endverbatim
!>
!> \param[in,out] x
!> \verbatim
!>          x is double precision array, dimension at least
!>           ( 1 + ( n - 1 )*abs( incx ) ).
!>           before entry, the incremented array x must contain the n
!>           element right-hand side vector b. on exit, x is overwritten
!>           with the solution vector x.
!> \endverbatim
!>
!> \param[in] incx
!> \verbatim
!>          incx is integer
!>           on entry, incx specifies the increment for the elements of
!>           x. incx must not be zero.
!>
!>  level 2 blas routine.
!>
!>  -- written on 22-october-1986.
!>     jack dongarra, argonne national lab.
!>     jeremy du croz, nag central office.
!>     sven hammarling, nag central office.
!>     richard hanson, sandia national labs.
!> \endverbatim
!
!  authors:
!  ========
!
!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.
!
!> \ingroup trsv
!
!  =====================================================================
subroutine dtrsv(uplo,trans,diag,n,a,lda,x,incx)
implicit none
!
!  -- reference blas level2 routine --
!  -- reference blas is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
!     .. scalar arguments ..
integer, intent(in) :: incx !! on entry, incx specifies the increment for the elements of x. incx
!! must not be zero. level 2 blas routine.
integer, intent(in) :: lda !! on entry, lda specifies the first dimension of a as declared in the
!! calling (sub) program. lda must be at least
integer, intent(in) :: n !! on entry, n specifies the order of the matrix a. n must be at least
!! zero.
character, intent(in) :: diag !! on entry, diag specifies whether or not a is unit triangular as
!! follows: diag = 'u' or 'u'   a is assumed to be unit triangular.
character, intent(in) :: trans !! on entry, trans specifies the equations to be solved as follows:
!! trans = 'n' or 'n'   a*x = b.
character, intent(in) :: uplo !! on entry, uplo specifies whether the matrix is an upper or lower
!! triangular matrix as follows:
!     ..
!     .. array arguments ..
real(dp), intent(in) :: a(lda,*) !! before entry with  uplo = 'u' or 'u', the leading n by n upper
!! triangular part of the array a must contain the upper
real(dp), intent(inout) :: x(*) !! ( 1 + ( n - 1 )*abs( incx ) ). before entry, the incremented
!! array x must contain the n element right-hand side vector b. on
!! exit, x is overwritten
!     ..
!
!  =====================================================================
!
!     .. parameters ..
real(dp) :: zero
parameter (zero=0.0d+0)
!     ..
!     .. local scalars ..
real(dp) :: temp
integer :: i, info, ix, j, jx, kx
logical :: nounit
!     ..
!     .. external functions ..
!     ..
!     .. external subroutines ..
!     ..
!     .. intrinsic functions ..
intrinsic max
!     ..
!
!     test the input parameters.
!
info = 0
if (.not.lsame(uplo,'u') .and. .not.lsame(uplo,'l')) then
info = 1
else if (.not.lsame(trans,'n') .and. .not.lsame(trans,'t') .and. .not.lsame(trans,'c')) then
info = 2
else if (.not.lsame(diag,'u') .and. .not.lsame(diag,'n')) then
info = 3
else if (n.lt.0) then
info = 4
else if (lda.lt.max(1,n)) then
info = 6
else if (incx.eq.0) then
info = 8
end if
if (info.ne.0) then
call xerbla('dtrsv ',info)
return
end if
!
!     quick return if possible.
!
if (n.eq.0) return
!
nounit = lsame(diag,'n')
!
!     set up the start point in x if the increment is not unity. this
!     will be  ( n - 1 )*incx  too small for descending loops.
!
if (incx.le.0) then
kx = 1 - (n-1)*incx
else if (incx.ne.1) then
kx = 1
end if
!
!     start the operations. in this version the elements of a are
!     accessed sequentially with one pass through a.
!
if (lsame(trans,'n')) then
!
!        form  x := inv( a )*x.
!
if (lsame(uplo,'u')) then
if (incx.eq.1) then
do 20 j = n,1,-1
if (x(j).ne.zero) then
if (nounit) x(j) = x(j)/a(j,j)
temp = x(j)
do 10 i = j - 1,1,-1
x(i) = x(i) - temp*a(i,j)
10 continue
end if
20 continue
else
jx = kx + (n-1)*incx
do 40 j = n,1,-1
if (x(jx).ne.zero) then
if (nounit) x(jx) = x(jx)/a(j,j)
temp = x(jx)
ix = jx
do 30 i = j - 1,1,-1
ix = ix - incx
x(ix) = x(ix) - temp*a(i,j)
30 continue
end if
jx = jx - incx
40 continue
end if
else
if (incx.eq.1) then
do 60 j = 1,n
if (x(j).ne.zero) then
if (nounit) x(j) = x(j)/a(j,j)
temp = x(j)
do 50 i = j + 1,n
x(i) = x(i) - temp*a(i,j)
50 continue
end if
60 continue
else
jx = kx
do 80 j = 1,n
if (x(jx).ne.zero) then
if (nounit) x(jx) = x(jx)/a(j,j)
temp = x(jx)
ix = jx
do 70 i = j + 1,n
ix = ix + incx
x(ix) = x(ix) - temp*a(i,j)
70 continue
end if
jx = jx + incx
80 continue
end if
end if
else
!
!        form  x := inv( a**t )*x.
!
if (lsame(uplo,'u')) then
if (incx.eq.1) then
do 100 j = 1,n
temp = x(j)
do 90 i = 1,j - 1
temp = temp - a(i,j)*x(i)
90 continue
if (nounit) temp = temp/a(j,j)
x(j) = temp
100 continue
else
jx = kx
do 120 j = 1,n
temp = x(jx)
ix = kx
do 110 i = 1,j - 1
temp = temp - a(i,j)*x(ix)
ix = ix + incx
110 continue
if (nounit) temp = temp/a(j,j)
x(jx) = temp
jx = jx + incx
120 continue
end if
else
if (incx.eq.1) then
do 140 j = n,1,-1
temp = x(j)
do 130 i = n,j + 1,-1
temp = temp - a(i,j)*x(i)
130 continue
if (nounit) temp = temp/a(j,j)
x(j) = temp
140 continue
else
kx = kx + (n-1)*incx
jx = kx
do 160 j = n,1,-1
temp = x(jx)
ix = kx
do 150 i = n,j + 1,-1
temp = temp - a(i,j)*x(ix)
ix = ix - incx
150 continue
if (nounit) temp = temp/a(j,j)
x(jx) = temp
jx = jx - incx
160 continue
end if
end if
end if
!
return
!
!     end of dtrsv
!
end subroutine dtrsv
logical function lsame(a,b) result(same)
implicit none
character, intent(in) :: a !! First option letter, ASCII case insensitive.
character, intent(in) :: b !! Second option letter, ASCII case insensitive.
integer :: ia,ib
ia=iachar(a)
ib=iachar(b)
if (ia>=65 .and. ia<=90) ia=ia+32
if (ib>=65 .and. ib<=90) ib=ib+32
same=ia==ib
end function lsame
subroutine xerbla(name,index)
implicit none
character(*), intent(in) :: name !! Internal routine reporting an invalid argument.
integer, intent(in) :: index !! Invalid argument position; caller already sets INFO.
! All entry points are private and validated by the SciFort wrapper.
! LAPACK callers set INFO and return after this deliberately silent hook.
end subroutine xerbla
end module scifort_linalg_qz_lapack
