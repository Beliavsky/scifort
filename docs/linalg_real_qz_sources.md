# Real QZ source catalog, 2026-11-15

The private `src/linalg_qz_lapack.f90` adapts the following exact
Reference-LAPACK **v3.12.1** sources. Dependencies are kept private and
self-contained; they include the real block-swap Sylvester solves and their
conditioning support paths. No LAPACK/BLAS DLL or external dependency is
required. The complete upstream license remains in
[`LICENSES/lapack-BSD-3-Clause.txt`](../LICENSES/lapack-BSD-3-Clause.txt).

License SHA-256: `978539245c405775bf347b1256c36df4daf8b69a64b099e94b1d53742201621e`.

| Exact upstream source | SHA-256 |
| --- | --- |
| [SRC/dhgeqz.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dhgeqz.f) | `263c768401d3a2eaf472d34ae2097a168e252a998b612879d6627a3fe4f7af14` |
| [SRC/dgghrd.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dgghrd.f) | `3e80ba8d0c9435c0d7e96296d8c5623a0ce56f32c94601eb7cc06cf25e69ba54` |
| [SRC/dtgexc.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dtgexc.f) | `3d626ac2577b8ac8b14ccceabd0f256ef3cf9c7ac992a0dedb5e4aca49c20c71` |
| [INSTALL/dlamch.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/INSTALL/dlamch.f) | `d532bb8ee8e3b9ffe29bf91eba98f4a304dfb599009dc105279ba1e028b96db0` |
| [SRC/dlanhs.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlanhs.f) | `68531c728815c691fa06f6f7e6611764895d47e6279da872f459e87d306e2102` |
| [SRC/dlapy2.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlapy2.f) | `274638b08ffe1c6f7bde140f6c9bf8bfab9087731c761f35a0b668f6d5251565` |
| [SRC/dlapy3.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlapy3.f) | `468594c6a79215e5379699fd849bfbeff7915404ee68880b8f8b97c6b4a0fa10` |
| [SRC/dlag2.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlag2.f) | `909b9f549b7c0af6b350c7caffd6e5000b41ad6e7a968e95b0513d7e050b1a2f` |
| [SRC/dlarfg.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlarfg.f) | `2e6f99fe4c695860699362b5206f205952d598ffcb0e6ef62eb0f31954e311cb` |
| [SRC/dlartg.f90](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlartg.f90) | `c5fbb9c2035f0314933b47e4247dfd8f1944f671d7764be2cc8fb15439e4be0c` |
| [SRC/dlaset.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlaset.f) | `85d603b4ace0f14d757548eb1f2eff8d51b4111b81a2a629c31606d4a8b7db73` |
| [SRC/dlasv2.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlasv2.f) | `a388627977f9c7d79eddb235ffa33e69b407ed48080c5d388c84ceba93507dd3` |
| [BLAS/SRC/drot.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/BLAS/SRC/drot.f) | `973468ddcb08bd3ec04aed00fafae46915a7cce71ad17925516e29a9ea281fd9` |
| [SRC/dtgex2.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dtgex2.f) | `58a50c2ce53f333a58b112bf1d79135f3c0b3a09ee97c420296fab8db97acf2a` |
| [SRC/dlassq.f90](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlassq.f90) | `84cd1d0035b3b82374f8a591fce1a2a4a7ef7a11e57cae15a84a22de5738a3d5` |
| [SRC/disnan.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/disnan.f) | `781881a6f42e246bcf4d4fd588cf46a784c476e697659ae9019458097d15c897` |
| [BLAS/SRC/dnrm2.f90](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/BLAS/SRC/dnrm2.f90) | `2a5bc5ec9ab8077da1f81e2169a3868df1629e39f36ba969ada1f491a3b94979` |
| [BLAS/SRC/dscal.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/BLAS/SRC/dscal.f) | `ea3794ad26c4f1dc852a52c0aee10edb25b0e7f0d12a61c91656144f139c221c` |
| [BLAS/SRC/dgemm.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/BLAS/SRC/dgemm.f) | `aa22f67c2c3ef9b0b0cf42fcd6c4080076948d5f4b8b4d982a78eb6e33f9e191` |
| [SRC/dgeqr2.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dgeqr2.f) | `49c9dcb928683767041f1434b308e5b7b27bbe84f4d573cb540140f4b8e234f8` |
| [SRC/dgerq2.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dgerq2.f) | `615934f2d61c7f4b6b24a0c351caefb4d81bd7fc9935ca00a8cad8e8871b4092` |
| [SRC/dlacpy.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlacpy.f) | `00a07a6244d34435e247cb4b037dcf512d10c27177c8f96410d325e9554c7c18` |
| [SRC/dlagv2.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlagv2.f) | `23877b077f5f0452bff68eb8e152531d87616822e76933bc8f01504b03fd5546` |
| [SRC/dorg2r.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dorg2r.f) | `15b212557a43409d7f31a55473052054f0f5da9840ee3f7c32c90ded4072ede2` |
| [SRC/dorgr2.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dorgr2.f) | `6ea6c50aa55a96044495775dda473d44a81e7862d83169a9caab23ed92d9a51b` |
| [SRC/dorm2r.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dorm2r.f) | `aa710e75be766dd3500b971808568ec9af3f32cbb0333bceb09c47d3975a83bb` |
| [SRC/dormr2.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dormr2.f) | `61e89c43d48e932f49ce88ab4bcb84bfbd97886e28f5ff4548ebad6ced2a1965` |
| [SRC/dtgsy2.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dtgsy2.f) | `feab4f13f7020cc31825aeddbbdef94bfad7be70ad088d1819f4d96ccd43843e` |
| [SRC/dlaisnan.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlaisnan.f) | `2d704d97de2f31446e8615602e5e4d1d5fd7a688695b4c1775ffd9cedd0f230c` |
| [SRC/dlarf1f.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlarf1f.f) | `a09860405adf9e56c42d51a107abc764b51367c3864102e83c0710d5cec48738` |
| [SRC/dlarf1l.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlarf1l.f) | `aa82686d354f1114acf9b22d53dc04b54dc9e1d171f6269d20d02b04302fcc8c` |
| [SRC/dlarf.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlarf.f) | `c1a4990c808ff3381226863350c16e25e54fd6833bb320d3d6c335208eece011` |
| [BLAS/SRC/daxpy.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/BLAS/SRC/daxpy.f) | `6d978fe749ad22249d05a230b9c0ba07e6afccbcc59b746f198e285487811af1` |
| [BLAS/SRC/dcopy.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/BLAS/SRC/dcopy.f) | `fc821bf4413f81d371f6678eb9888622008ed44819576cbc6d13fa66df3e1672` |
| [BLAS/SRC/dgemv.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/BLAS/SRC/dgemv.f) | `2320b9024e296f2a78882dd84c1b96be0af985dc177f2f373cb122f86f5d925d` |
| [BLAS/SRC/dger.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/BLAS/SRC/dger.f) | `6da9a924d4eed5fca0ea42d34237c2ff8b75296ce17bcacfd2b9a829a9e47f4f` |
| [SRC/dgesc2.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dgesc2.f) | `f1fe851e55aab6bcdaa859808bb2bb9201b5cdc1ca064d3da505a0a54bc5326d` |
| [SRC/dgetc2.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dgetc2.f) | `206c7c49541cf3c682d250eeb622c74a9b2531cdd44a7977c7ad923e8d49e3c7` |
| [SRC/dlatdf.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlatdf.f) | `1c0ab10d26a591ffb815834db51a68ade842d5d63618a0e4a5e9437acb344f2b` |
| [SRC/iladlr.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/iladlr.f) | `0f3f0ed30ba918673589cdb963f23719903b91b6e1ac0e1e4535fd978502d003` |
| [SRC/iladlc.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/iladlc.f) | `48d06e17c656c37a50e223634d188ecab415fac39f3f3e21cecf929b1fc10fab` |
| [SRC/dlaswp.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlaswp.f) | `8f50f034f7c1c075db12dcb277e743fdd4938d788fde864a1adcca70fedbc102` |
| [BLAS/SRC/idamax.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/BLAS/SRC/idamax.f) | `05e6970ea232fceb382f20be73f222c2a83e9cd8872e930a210e7e208d5ca3e4` |
| [BLAS/SRC/dswap.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/BLAS/SRC/dswap.f) | `cb29b27f67c436340bfc6be497de96aa9ccfaf4a8ae4d7e63d7a7603fdc100ed` |
| [SRC/dgecon.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dgecon.f) | `ac53214ea31c342aab4dd2cc4f5d1d4c8937eb51805eda864d13364b3a9d1141` |
| [BLAS/SRC/dasum.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/BLAS/SRC/dasum.f) | `4a315811bfdd7870f49c540b87b1040a6fa1476ab3ad02851477c7354893c275` |
| [BLAS/SRC/ddot.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/BLAS/SRC/ddot.f) | `8983346c864ecb34ef9da011a5d141c68a7329d9bc71fceb5acf2c1949b5bb4f` |
| [SRC/dlacn2.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlacn2.f) | `3aabe9c6a0f01a9e9b41076c63b4c2cfde5f4c1c466c4863fcee2cee8a033afa` |
| [SRC/dlatrs.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlatrs.f) | `b7135a1ef9bab94f07ada5ddbfe3e2aba93e8516df9ed0de9df24217d2fd54cd` |
| [SRC/drscl.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/drscl.f) | `82a581ea1a391f35f85715a405682b506952c8d2f6c5632354281ecb725e2071` |
| [SRC/dlange.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/SRC/dlange.f) | `a438a2c4ef52e8b45177b46772879abe37ef38853c62e53bfe17d8c0548a7a96` |
| [BLAS/SRC/dtrsv.f](https://raw.githubusercontent.com/Reference-LAPACK/lapack/v3.12.1/BLAS/SRC/dtrsv.f) | `f111927caedd1f2cdde08948fdf3ed21df26b99b1d7265dd5435028a60c661f6` |
| [SRC/la_constants.f90](https://github.com/Reference-LAPACK/lapack/blob/v3.12.1/SRC/la_constants.f90) (selected binary64 safe-scaling definitions) | `1d9bdd0fad20e76f6c418b0d5405aa7a256d89cfe57eba2a36a9c301c5ce7d1b` |

Modifications: fixed form to lowercase free form; `dp` from `scifort_kinds`;
private module interfaces instead of implicit external calls; explicit dummy
intents and trailing documentation; function `result` names; safe-scaling
constants restricted to binary64; intrinsic `ieee_is_nan` replaces the
upstream `LA_ISNAN` generic. The private ASCII option comparator is original
SciFort code. The private error hook returns silently after callers set
`INFO`, so library code never prints or stops. `DORMR2` temporarily changes
and restores its reflector storage, so its private `a` dummy is correctly
declared `intent(inout)` despite the upstream descriptive input annotation.

No numerical coefficients, QZ iteration steps, exceptional shifts, block
swap acceptance tests or condition-estimation branches were changed.
All original author/reference comments, NAG attribution and the applicable
copyright notice are retained. Original adapters use MIT; the mixed source
is marked `MIT AND BSD-3-Clause`. No conversion generator or upstream archive
is committed.

The separately compiled unchanged reference also uses original `LSAME`,
`XERBLA`, `LA_CONSTANTS` and `LA_XISNAN`; they are validation-only cached source.
The full upstream LAPACK test suite was not run. See the
[validation record](../validation/real_qz_2026-11-15.md) and
[provenance](../CODE_PROVENANCE.md).
