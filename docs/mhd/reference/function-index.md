# MHD Solver Function Index

This document is a mechanically generated function index for the C MHD solver (`/workspace/src/mhd`). It was produced by a Python script (`/tmp/gen_function_index.py`, not checked into the repo) that parses top-level function definitions out of the 8 files in scope, builds a call graph with comments and string literals stripped, and determines liveness transitively from a root set (the public API in `mhd.h`, plus every name referenced from `src/kinetic/*.cpp`, `src/tasks/*.cpp`, and `tests/regression/*.c`).

It reflects git commit `7973092` (measured on branch `refactor/p0-harness`). A second agent was concurrently editing `src/mhd/*.c` and committing per file; this document was generated from the tree exactly as it stood at that commit.

**This document must be regenerated after any change to the source in `src/mhd`.** Line ranges, lengths, and LIVE/DEAD status will drift the moment a function is added, removed, reordered, or gains/loses a caller.

**Line numbers are PHYSICAL line numbers** (what an editor shows, i.e. `sed -n`/`wc -l` line counting), not compiler/debugger line numbers. Several of these files carry `#line` directives (used to remap diagnostics to a different numbering, e.g. `mass_matrix_coefficients.c` has 26 of them, `ts_functions.c` has 72), so a line number reported by a compiler warning or `gdb` for one of these files will NOT match the line numbers in this document.

**Method notes.** A function is LIVE if transitively reachable from a root. Roots are: (a) every name declared in `src/mhd/mhd.h` (marked ENTRY below), and (b) every bare identifier found in `src/kinetic/*.cpp`, `src/tasks/*.cpp`, and `tests/regression/*.c` after stripping comments, that matches a function defined in these 8 files. Matching uses word-boundary identifier search, not just `name(`, so a function passed purely as a pointer (e.g. the `Fn` argument of `TSSetIFunction(ts, NULL, Fn, user)`) still counts as a call/reference. Comments (`//` and `/* */`) and string literal contents (so a `PetscLogEventRegister("Name", ...)` label is not a call) are stripped before matching, in both the `.c` files and the external root-reference files.


## src/mhd/ts_functions.c

5797 lines, 22 functions.

| Function | Lines | Length | Linkage | Status | Called from |
|---|---|---|---|---|---|
| `MFD_GetSlotsSolution` | 42-80 | 39 | static | LIVE | `src/mhd/ts_functions.c:326`, `src/mhd/ts_functions.c:1328` |
| `MFD_GetSlotsCoords` | 82-120 | 39 | static | LIVE | `src/mhd/ts_functions.c:334`, `src/mhd/ts_functions.c:1336` |
| `MFD_CellVolume` | 151-164 | 14 | static | LIVE | `src/mhd/ts_functions.c:514`, `src/mhd/ts_functions.c:1392` |
| `MFD_CellEdgeLengths` | 168-214 | 47 | static | LIVE | `src/mhd/ts_functions.c:517`, `src/mhd/ts_functions.c:846`, `src/mhd/ts_functions.c:1395` +2 more |
| `FormIJacobian_BImplicit` | 217-220 | 4 | extern | LIVE | `src/mhd/mhd.c:297` |
| `vperp_residual` | 267-1223 | 957 | static | LIVE | `src/mhd/ts_functions.c:1240`, `src/mhd/ts_functions.c:1241`, `src/mhd/ts_functions.c:1718` |
| `FormIFunction_Vperp_viscosity` | 1226-1242 | 17 | extern | LIVE | `src/mhd/mhd.c:316`, `tests/regression/t1_residuals.c:51`, `tests/regression/t1_residuals.c:245` |
| `initialize_ep` | 1273-1695 | 423 | static | LIVE | `src/mhd/ts_functions.c:1697`, `src/mhd/ts_functions.c:1702` |
| `FormIFunction_InitializeEP` | 1696-1698 | 3 | extern | LIVE | `src/mhd/ts_functions.c:3030`, `tests/regression/t1_residuals.c:53` |
| `FormIFunction_InitializeEP_halo` | 1701-1703 | 3 | extern | LIVE | `src/mhd/ts_functions.c:5409`, `tests/regression/t1_residuals.c:54` |
| `FormIFunction_newequilibrium_Vperp` | 1714-1719 | 6 | extern | LIVE | `src/mhd/ts_functions.c:5118`, `tests/regression/t1_residuals.c:52` |
| `FormRHSFunction_BImplicit` | 1723-1989 | 267 | extern | LIVE | `src/mhd/mhd.c:324`, `tests/regression/t1_residuals.c:265` |
| `FormInitialSolution` | 1991-3302 | 1312 | extern | LIVE | `src/mhd/mhd.c:679` |
| `FormExactSolution` | 3304-4323 | 1020 | extern | LIVE | `src/mhd/ts_functions.c:347`, `src/mhd/ts_functions.c:1349`, `src/mhd/ts_functions.c:4376` |
| `Monitor` | 4327-4532 | 206 | extern | LIVE | `src/mhd/ts_functions.c:5103`, `src/mhd/mhd.c:244`, `src/mhd/mhd.c:733` |
| `FormDummyIJacobian4` | 4534-4700 | 167 | extern | LIVE | `src/mhd/mhd.c:359` |
| `SampleShellPCSetUp` | 4702-4771 | 70 | extern | LIVE | `src/mhd/mhd.c:430` |
| `SampleShellPCApply` | 4775-4822 | 48 | extern | LIVE | `src/mhd/mhd.c:433` |
| `SampleShellPCDestroy` | 4824-4838 | 15 | extern | LIVE | `src/mhd/mhd.c:436` |
| `ReadInitialData` | 4852-4888 | 37 | extern | LIVE | `src/mhd/mhd.c:126`, `src/mhd/mhd.c:140`, `src/mhd/mhd.c:149` |
| `FormInitialSolution_psi` | 4896-5674 | 779 | extern | LIVE | `src/mhd/mhd.c:677` |
| `stag_vec_io` | 5697-5797 | 101 | extern | LIVE | `src/mhd/ts_functions.c:5084`, `src/mhd/ts_functions.c:5373`, `src/mhd/mhd.c:867` +3 more |

**Summary:** 22 functions, 0 ENTRY, 22 LIVE, 0 DEAD, 6 static.


## src/mhd/geometry.c

4301 lines, 17 functions.

| Function | Lines | Length | Linkage | Status | Called from |
|---|---|---|---|---|---|
| `cyldistance` | 29-33 | 5 | extern | LIVE | `src/mhd/ts_functions.c:181`, `src/mhd/ts_functions.c:189`, `src/mhd/ts_functions.c:191` +61 more |
| `surface` | 35-132 | 98 | extern | LIVE | `src/mhd/ts_functions.c:585`, `src/mhd/ts_functions.c:648`, `src/mhd/ts_functions.c:660` +130 more |
| `SaveSolution` | 142-333 | 192 | extern | LIVE | `src/mhd/mhd.c:703` |
| `SaveCoordinates` | 335-756 | 422 | extern | LIVE | `src/mhd/mhd.c:237` |
| `CellToVertexProjectionScalar` | 760-853 | 94 | extern | LIVE | `src/mhd/ts_functions.c:396`, `tests/regression/t0_operators.c:47` |
| `VertexToEdgeReconstruction` | 861-1097 | 237 | extern | LIVE | `src/mhd/ts_functions.c:472`, `tests/regression/t0_operators.c:48` |
| `VertexToFaceReconstruction` | 1101-1301 | 201 | extern | LIVE | `src/mhd/ts_functions.c:462`, `tests/regression/t0_operators.c:49` |
| `EdgeToCellReconstruction_r` | 1305-1376 | 72 | extern | LIVE | `src/mhd/geometry.c:3814` |
| `EdgeToCellReconstruction_phi` | 1378-1449 | 72 | extern | LIVE | `src/mhd/geometry.c:3860` |
| `EdgeToCellReconstruction_z` | 1451-1522 | 72 | extern | LIVE | `src/mhd/geometry.c:3890` |
| `FaceToVertexProjection` | 1528-2432 | 905 | extern | LIVE | `src/mhd/ts_functions.c:412`, `tests/regression/t0_operators.c:50` |
| `EdgeToVertexProjection` | 2438-3069 | 632 | extern | LIVE | `src/mhd/ts_functions.c:424`, `src/mhd/ts_functions.c:435`, `src/mhd/ts_functions.c:443` +4 more |
| `CellToFaceProjection` | 3075-3304 | 230 | extern | LIVE | `src/mhd/ts_functions.c:405`, `tests/regression/t0_operators.c:52` |
| `VertexCrossProduct` | 3306-3781 | 476 | extern | LIVE | `src/mhd/ts_functions.c:470` |
| `getEJArray` | 3783-4022 | 240 | extern | LIVE | `src/mhd/mhd.c:790`, `src/mhd/mhd.c:806` |
| `getVArray` | 4034-4155 | 122 | extern | LIVE | `src/mhd/mhd.c:812` |
| `getBArray` | 4160-4287 | 128 | extern | LIVE | `src/mhd/mhd.c:788`, `src/mhd/mhd.c:814` |

**Summary:** 17 functions, 0 ENTRY, 17 LIVE, 0 DEAD, 0 static.


## src/mhd/mimetic_operators.c

2339 lines, 13 functions.

| Function | Lines | Length | Linkage | Status | Called from |
|---|---|---|---|---|---|
| `FormDiscreteDivergence` | 34-175 | 142 | extern | LIVE | `src/mhd/ts_functions.c:4458`, `src/mhd/ts_functions.c:4497` |
| `ApplyDerivedDivergence` | 186-376 | 191 | extern | LIVE | `src/mhd/ts_functions.c:1191`, `src/mhd/ts_functions.c:1665` |
| `ApplyVectorLaplacian` | 382-638 | 257 | extern | LIVE | `src/mhd/ts_functions.c:385`, `tests/regression/t0_operators.c:44` |
| `FormPrimaryCurl` | 644-843 | 200 | extern | LIVE | `src/mhd/ts_functions.c:4438`, `src/mhd/ts_functions.c:5066`, `tests/regression/t0_operators.c:40` |
| `derived_curl` | 855-1129 | 275 | static | LIVE | `src/mhd/mimetic_operators.c:1132`, `src/mhd/mimetic_operators.c:1137`, `src/mhd/mimetic_operators.c:1142` |
| `FormDerivedCurl` | 1131-1133 | 3 | extern | LIVE | `src/mhd/ts_functions.c:4313`, `tests/regression/t0_operators.c:41` |
| `FormDerivedCurlnores` | 1136-1138 | 3 | extern | LIVE | `src/mhd/monitor_functions.c:1667`, `tests/regression/t0_operators.c:42` |
| `FormDerivedCurlnomp` | 1141-1143 | 3 | extern | LIVE | `src/mhd/ts_functions.c:423`, `src/mhd/geometry.c:3811`, `src/mhd/monitor_functions.c:289` +1 more |
| `FormSourceTermPotential` | 1148-1469 | 322 | extern | LIVE | `src/mhd/ts_functions.c:339`, `src/mhd/ts_functions.c:1341` |
| `FormDiscreteGradientEP` | 1471-1882 | 412 | extern | DEAD | - |
| `FormDiscreteGradientEP_noMat` | 1884-2061 | 178 | extern | LIVE | `src/mhd/ts_functions.c:357`, `src/mhd/ts_functions.c:1359`, `src/mhd/mimetic_operators.c:2332` +1 more |
| `FormDiscreteGradientVectorField` | 2065-2318 | 254 | extern | LIVE | `src/mhd/ts_functions.c:379`, `src/mhd/mimetic_operators.c:431` |
| `FormElectricField` | 2320-2339 | 20 | extern | LIVE | `src/mhd/geometry.c:3809`, `tests/regression/t0_operators.c:46` |

**Summary:** 13 functions, 0 ENTRY, 12 LIVE, 1 DEAD, 1 static.


## src/mhd/mass_matrix_coefficients.c

1271 lines, 24 functions.

| Function | Lines | Length | Linkage | Status | Called from |
|---|---|---|---|---|---|
| `betaf` | 26-78 | 53 | extern | LIVE | `src/mhd/ts_functions.c:856`, `src/mhd/ts_functions.c:857`, `src/mhd/ts_functions.c:858` +92 more |
| `betae_sum` | 93-291 | 199 | static | LIVE | `src/mhd/mass_matrix_coefficients.c:294`, `src/mhd/mass_matrix_coefficients.c:299`, `src/mhd/mass_matrix_coefficients.c:304` +3 more |
| `betae` | 293-295 | 3 | extern | LIVE | `src/mhd/mimetic_operators.c:1132`, `tests/regression/t0_coefficients.c:67` |
| `betaenores` | 298-300 | 3 | extern | LIVE | `src/mhd/mimetic_operators.c:1137`, `tests/regression/t0_coefficients.c:68` |
| `betaenomp` | 303-305 | 3 | extern | LIVE | `src/mhd/mimetic_operators.c:353`, `src/mhd/mimetic_operators.c:594`, `src/mhd/mimetic_operators.c:595` +3 more |
| `betae2` | 312-314 | 3 | extern | LIVE | `src/mhd/ts_functions.c:859`, `src/mhd/ts_functions.c:870`, `src/mhd/ts_functions.c:894` +19 more |
| `betaephi_isolcell` | 317-319 | 3 | extern | LIVE | `src/mhd/ts_functions.c:1702`, `tests/regression/t0_coefficients.c:71` |
| `betaeperp2` | 322-324 | 3 | extern | LIVE | `src/mhd/ts_functions.c:1702`, `tests/regression/t0_coefficients.c:72` |
| `betavnomp` | 331-552 | 222 | extern | LIVE | `src/mhd/mimetic_operators.c:353`, `src/mhd/mimetic_operators.c:594`, `src/mhd/mimetic_operators.c:595` +2 more |
| `alphaec_sum` | 569-650 | 82 | static | LIVE | `src/mhd/mass_matrix_coefficients.c:654`, `src/mhd/mass_matrix_coefficients.c:727`, `src/mhd/mass_matrix_coefficients.c:863` +2 more |
| `alphaec2` | 652-655 | 4 | extern | LIVE | `src/mhd/mass_matrix_coefficients.c:313`, `tests/regression/t0_coefficients.c:83` |
| `alphavcnomp` | 660-723 | 64 | extern | LIVE | `src/mhd/mass_matrix_coefficients.c:347`, `src/mhd/mass_matrix_coefficients.c:352`, `src/mhd/mass_matrix_coefficients.c:354` +87 more |
| `alphaec` | 725-728 | 4 | extern | LIVE | `src/mhd/mass_matrix_coefficients.c:294`, `tests/regression/t0_coefficients.c:82` |
| `alphaecnores` | 731-793 | 63 | extern | LIVE | `src/mhd/mass_matrix_coefficients.c:299`, `tests/regression/t0_coefficients.c:84` |
| `alphaecnomp` | 795-857 | 63 | extern | LIVE | `src/mhd/mass_matrix_coefficients.c:304`, `tests/regression/t0_coefficients.c:85` |
| `alphaecphi` | 861-864 | 4 | extern | LIVE | `tests/regression/t0_coefficients.c:87` |
| `alphaecperp2` | 867-870 | 4 | extern | LIVE | `src/mhd/mass_matrix_coefficients.c:323`, `tests/regression/t0_coefficients.c:86` |
| `alphaecphi_isolcell` | 875-882 | 8 | extern | LIVE | `src/mhd/mass_matrix_coefficients.c:318`, `tests/regression/t0_coefficients.c:88` |
| `alphafc` | 885-948 | 64 | extern | LIVE | `src/mhd/mass_matrix_coefficients.c:38`, `src/mhd/mass_matrix_coefficients.c:40`, `src/mhd/mass_matrix_coefficients.c:42` +9 more |
| `edge_average` | 958-1156 | 199 | static | LIVE | `src/mhd/mass_matrix_coefficients.c:1159`, `src/mhd/mass_matrix_coefficients.c:1268` |
| `rese` | 1158-1160 | 3 | extern | LIVE | `tests/regression/t0_coefficients.c:74` |
| `resec` | 1163-1213 | 51 | extern | LIVE | `src/mhd/mass_matrix_coefficients.c:1159`, `tests/regression/t0_coefficients.c:91` |
| `conduc` | 1215-1265 | 51 | extern | LIVE | `src/mhd/mass_matrix_coefficients.c:1268`, `tests/regression/t0_coefficients.c:92` |
| `condu` | 1267-1269 | 3 | extern | LIVE | `src/mhd/ts_functions.c:687`, `src/mhd/ts_functions.c:688`, `src/mhd/ts_functions.c:691` +73 more |

**Summary:** 24 functions, 0 ENTRY, 24 LIVE, 0 DEAD, 3 static.


## src/mhd/monitor_functions.c

1894 lines, 8 functions.

| Function | Lines | Length | Linkage | Status | Called from |
|---|---|---|---|---|---|
| `DumpVelocity_Cell` | 27-240 | 214 | extern | DEAD | - |
| `DumpSolution_Cell` | 242-918 | 677 | extern | LIVE | `src/mhd/ts_functions.c:4398` |
| `DumpError` | 922-1482 | 561 | extern | LIVE | `src/mhd/ts_functions.c:4422` |
| `DumpDivergence` | 1484-1555 | 72 | extern | LIVE | `src/mhd/ts_functions.c:4484` |
| `DumpLevelSet` | 1557-1618 | 62 | extern | LIVE | `src/mhd/ts_functions.c:4400` |
| `SaveIntermediateSolution` | 1620-1634 | 15 | extern | LIVE | `src/mhd/ts_functions.c:4384`, `src/mhd/ts_functions.c:5376`, `src/mhd/ts_functions.c:5383` |
| `ComputeCurrent` | 1636-1737 | 102 | extern | LIVE | `src/mhd/ts_functions.c:4508` |
| `DumpEdgeField` | 1741-1884 | 144 | extern | DEAD | - |

**Summary:** 8 functions, 0 ENTRY, 6 LIVE, 2 DEAD, 0 static.


## src/mhd/mhd.c

989 lines, 13 functions.

| Function | Lines | Length | Linkage | Status | Called from |
|---|---|---|---|---|---|
| `mhd_fail` | 64-72 | 9 | extern | ENTRY | (public API) |
| `mhd_PetscInit` | 74-87 | 14 | extern | ENTRY | `src/kinetic/avalanche.cpp:24`, `src/kinetic/hybrid.cpp:26`, `src/kinetic/profile.cpp:45` |
| `mhd_initialize` | 89-682 | 594 | extern | ENTRY | `src/kinetic/kinetic.cpp:249` |
| `mhd_step` | 684-740 | 57 | extern | ENTRY | `src/kinetic/kinetic.cpp:252`, `src/tasks/BackgroundFields.cpp:7` |
| `mhd_resetState` | 742-757 | 16 | extern | ENTRY | `src/tasks/BackgroundFields.cpp:13` |
| `mhd_getF` | 759-824 | 66 | extern | ENTRY | `src/tasks/Interpolate.cpp:17`, `src/tasks/Interpolate.cpp:18`, `src/tasks/Interpolate.cpp:19` |
| `mhd_destroy` | 826-849 | 24 | extern | ENTRY | `src/kinetic/hybrid.cpp:64`, `src/kinetic/hybrid.cpp:75` |
| `mhd_savesolution` | 852-872 | 21 | extern | ENTRY | `src/kinetic/kinetic.cpp:1196` |
| `mhd_loadsolution` | 874-899 | 26 | extern | ENTRY | `src/kinetic/kinetic.cpp:1208` |
| `mhd_save_hdf5` | 916-941 | 26 | extern | ENTRY | (public API) |
| `mhd_load_hdf5` | 943-968 | 26 | extern | ENTRY | (public API) |
| `view4d_zero` | 970-973 | 4 | extern | LIVE | `src/mhd/mhd.c:777` |
| `subview_exclude_d1` | 977-983 | 7 | extern | LIVE | `src/mhd/mhd.c:821` |

**Summary:** 13 functions, 11 ENTRY, 2 LIVE, 0 DEAD, 0 static.


## src/mhd/mfd_config.c

143 lines, 1 functions.

| Function | Lines | Length | Linkage | Status | Called from |
|---|---|---|---|---|---|
| `AppCtxView` | 10-142 | 133 | extern | LIVE | `src/mhd/mhd.c:112` |

**Summary:** 1 functions, 0 ENTRY, 1 LIVE, 0 DEAD, 0 static.


## src/mhd/default_petsc_options.c

150 lines, 1 functions.

| Function | Lines | Length | Linkage | Status | Called from |
|---|---|---|---|---|---|
| `default_petsc_options` | 107-150 | 44 | extern | LIVE | `src/mhd/mhd.c:110` |

**Summary:** 1 functions, 0 ENTRY, 1 LIVE, 0 DEAD, 0 static.


## Summary

| File | Lines | Functions | ENTRY | LIVE | DEAD | Static |
|---|---|---|---|---|---|---|
| `ts_functions.c` | 5797 | 22 | 0 | 22 | 0 | 6 |
| `geometry.c` | 4301 | 17 | 0 | 17 | 0 | 0 |
| `mimetic_operators.c` | 2339 | 13 | 0 | 12 | 1 | 1 |
| `mass_matrix_coefficients.c` | 1271 | 24 | 0 | 24 | 0 | 3 |
| `monitor_functions.c` | 1894 | 8 | 0 | 6 | 2 | 0 |
| `mhd.c` | 989 | 13 | 11 | 2 | 0 | 0 |
| `mfd_config.c` | 143 | 1 | 0 | 1 | 0 | 0 |
| `default_petsc_options.c` | 150 | 1 | 0 | 1 | 0 | 0 |
| **Total** | 16884 | 99 | 11 | 85 | 3 | 10 |

## Quarantine candidates (DEAD functions)

The following functions are unreachable from any root (public API in `mhd.h`, or any reference from `src/kinetic`, `src/tasks`, `tests/regression`) and are quarantine candidates:


### src/mhd/mimetic_operators.c

- `FormDiscreteGradientEP` (1471-1882, 412 lines, extern)

### src/mhd/monitor_functions.c

- `DumpVelocity_Cell` (27-240, 214 lines, extern)
- `DumpEdgeField` (1741-1884, 144 lines, extern)

## Phantom declarations

Names declared in `src/mhd/*.h` with no definition in any `build/src/CMakeFiles/mhd_core.dir/mhd/*.o` (checked via `nm --defined-only`):

None.

