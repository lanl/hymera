# MHD Solver Function Index

This document is a mechanically generated function index for the C MHD solver (`/workspace/src/mhd`). It was produced by scripted `grep`/`nm` extraction over the 9 files in scope (`ts_functions.c`, `geometry.c`, `mimetic_operators.c`, `mass_matrix_coefficients.c`, `monitor_functions.c`, `mhd.c`, `mfd_config.c`, `default_petsc_options.c`, and their matching `.h` headers), not by reading or judging the code. It reflects git commit `833eff2`.

**This document must be regenerated after any change to the source in `src/mhd`.** Line ranges, lengths, and LIVE/DEAD status will drift out of date the moment a function is added, removed, reordered, or gains/loses a caller.

Method: top-level function definitions were located with `grep -nE '^[A-Za-z_][A-Za-z0-9_ \*]*\([^;]*$'` in each `.c` file; a function's end line is the line before the next function's start line (or the file's last line, for the final function). Extracted counts matched the expected counts for every file (`ts_functions.c`=41, `geometry.c`=38, `mimetic_operators.c`=25, `mass_matrix_coefficients.c`=31, `monitor_functions.c`=14; `mhd.c`=14, `mfd_config.c`=1, `default_petsc_options.c`=1), and were cross-checked against `nm --defined-only` on the corresponding `build/src/CMakeFiles/mhd_core.dir/mhd/*.o` objects, which produced identical symbol counts and names for every file. LIVE/DEAD status for each function was determined by `grep -rn '\b<name>\b' src tests --include=*.c --include=*.h --include=*.cpp --include=*.hpp`, excluding the function's own definition line, declarations in `.h` files, and `PetscLogEventRegister("<name>", ...)` self-naming string literals inside the function's own body (these are log labels, not calls). A function with zero remaining ("genuine") references is DEAD.

**Correction applied after generation.** The reference filter above does not exclude commented-out code, so two functions were initially recorded LIVE on the strength of a `//`-prefixed call: `FormIFunction_DampingV` (`src/mhd/mhd.c:298`) and `Dump1stVertexField` (`src/mhd/ts_functions.c:17686`). Both are DEAD. Any regeneration of this document must skip lines whose first non-whitespace characters are `//` or `/*`.

## src/mhd/ts_functions.c

22656 lines, 41 functions.

| Function | Lines | Length | Status | Referenced from |
|---|---|---|---|---|
| `FormIJacobian_BImplicit` | 21-25 | 5 | LIVE | `src/mhd/mhd.c:280` |
| `FormIFunction_Inertia_V` | 26-973 | 948 | DEAD | - |
| `FormIFunction_Inertia_V_ni` | 974-1934 | 961 | DEAD | - |
| `FormIFunction_Inertia_viscosity` | 1935-2909 | 975 | DEAD | - |
| `FormIFunction_Vperp_viscosity` | 2910-3885 | 976 | LIVE | `src/mhd/mhd.c:299`, `src/mhd/ts_functions.c:17658`, `src/mhd/ts_functions.c:20968` +4 more |
| `FormIFunction_Vperp_viscosity_halo` | 3886-4866 | 981 | DEAD | - |
| `FormIFunction_Vperp_viscosity_halo_isolcell` | 4867-5847 | 981 | DEAD | - |
| `FormIFunction_Inertia` | 5848-6822 | 975 | DEAD | - |
| `FormIFunction2` | 6823-7516 | 694 | DEAD | - |
| `FormIFunction` | 7517-8192 | 676 | DEAD | - |
| `FormIFunction_BImplicit2` | 8193-9003 | 811 | DEAD | - |
| `FormIFunction_Initializepsi` | 9004-9399 | 396 | LIVE | `src/mhd/ts_functions.c:21944` |
| `FormIFunction_InitializeEP` | 9400-9955 | 556 | LIVE | `src/mhd/ts_functions.c:15212`, `src/mhd/ts_functions.c:22364` |
| `FormIFunction_InitializeEP_halo` | 9956-10515 | 560 | LIVE | `src/mhd/ts_functions.c:20718`, `src/mhd/ts_functions.c:21531` |
| `FormIFunction_InitializeEPV` | 10516-11303 | 788 | DEAD | - |
| `FormIFunction_newequilibrium` | 11304-12103 | 800 | DEAD | - |
| `FormIFunction_newequilibrium_Vperp` | 12104-12907 | 804 | LIVE | `src/mhd/ts_functions.c:20426`, `src/mhd/ts_functions.c:21236` |
| `FormIFunction_DampingV` | 12908-13805 | 898 | DEAD | - (only reference is the commented-out `src/mhd/mhd.c:298`) |
| `FormRHSFunction_BImplicit` | 13806-14153 | 348 | LIVE | `src/mhd/mhd.c:300`, `src/mhd/mhd.c:302`, `src/mhd/ts_functions.c:20427` +1 more |
| `FormInitialSolution` | 14154-15485 | 1332 | LIVE | `src/mhd/mhd.c:669`, `src/mhd/mhd.c:675`, `src/mhd/geometry.c:572` +1 more |
| `FormExactSolution` | 15486-16523 | 1038 | LIVE | `src/mhd/ts_functions.c:164`, `src/mhd/ts_functions.c:1112`, `src/mhd/ts_functions.c:2073` +16 more |
| `FormExactSolution_LargeData` | 16524-17551 | 1028 | LIVE | `src/mhd/ts_functions.c:19792` |
| `Monitor` | 17552-17823 | 272 | LIVE | `src/mhd/mhd.c:227`, `src/mhd/mhd.c:711`, `src/mhd/mhd.c:716` +9 more |
| `FormDummyIJacobian4` | 17824-17990 | 167 | LIVE | `src/mhd/mhd.c:337` |
| `SampleShellPCSetUp` | 17991-18087 | 97 | LIVE | `src/mhd/mhd.c:426` |
| `SampleShellPCSetUp_SuperLU` | 18088-18151 | 64 | DEAD | - |
| `SampleShellPCApply` | 18152-18199 | 48 | LIVE | `src/mhd/mhd.c:429` |
| `SampleShellPCDestroy` | 18200-18214 | 15 | LIVE | `src/mhd/mhd.c:432` |
| `SampleShellPCSetUp_Diag` | 18215-18273 | 59 | DEAD | - |
| `SampleShellPCApply_Diag` | 18274-18304 | 31 | DEAD | - |
| `SampleShellPCDestroy_Diag` | 18305-18317 | 13 | DEAD | - |
| `SampleShellPCSetUp_ApproximateDiag` | 18318-18403 | 86 | DEAD | - |
| `SampleShellPCApply_ApproximateDiag` | 18404-18445 | 42 | DEAD | - |
| `SampleShellPCDestroy_ApproximateDiag` | 18446-18460 | 15 | DEAD | - |
| `ReadInitialData` | 18461-18498 | 38 | LIVE | `src/mhd/mhd.c:109`, `src/mhd/mhd.c:123`, `src/mhd/mhd.c:132` +7 more |
| `ReadALine` | 18499-18548 | 50 | LIVE | `src/mhd/ts_functions.c:16960`, `src/mhd/ts_functions.c:16963`, `src/mhd/ts_functions.c:16968` +16 more |
| `FormInitialSolution_LargeData` | 18549-19664 | 1116 | DEAD | - |
| `FormIFunction_InitializeEP_LargeData` | 19665-20205 | 541 | LIVE | `src/mhd/ts_functions.c:19594` |
| `FormInitialSolution_psi` | 20206-20997 | 792 | LIVE | `src/mhd/mhd.c:117`, `src/mhd/mhd.c:665`, `src/mhd/mhd.c:673` |
| `FormInitialSolution_psi_fromNphi2` | 20998-21819 | 822 | DEAD | - |
| `FormInitialpsi` | 21820-22656 | 837 | DEAD | - |

## src/mhd/geometry.c

10111 lines, 38 functions.

| Function | Lines | Length | Status | Referenced from |
|---|---|---|---|---|
| `cyldistance` | 29-34 | 6 | LIVE | `src/mhd/ts_functions.c:319`, `src/mhd/ts_functions.c:327`, `src/mhd/ts_functions.c:329` +397 more |
| `surface` | 35-133 | 99 | LIVE | `src/mhd/ts_functions.c:422`, `src/mhd/ts_functions.c:477`, `src/mhd/ts_functions.c:489` +848 more |
| `EBoundaryAdjusters` | 134-555 | 422 | DEAD | - |
| `ComputeIsEBoundary` | 556-1025 | 470 | LIVE | `src/mhd/mimetic_operators.c:1359`, `src/mhd/mimetic_operators.c:1922` |
| `ComputeIsBBoundary` | 1026-1214 | 189 | DEAD | - |
| `ComputeIsCBoundary` | 1215-1306 | 92 | DEAD | - |
| `SaveSolution` | 1307-1498 | 192 | LIVE | `src/mhd/mhd.c:700` |
| `SaveCoordinates` | 1499-1920 | 422 | LIVE | `src/mhd/mhd.c:220` |
| `ReadDataInVec` | 1921-2264 | 344 | DEAD | - |
| `CellToVertexProjectionScalar` | 2265-2358 | 94 | LIVE | `src/mhd/ts_functions.c:1141`, `src/mhd/ts_functions.c:2112`, `src/mhd/ts_functions.c:3089` +7 more |
| `CellToVertexProjectionVector` | 2359-2543 | 185 | DEAD | - |
| `VertexToCellReconstruction` | 2544-2727 | 184 | DEAD | - |
| `VertexToEdgeReconstruction_scalar` | 2728-2884 | 157 | LIVE | `src/mhd/ts_functions.c:22224` |
| `VertexToEdgeReconstruction` | 2885-3121 | 237 | LIVE | `src/mhd/ts_functions.c:256`, `src/mhd/ts_functions.c:1211`, `src/mhd/ts_functions.c:2182` +11 more |
| `VertexToEdgeReconstructionMat` | 3122-3416 | 295 | DEAD | - |
| `VertexToFaceReconstruction` | 3417-3617 | 201 | LIVE | `src/mhd/ts_functions.c:246`, `src/mhd/ts_functions.c:1201`, `src/mhd/ts_functions.c:2172` +11 more |
| `VertexToFaceReconstructionMat` | 3618-3864 | 247 | DEAD | - |
| `EdgeToCellReconstruction_r` | 3865-3936 | 72 | LIVE | `src/mhd/geometry.c:8656` |
| `EdgeToCellReconstruction_phi` | 3937-4008 | 72 | LIVE | `src/mhd/geometry.c:8702` |
| `EdgeToCellReconstruction_z` | 4009-4080 | 72 | LIVE | `src/mhd/geometry.c:8732` |
| `EdgeToCellReconstructionMat` | 4081-4219 | 139 | DEAD | - |
| `FaceToCellReconstructionMat` | 4220-4316 | 97 | DEAD | - |
| `FaceToVertexProjection` | 4317-5221 | 905 | LIVE | `src/mhd/ts_functions.c:200`, `src/mhd/ts_functions.c:1155`, `src/mhd/ts_functions.c:2126` +12 more |
| `FaceToVertexProjectionMat` | 5222-6088 | 867 | DEAD | - |
| `EdgeToVertexProjection_Original` | 6089-6487 | 399 | DEAD | - |
| `EdgeToVertexProjection` | 6488-7119 | 632 | LIVE | `src/mhd/mimetic_operators.c:3441`, `src/mhd/mimetic_operators.c:3448`, `src/mhd/ts_functions.c:212` +44 more |
| `EdgeToVertexProjectionMat` | 7120-7714 | 595 | DEAD | - |
| `CellToFaceProjectionMat` | 7715-7919 | 205 | DEAD | - |
| `CellToFaceProjection` | 7920-8149 | 230 | LIVE | `src/mhd/ts_functions.c:193`, `src/mhd/ts_functions.c:1148`, `src/mhd/ts_functions.c:2119` +11 more |
| `VertexCrossProduct` | 8150-8625 | 476 | LIVE | `src/mhd/ts_functions.c:254`, `src/mhd/ts_functions.c:267`, `src/mhd/ts_functions.c:1209` +20 more |
| `getEJArray` | 8626-8875 | 250 | LIVE | `src/mhd/mhd.c:785`, `src/mhd/mhd.c:801` |
| `getVArray` | 8876-8997 | 122 | LIVE | `src/mhd/mhd.c:807` |
| `getJArray` | 8998-9119 | 122 | DEAD | - |
| `getBArray` | 9120-9248 | 129 | LIVE | `src/mhd/mhd.c:783`, `src/mhd/mhd.c:809` |
| `FromPetscVecToArray` | 9249-9679 | 431 | DEAD | - |
| `CellCoordArrays` | 9680-9810 | 131 | DEAD | - |
| `ScatterTest` | 9811-10053 | 243 | DEAD | - |
| `isInDomain` | 10054-10111 | 58 | DEAD | - |

## src/mhd/mimetic_operators.c

6985 lines, 25 functions.

| Function | Lines | Length | Status | Referenced from |
|---|---|---|---|---|
| `FormMaterialPropertiesMatrix` | 26-96 | 71 | LIVE | `src/mhd/mimetic_operators.c:1915` |
| `FormFaceMassMatrix` | 97-314 | 218 | DEAD | - |
| `FormEdgeMassMatrix` | 315-623 | 309 | DEAD | - |
| `FormDiscreteDivergence` | 624-765 | 142 | LIVE | `src/mhd/mimetic_operators.c:1331`, `src/mhd/mimetic_operators.c:1332`, `src/mhd/mimetic_operators.c:1333` +5 more |
| `FormGradDerivedDivergence` | 766-1301 | 536 | DEAD | - |
| `FormDerivedGradDivergence` | 1302-1842 | 541 | DEAD | - |
| `FormDerivedGradEtaDivergence` | 1843-2408 | 566 | DEAD | - |
| `FormDerivedDivergence` | 2409-2965 | 557 | LIVE | `src/mhd/mimetic_operators.c:789`, `src/mhd/mimetic_operators.c:790`, `src/mhd/mimetic_operators.c:791` |
| `ApplyDerivedDivergence` | 2966-3156 | 191 | LIVE | `src/mhd/mimetic_operators.c:3175`, `src/mhd/ts_functions.c:940`, `src/mhd/ts_functions.c:1901` +13 more |
| `ApplyDeltastar` | 3157-3184 | 28 | DEAD | - |
| `ApplyDeltastar2` | 3185-3292 | 108 | LIVE | `src/mhd/ts_functions.c:9318`, `src/mhd/ts_functions.c:22216` |
| `ApplyVectorLaplacian` | 3293-3549 | 257 | LIVE | `src/mhd/ts_functions.c:2103`, `src/mhd/ts_functions.c:3080`, `src/mhd/ts_functions.c:4054` +6 more |
| `FormDerivedGradient` | 3550-3855 | 306 | LIVE | `src/mhd/mimetic_operators.c:1344`, `src/mhd/mimetic_operators.c:1345`, `src/mhd/mimetic_operators.c:1346` +3 more |
| `FormDiscreteGradient` | 3856-4249 | 394 | LIVE | `src/mhd/mimetic_operators.c:782`, `src/mhd/mimetic_operators.c:783`, `src/mhd/mimetic_operators.c:784` |
| `FormPrimaryCurl` | 4250-4449 | 200 | LIVE | `src/mhd/ts_functions.c:17697`, `src/mhd/ts_functions.c:17728`, `src/mhd/ts_functions.c:20375` +2 more |
| `FormDerivedCurl` | 4450-4724 | 275 | LIVE | `src/mhd/ts_functions.c:15186`, `src/mhd/ts_functions.c:16494`, `src/mhd/ts_functions.c:17524` +1 more |
| `FormDerivedCurlnores` | 4725-4999 | 275 | LIVE | `src/mhd/monitor_functions.c:2505` |
| `FormDerivedCurlnomp` | 5000-5274 | 275 | LIVE | `src/mhd/ts_functions.c:211`, `src/mhd/ts_functions.c:1166`, `src/mhd/ts_functions.c:2137` +15 more |
| `FormDerivedCurlExt` | 5275-5608 | 334 | DEAD | - |
| `FormSourceTermPotential` | 5609-5930 | 322 | LIVE | `src/mhd/ts_functions.c:156`, `src/mhd/ts_functions.c:1104`, `src/mhd/ts_functions.c:2065` +15 more |
| `FormDiscreteGradientEP` | 5931-6342 | 412 | LIVE | `src/mhd/mimetic_operators.c:2903`, `src/mhd/monitor_functions.c:2527` |
| `FormDiscreteGradientEP_noMat` | 6343-6520 | 178 | LIVE | `src/mhd/mimetic_operators.c:6978`, `src/mhd/ts_functions.c:174`, `src/mhd/ts_functions.c:1122` +13 more |
| `FormDiscreteGradientEP_tilde` | 6521-6712 | 192 | LIVE | `src/mhd/mimetic_operators.c:3173` |
| `FormDiscreteGradientVectorField` | 6713-6966 | 254 | LIVE | `src/mhd/mimetic_operators.c:3341`, `src/mhd/ts_functions.c:188`, `src/mhd/ts_functions.c:1136` +8 more |
| `FormElectricField` | 6967-6985 | 19 | LIVE | `src/mhd/geometry.c:8651` |

## src/mhd/mass_matrix_coefficients.c

4231 lines, 31 functions.

| Function | Lines | Length | Status | Referenced from |
|---|---|---|---|---|
| `betaf` | 26-79 | 54 | LIVE | `src/mhd/mass_matrix_coefficients.c:75`, `src/mhd/mass_matrix_coefficients.c:129`, `src/mhd/ts_functions.c:705` +632 more |
| `betaf_wmp` | 80-133 | 54 | LIVE | `src/mhd/mimetic_operators.c:137`, `src/mhd/mimetic_operators.c:155`, `src/mhd/mimetic_operators.c:173` +9 more |
| `betae` | 134-333 | 200 | LIVE | `src/mhd/mass_matrix_coefficients.c:329`, `src/mhd/mimetic_operators.c:355`, `src/mhd/mimetic_operators.c:373` +194 more |
| `betaenores` | 334-529 | 196 | LIVE | `src/mhd/mass_matrix_coefficients.c:525`, `src/mhd/mimetic_operators.c:4885`, `src/mhd/mimetic_operators.c:4886` +29 more |
| `betaenomp` | 530-725 | 196 | LIVE | `src/mhd/mass_matrix_coefficients.c:721`, `src/mhd/mimetic_operators.c:2573`, `src/mhd/mimetic_operators.c:2587` +45 more |
| `betaeperp` | 726-925 | 200 | LIVE | `src/mhd/mass_matrix_coefficients.c:921` |
| `betaephi` | 926-1125 | 200 | LIVE | `src/mhd/mass_matrix_coefficients.c:1121` |
| `betae2` | 1126-1325 | 200 | LIVE | `src/mhd/mass_matrix_coefficients.c:1321`, `src/mhd/ts_functions.c:3601`, `src/mhd/ts_functions.c:3618` +18 more |
| `betaephi_isolcell` | 1326-1525 | 200 | LIVE | `src/mhd/mass_matrix_coefficients.c:1521`, `src/mhd/ts_functions.c:5582`, `src/mhd/ts_functions.c:10383` |
| `betaeperp2` | 1526-1725 | 200 | LIVE | `src/mhd/mass_matrix_coefficients.c:1721`, `src/mhd/ts_functions.c:4612`, `src/mhd/ts_functions.c:4622` +22 more |
| `betaephi2` | 1726-1925 | 200 | LIVE | `src/mhd/mass_matrix_coefficients.c:1921`, `src/mhd/ts_functions.c:4601` |
| `betav` | 1926-2148 | 223 | LIVE | `src/mhd/mass_matrix_coefficients.c:2144` |
| `betavnomp` | 2149-2371 | 223 | LIVE | `src/mhd/mass_matrix_coefficients.c:2367`, `src/mhd/mimetic_operators.c:2743`, `src/mhd/mimetic_operators.c:2763` +10 more |
| `alphaec2` | 2372-2486 | 115 | LIVE | `src/mhd/mass_matrix_coefficients.c:1143`, `src/mhd/mass_matrix_coefficients.c:1146`, `src/mhd/mass_matrix_coefficients.c:1149` +53 more |
| `alphavc` | 2487-2592 | 106 | LIVE | `src/mhd/mass_matrix_coefficients.c:1942`, `src/mhd/mass_matrix_coefficients.c:1947`, `src/mhd/mass_matrix_coefficients.c:1949` +86 more |
| `alphavcnomp` | 2593-2657 | 65 | LIVE | `src/mhd/mass_matrix_coefficients.c:2165`, `src/mhd/mass_matrix_coefficients.c:2170`, `src/mhd/mass_matrix_coefficients.c:2172` +86 more |
| `alphaec` | 2658-2772 | 115 | LIVE | `src/mhd/mass_matrix_coefficients.c:151`, `src/mhd/mass_matrix_coefficients.c:154`, `src/mhd/mass_matrix_coefficients.c:157` +53 more |
| `alphaecnores` | 2773-2836 | 64 | LIVE | `src/mhd/mass_matrix_coefficients.c:347`, `src/mhd/mass_matrix_coefficients.c:350`, `src/mhd/mass_matrix_coefficients.c:353` +52 more |
| `alphaecnomp` | 2837-2900 | 64 | LIVE | `src/mhd/mass_matrix_coefficients.c:543`, `src/mhd/mass_matrix_coefficients.c:546`, `src/mhd/mass_matrix_coefficients.c:549` +52 more |
| `alphaecperp` | 2901-3015 | 115 | LIVE | `src/mhd/mass_matrix_coefficients.c:743`, `src/mhd/mass_matrix_coefficients.c:746`, `src/mhd/mass_matrix_coefficients.c:749` +53 more |
| `alphaecphi` | 3016-3130 | 115 | LIVE | `src/mhd/mass_matrix_coefficients.c:943`, `src/mhd/mass_matrix_coefficients.c:946`, `src/mhd/mass_matrix_coefficients.c:949` +54 more |
| `alphaecperp2` | 3131-3245 | 115 | LIVE | `src/mhd/mass_matrix_coefficients.c:1543`, `src/mhd/mass_matrix_coefficients.c:1546`, `src/mhd/mass_matrix_coefficients.c:1549` +53 more |
| `alphaecphi2` | 3246-3360 | 115 | LIVE | `src/mhd/mass_matrix_coefficients.c:1743`, `src/mhd/mass_matrix_coefficients.c:1746`, `src/mhd/mass_matrix_coefficients.c:1749` +53 more |
| `alphaecphi_isolcell` | 3361-3493 | 133 | LIVE | `src/mhd/mass_matrix_coefficients.c:1343`, `src/mhd/mass_matrix_coefficients.c:1346`, `src/mhd/mass_matrix_coefficients.c:1349` +52 more |
| `alphafc` | 3494-3558 | 65 | LIVE | `src/mhd/mass_matrix_coefficients.c:38`, `src/mhd/mass_matrix_coefficients.c:40`, `src/mhd/mass_matrix_coefficients.c:42` +8 more |
| `alphafc_wmp` | 3559-3663 | 105 | LIVE | `src/mhd/mass_matrix_coefficients.c:92`, `src/mhd/mass_matrix_coefficients.c:94`, `src/mhd/mass_matrix_coefficients.c:96` +8 more |
| `alphac` | 3664-3727 | 64 | LIVE | `src/mhd/mimetic_operators.c:3709` |
| `rese` | 3728-3927 | 200 | LIVE | `src/mhd/mass_matrix_coefficients.c:3923`, `src/mhd/monitor_functions.c:2539`, `src/mhd/monitor_functions.c:2540` +10 more |
| `resec` | 3928-3979 | 52 | LIVE | `src/mhd/mass_matrix_coefficients.c:3745`, `src/mhd/mass_matrix_coefficients.c:3748`, `src/mhd/mass_matrix_coefficients.c:3751` +52 more |
| `conduc` | 3980-4031 | 52 | LIVE | `src/mhd/mass_matrix_coefficients.c:4049`, `src/mhd/mass_matrix_coefficients.c:4052`, `src/mhd/mass_matrix_coefficients.c:4055` +52 more |
| `condu` | 4032-4231 | 200 | LIVE | `src/mhd/mass_matrix_coefficients.c:4227`, `src/mhd/mimetic_operators.c:5884`, `src/mhd/mimetic_operators.c:5885` +180 more |

## src/mhd/monitor_functions.c

3352 lines, 14 functions.

| Function | Lines | Length | Status | Referenced from |
|---|---|---|---|---|
| `DumpVelocity_Cell` | 29-242 | 214 | LIVE | `src/mhd/ts_functions.c:17637`, `src/mhd/ts_functions.c:17643` |
| `DumpSolution_Cell` | 243-919 | 677 | LIVE | `src/mhd/ts_functions.c:17622`, `src/mhd/ts_functions.c:17662`, `src/mhd/ts_functions.c:17698` +15 more |
| `DumpSolution` | 920-1764 | 845 | LIVE | `src/mhd/ts_functions.c:15469`, `src/mhd/ts_functions.c:19651`, `src/mhd/ts_functions.c:20691` +5 more |
| `DumpError` | 1765-2325 | 561 | LIVE | `src/mhd/mimetic_operators.c:270`, `src/mhd/ts_functions.c:17712` |
| `DumpDivergence` | 2326-2397 | 72 | LIVE | `src/mhd/ts_functions.c:17774` |
| `DumpLevelSet` | 2398-2459 | 62 | LIVE | `src/mhd/ts_functions.c:17624` |
| `SaveIntermediateSolution` | 2460-2474 | 15 | LIVE | `src/mhd/ts_functions.c:17608`, `src/mhd/ts_functions.c:20685`, `src/mhd/ts_functions.c:20692` |
| `ComputeCurrent` | 2475-2643 | 169 | LIVE | `src/mhd/ts_functions.c:17798` |
| `Dump1stVertexField` | 2644-2998 | 355 | DEAD | - (only reference is the commented-out `src/mhd/ts_functions.c:17686`) |
| `DumpEdgeField` | 2999-3142 | 144 | LIVE | `src/mhd/ts_functions.c:17678`, `src/mhd/geometry.c:8654`, `src/mhd/geometry.c:8655` |
| `multiplybyR` | 3143-3149 | 7 | DEAD | - |
| `createHermiteFD` | 3150-3183 | 34 | DEAD | - |
| `getHermiteDataFD` | 3184-3267 | 84 | LIVE | `src/mhd/monitor_functions.c:3173` |
| `DumpPsi_Cell` | 3268-3352 | 85 | DEAD | - |

## src/mhd/mhd.c

1041 lines, 14 functions.

| Function | Lines | Length | Status | Referenced from |
|---|---|---|---|---|
| `mhd_PetscInit` | 56-67 | 12 | ENTRY | `src/mhd/mhd.cpp:26`, `src/kinetic/hybrid.cpp:26`, `src/kinetic/profile.cpp:45` +1 more |
| `mhd_initialize` | 68-679 | 612 | ENTRY | `src/mhd/default_petsc_options.c:7`, `src/mhd/mhd.c:60`, `src/mhd/mhd.c:93` +7 more |
| `mhd_step` | 680-738 | 59 | ENTRY | `src/tasks/BackgroundFields.cpp:7`, `src/kinetic/kinetic.cpp:249` |
| `mhd_resetState` | 739-754 | 16 | ENTRY | `src/tasks/BackgroundFields.cpp:13` |
| `mhd_getF` | 755-820 | 66 | ENTRY | `src/tasks/Interpolate.cpp:17`, `src/tasks/Interpolate.cpp:18`, `src/tasks/Interpolate.cpp:19` |
| `mhd_destroy` | 821-863 | 43 | ENTRY | `src/mhd/mhd.cpp:67`, `src/kinetic/hybrid.cpp:64`, `src/kinetic/hybrid.cpp:75` |
| `stag_vec_io` | 864-929 | 66 | LIVE | `src/mhd/mhd.c:944`, `src/mhd/mhd.c:965`, `src/mhd/mhd.c:991` +1 more |
| `mhd_savesolution` | 930-950 | 21 | ENTRY | `src/kinetic/kinetic.cpp:1193` |
| `mhd_loadsolution` | 951-976 | 26 | ENTRY | `src/kinetic/kinetic.cpp:1205` |
| `mhd_save_hdf5` | 977-997 | 21 | DEAD | - |
| `mhd_load_hdf5` | 998-1018 | 21 | DEAD | - |
| `view4d_zero` | 1019-1023 | 5 | LIVE | `src/mhd/mhd.c:51`, `src/mhd/mhd.c:772` |
| `view3d_zero` | 1024-1028 | 5 | LIVE | `src/mhd/mhd.c:52` |
| `subview_exclude_d1` | 1029-1041 | 13 | LIVE | `src/mhd/mhd.c:54`, `src/mhd/mhd.c:816` |

## src/mhd/mfd_config.c

141 lines, 1 functions.

| Function | Lines | Length | Status | Referenced from |
|---|---|---|---|---|
| `AppCtxView` | 10-141 | 132 | LIVE | `src/mhd/mhd.c:91`, `src/mhd/mhd.c:93` |

## src/mhd/default_petsc_options.c

150 lines, 1 functions.

| Function | Lines | Length | Status | Referenced from |
|---|---|---|---|---|
| `default_petsc_options` | 107-150 | 44 | LIVE | `src/mhd/default_petsc_options.c:2`, `src/mhd/mhd.c:47`, `src/mhd/mhd.c:89` |

## Phantom declarations

These 11 names are declared in a header under `src/mhd/` but defined nowhere in the repository (no matching symbol in `nm --defined-only` on any `mhd_core` object file, and no function body found by grep across `src/` and `tests/`). Anything that calls them would fail to link if it were ever compiled in; nothing currently does.

| Declared name | Header:line |
|---|---|
| `Computepsi` | `src/mhd/geometry.h:61` |
| `FromPetscVecToArray_EfieldCell` | `src/mhd/geometry.h:62` |
| `interpolateHermite` | `src/mhd/monitor_functions.h:47` |
| `interpolate2D` | `src/mhd/monitor_functions.h:48` |
| `interpolatey` | `src/mhd/monitor_functions.h:49` |
| `interpolatex` | `src/mhd/monitor_functions.h:50` |
| `integrateHermite_1D_r` | `src/mhd/monitor_functions.h:51` |
| `integrateHermite_2D_z` | `src/mhd/monitor_functions.h:52` |
| `evalHermite1D` | `src/mhd/monitor_functions.h:53` |
| `evalHermite2D` | `src/mhd/monitor_functions.h:54` |
| `TSAdaptChoose_user` | `src/mhd/monitor_functions.h:56` |
| `Update_J_RE` | `src/mhd/ts_functions.h:52` |

Note: this is 12 names, one more than the ~11 anticipated; all 12 were independently confirmed absent from both `nm --defined-only` and a repo-wide grep for a function body.

## Summary

| File | Lines | Functions | LIVE | ENTRY | DEAD | Dead lines |
|---|---|---|---|---|---|---|
| `ts_functions.c` | 22656 | 41 | 19 | 0 | 22 | 13573 |
| `geometry.c` | 10111 | 38 | 19 | 0 | 19 | 5245 |
| `mimetic_operators.c` | 6985 | 25 | 18 | 0 | 7 | 2532 |
| `mass_matrix_coefficients.c` | 4231 | 31 | 31 | 0 | 0 | 0 |
| `monitor_functions.c` | 3352 | 14 | 10 | 0 | 4 | 481 |
| `mhd.c` | 1041 | 14 | 4 | 8 | 2 | 42 |
| `mfd_config.c` | 141 | 1 | 1 | 0 | 0 | 0 |
| `default_petsc_options.c` | 150 | 1 | 1 | 0 | 0 | 0 |
| **Total** | **48667** | **165** | **103** | **8** | **54** | **21873** |
