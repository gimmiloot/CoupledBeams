# NLSP: local 3D FEM infrastructure readiness

2026-10-08, local Windows checkout `D:\PHD\CoupledBeams\CoupledBeams`.
Initial state: main, HEAD `7d1363cc47d9cea426515fa0783a0643b47a5361`, clean.
**AUDIT_COMPLETE_WITH_QUALIFICATIONS.** Gmsh and an already unpacked CCX2.22
execute. Linear solid eigenanalysis was previously used. Geometric nonlinear
static/transient are documented, but new project workflows need adaptation and
`NEEDS_EXECUTION_TEST`. This audit ran no FEM job, mesh generation or installation.

## A. Programs and bounded discovery

D:\PHD has `CalculiX-Windows-master`, `CoupledBeams`, `gmsh-4.15.2-Windows64`.
Known paths first; then software depth6/project depth5,628 directories, no link
descent. The initial depth limit missed a deep historical cache, subsequently
checked at its exact known path. No unrestricted disk traversal occurred.

Full shared path used below:
`CCX_BIN=D:\PHD\CoupledBeams\CoupledBeams\results\_smoke\3d_fem_environment_check\calculix_2p22\CalculiX-2.22.0-win-x64\bin`.

| Program | Path | Version/source | Current evidence | Purpose |
|---|---|---|---|---|
| Gmsh | `D:\PHD\gmsh-4.15.2-Windows64\gmsh.exe` | 4.15.2 via --version | PE x64; rc0,0.109s; FOUND_FILE/EXECUTABLE_CONFIRMED/VERSION_CONFIRMED | Solid meshing/export |
| CCX | `CCX_BIN\ccx.exe` | 2.22 via -v, build/source | PE x64; rc201,0.180s; FOUND_FILE/EXECUTABLE_CONFIRMED/VERSION_CONFIRMED | Single-threaded solver |
| CCX MT | `CCX_BIN\ccx_MT.exe` | 2.22 via -v, build/source | PE x64; rc201,0.130s; FOUND_FILE/EXECUTABLE_CONFIRMED/VERSION_CONFIRMED | Multithreaded solver |
| CGX | `CCX_BIN\cgx.exe` | Package buildInfo2.22, not runtime-probed | PE x64; FOUND_FILE/CAPABILITY_NOT_TESTED; no GUI | Optional pre/postprocessor |
| Registered ANSYS Discovery/SpaceClaim | `D:\Ansys\ANSYS Inc\v202\SCDM\SpaceClaim.exe` | Registry20.2; file2020.2.40960.2818/product2020.2.0.0 | FOUND_FILE; launch/license/solver capabilities UNKNOWN | Detected CAD/frontend, not confirmed FEM solver |

CCX -v prints version then calls Fortran stop: local ccx_2.22.c:129-131 and
stop.f (read from source tar stream) establish expected exit201. This is not a
loader failure. Probes used timeout10s, CREATE_NO_WINDOW, empty temporary cwd;
CCX created zero files. No archive executable was extracted or launched.

Absent historical paths: `D:\PHD\gmsh\gmsh-4.15.2-Windows64\gmsh.exe`,
`D:\PHD\calculix`, `D:\PHD\calculix\calculix_2.22_4win\ccx_static.exe`.
No ccx_static/ccx_dynamic found in inspected files/catalogs. The actual CCX/MT
suffix distinction is threading, confirmed by README/Makefiles, not static vs
dynamic capability. Unavailable historical names do not prove a solver absent.

GMSH_EXE/CCX_EXE unset at Process/User/Machine levels; neither executable on PATH.
Absolute paths work; no persistent environment change is necessary. Code_Aster,
Salome-Meca and Elmer were not found in bounded D:\PHD/PATH/filtered app checks.
This is not an exhaustive search of every disk/environment. Windows filesystem
is accessible; LOCAL_WINDOWS_FILESYSTEM_UNAVAILABLE does not apply.

Release ZIP catalogs under `D:\PHD\CalculiX-Windows-master\releases` contain
2.20/2.21/2.22/2.23 and GE-OSS2.9/2.10, CCX/MT/CGX and runtime files. These are
distributions, not extra execution-confirmed installations; none unpacked/copied.
Existing CCX/MT binaries match corresponding2.22 ZIP streams exactly.

| Artifact | SHA256 |
|---|---|
| gmsh.exe | `317c43391e5b1fab3a1dd80dc5245dad6e2d087910f4b8ebc234bd6d4b8f41a1` |
| ccx.exe | `b29db0909b36c01abd876015c62e793efe73e158bfd417deb86ac48f771c9cb8` |
| ccx_MT.exe | `4d10c3025ebb876c333e27ae46c14e515867da8932aab259ee1676f4f7b38264` |
| ccx_2.22.pdf | `56963f827422ec7663cf218b60fffded19fd6ccebab793d2ccba667227d19d39` |

CCX_BIN/buildInfo.log:11Sep2024 19:48:39, CCX/CGX2.22, MinGW x86_64-win32-seh
GNU4.8.2. Both builds use SPOOLES/ARPACK; MT adds USE_MT/spoolesMT.
Imported non-system DLLs are present: pthreadGC2, libgcc_s_seh-1,
libwinpthread-1, libgfortran-3; ccx.exe also imports libgomp-1. Bin additionally
has libquadmath-0, libstdc++-6, glut64. KERNEL32/msvcrt are Windows system DLLs.
Successful probes establish loader readiness, not every analysis procedure.
The runtime resides in ignored old results: it is not guaranteed by a fresh clone.
No relocation/installation was done. Gmsh's probe needs no extra DLL setup here.

## Python

Historical interpreter exists: `D:\python\Pycharm\pythonProject\.venv\Scripts\python.exe`.
Python3.12.4 AMD64, NumPy2.1.3/SciPy1.15.2 import successfully.

| Package in THIS interpreter | Status |
|---|---|
| gmsh | NOT_INSTALLED |
| meshio | NOT_INSTALLED |
| pyvista | NOT_INSTALLED |
| sfepy | NOT_INSTALLED |
| dolfinx | NOT_INSTALLED |
| fenics | NOT_INSTALLED |
| petsc4py | NOT_INSTALLED |
| slepc4py | NOT_INSTALLED |

Module/distribution discovery does not establish a nonlinear solver by import.
Other interpreters were not inventoried. These missing optional packages do not
block the existing CLI workflow/custom INP/DAT/FRD readers. Nothing installed.
Existing `D:\texstudio\miktex\miktex\bin\x64\pdftotext.exe` reads local manuals.

## Historical calculations and angular difficulties

Tracked sources: [FEM status](../thickness_mismatch/fem_validation_status.md),
[May26 audit](../thickness_mismatch/solid_fem_audit.md),
[orthotropic anchor](../anisotropic_rods/yartsev_ch2_limited_3d_fem_anchor_design.md).
Their historical statuses are not rewritten. The old absent-tools snapshot is
not present availability; missing E3/nu13/nu23 in the anisotropic design does not
block the specified isotropic G20.

Actual retained inputs/results/logs read now:

- `_smoke/3d_fem_environment_check`: old Gmsh4.15.2/CCX2.22 jobs rc0,
  circular L2/r.04/E=rho1/nu.3, C3D10,2404 nodes/1108 solids,40 DAT frequencies
  and40 FRD1PMODE/DISP blocks. One negative-Jacobian warning: execution success
  is NOT mesh-quality certification.
- `eb_vs_timoshenko_3d_validation/uniform_beta0_eps0p05`: circular L2/r.1,
  30673 nodes/19239 solids,one component,60 modes,zero recorded warnings;
  four classified circular bending pairs closer to Timoshenko.
- Retained stepped-cylinder reports: mu.5/eta.1/epsilon.0025/.01/.05;
  eps.01 refinement has92150 nodes/56776 solids,60 modes,zero recorded warnings.

Fused/point-joint and some convergence history is documented in tracked reports,
but old raw paths are now absent and were NOT reverified. No old result changed.
Angular problems combine distinct causes: finite fused overlap/spherical fallback
changes the physical joint; rigid end faces suppress deformation differently from
an ideal1D point joint; planar U3 clamps could conflict with rigid-body dependent
nodes; close clusters/duplicate or weak MAC complicate mode identification;
Jacobian quality/doublet splitting are separate numerical issues. Straight beta0
checks do not suggest generic solver failure. These joint issues are not blockers
for one monolithic rectangular straight rod.

## B. Capability matrix

Primary source: [local CalculiX2.22 manual](D:/PHD/CalculiX-Windows-master/src/downloads/ccx_2.22.pdf),
Guido Dhondt,5Aug2024, verified hash above; printed=one-based PDF page numbers.
Shipped example input decks were read only. No web/new-version extrapolation.

| Capability | Status | Evidence / gap |
|---|---|---|
| Solid mesh/export | SUPPORTED_AND_PREVIOUSLY_USED | Quadratic tetrahedra, MSH4.1/INP, current Gmsh probe |
| Rectangular Box | SUPPORTED_BY_LOCAL_DOCUMENTATION_NOT_TESTED | Gmsh tutorials/t16.geo:13,19 OCC/Box; no G20 mesh |
| Quadratic solids/E/nu/rho/full-face clamps | SUPPORTED_AND_PREVIOUSLY_USED | C3D10 old inputs; manual §6.2.7 pp99-100 |
| Linear frequencies/eigenvectors/nodal export | SUPPORTED_AND_PREVIOUSLY_USED | Old DAT/FRD; *FREQUENCY §7.62 pp511-513; *NODE FILE pp557-560 |
| Multiple mesh levels | SUPPORTED_AND_PREVIOUSLY_USED | Old straight refinement; reusable convergence code, not new G20 convergence |
| Geometric nonlinear static/increments | SUPPORTED_BY_LOCAL_DOCUMENTATION_NOT_TESTED | *STATIC pp594-595, *STEP NLGEOM p600; project deck/parser missing |
| Direct nonlinear transient+NLGEOM | SUPPORTED_BY_LOCAL_DOCUMENTATION_NOT_TESTED | p600 explicitly allows DYNAMIC; pp477-479; local beamnldy.inp uses this combination |
| Initial nodal U/V | SUPPORTED_BY_LOCAL_DOCUMENTATION_NOT_TESTED | *INITIAL CONDITIONS pp534-537, global components |
| Preload-to-free-motion protocol | LIKELY_SUPPORTED_NEEDS_SMALL_TEST | Local STATIC→DYNAMIC examples; CLOAD OP=NEW p423; exact state/release NEEDS_EXECUTION_TEST |
| Full transient fields/timepoints | SUPPORTED_BY_LOCAL_DOCUMENTATION_NOT_TESTED | NODE FILE U,V/TIME POINTS pp558-560; increment-aware reader missing |
| Strains/stresses/reactions | SUPPORTED_BY_LOCAL_DOCUMENTATION_NOT_TESTED | EL FILE pp490-494, NODE FILE/PRINT pp558,562; nonlinear readers missing |
| Internal/kinetic energy export | SUPPORTED_BY_LOCAL_DOCUMENTATION_NOT_TESTED | EL PRINT ELSE/ELKE/TOTALS pp495-496; energy-control workflow missing |
| NLSP-specific theta/c recovery | UNKNOWN | Requires physical section-averaging definition, not a solver keyword |

Documented cautions for FUTURE tests:

- *ELASTIC+NLGEOM uses Lagrangian strain/PK2 stress (St-Venant–Kirchhoff,
  §6.8.1 p255). Declare the3D constitutive contract/small-strain regime;
  it is not automatically V0 or its full3D derivation. No material changed here.
- Implicit DYNAMIC defaults ALPHA=-.05 numerical damping (p477): free-motion
  energy comparison needs a declared integration contract and temporal check.
- ENER/ELSE requests start at the FIRST nonlinear step. ELSE internal and
  ELKE kinetic energy separately are not a full energy balance. RF includes
  reactions AND applied nodal/distributed loads; do not sum it blindly as reactions.
- Loads persist across nonperturbative steps unless removed. CLOAD OP=NEW
  removes prior concentrated loads; other loads need their own release.
  Initial displacement alone does not establish consistent initial stresses.

## C. Reusable code and D. Required work

| File/functions | Reuse | Gap |
|---|---|---|
| [single rod](../../scripts/analysis/solid_fem_single_rod_fixed_fixed.py): CasePaths, generate_mesh_with_gmsh_cli | Paths/CLI MSH+INP | Fresh output, explicit current exe paths |
| Same: write_gmsh_geo | OCC pattern | Cylinder only; small Box/config adapter |
| Same: read_gmsh_inp_mesh_data, coordinate_end_node_sets | Nodes/connectivity, x=0,L faces | Rectangular geometry contract |
| Same: write_calculix_mesh_include, write_calculix_template, write_calculix_inputs | C3D10 modal material/density/clamps | New rectangular case/config/normalization |
| Same: run_calculix | CLI/log/timeout | Modal-only; deletes stale OWN case outputs, never target old results |
| Same: parse_calculix_frd_mode_shapes | 1PMODE/DISP eigenvectors | Does not cover static/transient increments |
| [spectrum extraction](../../scripts/analysis/thickness_mismatch/audits/audit_full_spectrum_3d_fem_smoke_extraction.py): parse_calculix_frequency_table | Bounded eigen table, eigenvalue/RAD-TIME/CYCLES-TIME | Prefer over broad numeric fallback |
| [convergence](../../scripts/analysis/thickness_mismatch/audits/run_straight_uniform_3d_mesh_convergence.py): parse_mesh_quality, configure_single_module, run_one_mesh | Multi-h/report/log pattern | Circular defaults; quality not full independent Jacobian audit |
| [beta0 validator](../../scripts/analysis/thickness_mismatch/audits/validate_eb_timo_3d_beta0_uniform_eps0p05.py) | Old mesh/spectrum reports | Circular EB/Timo, not new rectangular NLSP |
| [point-joint](../../scripts/analysis/solid_fem_coupled_equal_rods_point_joint.py): solid_centerline_vector, centerline_mac, build_hungarian_assignment_rows | Shape-overlap ideas | Nodal-count means/nearest fills; new section recovery needed |
| [python_fem.py](../../src/my_project/fem/python_fem.py) | 1D axial/EB frame, DOFs u,v,theta | NOT solid/nonlinear3D solver |

INP reader accepts C3D* but drops per-element type; writer hardcodes C3D10.
Reuse quadratic tetrahedra, not silently C3D20 connectivity. Circular radius
summaries/doublet pairing must be replaced for a rectangle. Log warnings do not
provide complete determinant/volume/mass quality certification. Inspected FEM
writers contain no NLGEOM/STATIC/DYNAMIC/INITIAL CONDITIONS workflows.
The [environment CLI](../../scripts/analysis/thickness_mismatch/audits/check_3d_fem_environment.py)
is ACTIVE mesh/modal smoke, not passive inventory; it was NOT run. The
[fused-solid script](../../scripts/analysis/solid_fem_coupled_equal_rods.py) is
an old joint diagnostic, not the straight-rod route. No solver code was edited.

## Future comparison contract

Future monolithic rectangle: L=1, b=0.20, h=0.05, E=rho=1, nu=0.3; x along length, both full faces
fixed, no internal joint. No mesh/input was created. Solver units are user-defined
(§3 pp14-15); G20 numbers require a consistent normalized time/material contract.
DAT frequency is cycles/time, omega=2pi f; not experimental Hz by default.
Do not reuse circular Lambda conversion or renormalize each solid separately.

Linear observables: first bending/axial frequencies, nodal forms, mesh convergence.
Static: deflection, axial motion, section rotation, clamp reactions and nonlinear
correction. Transient: free histories/full fields, nonlinear bending corrections,
work/energy. No future eigenvector-amplitude fitting to disguise differences. Section-axis
and bending-plane mapping must follow the existing NLSP local-frame contract
before generating the rectangle; circular symmetry cannot hide an axis interchange.
M-H c is NOT a solid nodal DOF: physical section averaging for c/theta is still
needed and not implemented here. Full3D face clamp and1D c=0 are not identical
constraints on every local strain/warping; end regions can differ. V0 remains
an adopted reduced nonlinear law, not full3D truth.

## E. Minimal next test and effort

No program installation is presently necessary for Gmsh+CCX. Rough development
effort, NOT runtime guarantees: FEM-1 about0.5-1 focused working day for a Box/config
adapter, rectangular checks and reused modal/parsing workflow, then convergence;
FEM-2 an additional1-2 days for nonlinear load/increment/output readers and tests;
FEM-3 several additional days (roughly3-5 initially), potentially longer for
consistent3D initial-state recovery and verified spatial/temporal/energy comparisons.

ONE proposed next benchmark, NOT run: linear fixed-fixed G20 rectangular C3D10
solid with three prespecified mesh levels. New output namespace; inspect geometry,
connectivity/clamped faces/Jacobians; retain full raw frequencies/nodal forms;
identify bending planes/axial modes without circular doublet averaging and report
mesh convergence. Exact policies belong to a separately authorized task.
No FEM-2/FEM-3 execution is automatically authorized by readiness.

## Safe audit reproduction

Only inventory/version/doc reading; timeout10s, no GUI or solver jobs:

```powershell
Get-ChildItem -LiteralPath 'D:\PHD' -Directory
Test-Path -LiteralPath 'D:\PHD\gmsh-4.15.2-Windows64\gmsh.exe'
Get-Command gmsh,ccx -ErrorAction SilentlyContinue
[Environment]::GetEnvironmentVariable('GMSH_EXE','Process')
[Environment]::GetEnvironmentVariable('CCX_EXE','Process')

@'
import subprocess,tempfile
from pathlib import Path
b=Path(r'D:\PHD\CoupledBeams\CoupledBeams\results\_smoke\3d_fem_environment_check\calculix_2p22\CalculiX-2.22.0-win-x64\bin')
commands=[[r'D:\PHD\gmsh-4.15.2-Windows64\gmsh.exe','--version'],[str(b/'ccx.exe'),'-v'],[str(b/'ccx_MT.exe'),'-v']]
with tempfile.TemporaryDirectory(prefix='fem_version_only_') as cwd:
    for command in commands:
        p=subprocess.run(command,cwd=cwd,capture_output=True,text=True,timeout=10,creationflags=subprocess.CREATE_NO_WINDOW)
        print(command,p.returncode,p.stdout.strip(),p.stderr.strip())
    print('Created files:',[p.name for p in Path(cwd).iterdir()])
'@ | & 'D:\python\Pycharm\pythonProject\.venv\Scripts\python.exe' -

& 'D:\texstudio\miktex\miktex\bin\x64\pdftotext.exe' -layout -f 600 -l 600 `
 'D:\PHD\CalculiX-Windows-master\src\downloads\ccx_2.22.pdf' -
```

Archive catalogs were read via ZIP/tar streams without extraction. Python
package checks used importlib.util.find_spec/importlib.metadata in the stated
interpreter. Verification: safe probes/source and old-artifact inspection,
scoped links, frozen-file/staging preservation, git diff --check. No scientific
pytest/benchmark was run. Existing solvers/tests/results/manifests and nonlinear
energies/equations are unchanged; old D/K records preserved.
[NLSP-K09](../memory/knowledge.md#nlsp-k09) stores the substantive readiness fact.
The audit stops here: no FEM-1, static/transient job, mesh, installation, archive
expansion, persistent PATH change or new solver was authorized or performed.
