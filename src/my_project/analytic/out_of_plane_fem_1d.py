from __future__ import annotations

from dataclasses import dataclass
from math import cos, isfinite, sin

import numpy as np

from .formulas_thickness_mismatch import thickness_mismatch_factors


@dataclass(frozen=True)
class OutOfPlaneFem1DRodProperties:
    """Physical nondimensional properties for one local rod."""

    length: float
    area: float
    bending_inertia: float
    polar_inertia: float
    bending_stiffness: float
    torsional_stiffness: float
    mass_per_length: float
    polar_mass_per_length: float


@dataclass(frozen=True)
class OutOfPlaneFem1DModeResult:
    """Sorted FEM modes in Lambda scale with stiffness-energy fractions."""

    lambdas: np.ndarray
    eigenvalues: np.ndarray
    eigenvectors: np.ndarray
    bending_energy_fraction: np.ndarray
    torsion_energy_fraction: np.ndarray


def _finite_float(name: str, value: float) -> float:
    value_f = float(value)
    if not isfinite(value_f):
        raise ValueError(f"{name} must be finite.")
    return value_f


def _positive_float(name: str, value: float) -> float:
    value_f = _finite_float(name, value)
    if value_f <= 0.0:
        raise ValueError(f"{name} must be positive.")
    return value_f


def _positive_int(name: str, value: int) -> int:
    value_i = int(value)
    if value_i <= 0:
        raise ValueError(f"{name} must be positive.")
    return value_i


def bending_element_matrices_phi(
    length: float,
    bending_stiffness: float,
    mass_per_length: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return EB bending element matrices in DOFs ``[z1, phi1, z2, phi2]``.

    Standard Hermite Euler--Bernoulli matrices are written for
    ``theta = dz/ds``. The out-of-plane determinant uses ``phi = -dz/ds``.
    Therefore ``q_theta = diag(1, -1, 1, -1) q_phi`` and
    ``K_phi = T.T @ K_theta @ T``, with the same transformation for mass.
    """

    h = _positive_float("length", length)
    ej = _positive_float("bending_stiffness", bending_stiffness)
    rho_s = _positive_float("mass_per_length", mass_per_length)
    h2 = h * h

    k_theta = (ej / h**3) * np.array(
        [
            [12.0, 6.0 * h, -12.0, 6.0 * h],
            [6.0 * h, 4.0 * h2, -6.0 * h, 2.0 * h2],
            [-12.0, -6.0 * h, 12.0, -6.0 * h],
            [6.0 * h, 2.0 * h2, -6.0 * h, 4.0 * h2],
        ],
        dtype=float,
    )
    m_theta = (rho_s * h / 420.0) * np.array(
        [
            [156.0, 22.0 * h, 54.0, -13.0 * h],
            [22.0 * h, 4.0 * h2, 13.0 * h, -3.0 * h2],
            [54.0, 13.0 * h, 156.0, -22.0 * h],
            [-13.0 * h, -3.0 * h2, -22.0 * h, 4.0 * h2],
        ],
        dtype=float,
    )
    theta_to_phi = np.diag([1.0, -1.0, 1.0, -1.0])
    return theta_to_phi.T @ k_theta @ theta_to_phi, theta_to_phi.T @ m_theta @ theta_to_phi


def torsion_element_matrices(
    length: float,
    torsional_stiffness: float,
    polar_mass_per_length: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return Saint-Venant torsion matrices in DOFs ``[psi1, psi2]``."""

    h = _positive_float("length", length)
    g_jp = _positive_float("torsional_stiffness", torsional_stiffness)
    rho_jp = _positive_float("polar_mass_per_length", polar_mass_per_length)
    stiffness = (g_jp / h) * np.array([[1.0, -1.0], [-1.0, 1.0]], dtype=float)
    mass = (rho_jp * h / 6.0) * np.array([[2.0, 1.0], [1.0, 2.0]], dtype=float)
    return stiffness, mass


def out_of_plane_fem_1d_rod_properties(
    mu: float,
    epsilon: float,
    eta: float = 0.0,
    poisson: float = 0.3,
) -> tuple[OutOfPlaneFem1DRodProperties, OutOfPlaneFem1DRodProperties]:
    """Return nondimensional rod properties used by the 1D FEM validation.

    The scaling is ``l0=E=rho=S0=1``, ``J0=epsilon**2``, ``Jp0=2*J0`` and
    ``G=E/(2*(1+poisson))``. The generalized eigenproblem is
    ``K q = omega**2 M q`` and is converted to the analytic frequency
    parameter by ``Lambda = (omega**2 / epsilon**2)**0.25``.
    """

    epsilon_f = _positive_float("epsilon", epsilon)
    poisson_f = _finite_float("poisson", poisson)
    if poisson_f <= -1.0:
        raise ValueError("poisson must be greater than -1.")
    factors = thickness_mismatch_factors(mu, eta)
    for name, value in (("tau1", factors.tau1), ("tau2", factors.tau2)):
        if not isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be positive and finite.")

    shear_modulus = 1.0 / (2.0 * (1.0 + poisson_f))
    base_bending_inertia = epsilon_f**2
    base_polar_inertia = 2.0 * base_bending_inertia

    def rod(length: float, tau: float) -> OutOfPlaneFem1DRodProperties:
        area = tau**2
        bending_inertia = base_bending_inertia * tau**4
        polar_inertia = base_polar_inertia * tau**4
        return OutOfPlaneFem1DRodProperties(
            length=float(length),
            area=float(area),
            bending_inertia=float(bending_inertia),
            polar_inertia=float(polar_inertia),
            bending_stiffness=float(bending_inertia),
            torsional_stiffness=float(shear_modulus * polar_inertia),
            mass_per_length=float(area),
            polar_mass_per_length=float(polar_inertia),
        )

    return rod(1.0 - factors.mu, factors.tau1), rod(1.0 + factors.mu, factors.tau2)


def _node_dof(rod_index: int, node_index: int, component_index: int, nodes_per_rod: int) -> int:
    rod_offset = rod_index * nodes_per_rod * 3
    return rod_offset + 3 * node_index + component_index


def _assemble_full_split_matrices(
    *,
    mu: float,
    epsilon: float,
    eta: float,
    poisson: float,
    n_elements_per_rod: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rods = out_of_plane_fem_1d_rod_properties(mu=mu, epsilon=epsilon, eta=eta, poisson=poisson)
    n_elem = _positive_int("n_elements_per_rod", n_elements_per_rod)
    nodes_per_rod = n_elem + 1
    ndof = 2 * nodes_per_rod * 3
    bending_stiffness = np.zeros((ndof, ndof), dtype=float)
    torsion_stiffness = np.zeros((ndof, ndof), dtype=float)
    mass = np.zeros((ndof, ndof), dtype=float)

    def assemble(dofs: list[int], local_k_b: np.ndarray, local_k_t: np.ndarray, local_m: np.ndarray) -> None:
        for row_index, row_dof in enumerate(dofs):
            for col_index, col_dof in enumerate(dofs):
                bending_stiffness[row_dof, col_dof] += local_k_b[row_index, col_index]
                torsion_stiffness[row_dof, col_dof] += local_k_t[row_index, col_index]
                mass[row_dof, col_dof] += local_m[row_index, col_index]

    for rod_index, rod in enumerate(rods):
        h = rod.length / n_elem
        k_b, m_b = bending_element_matrices_phi(h, rod.bending_stiffness, rod.mass_per_length)
        k_t, m_t = torsion_element_matrices(h, rod.torsional_stiffness, rod.polar_mass_per_length)
        local_k_b = np.zeros((6, 6), dtype=float)
        local_k_t = np.zeros((6, 6), dtype=float)
        local_m = np.zeros((6, 6), dtype=float)
        for local_row, row in enumerate([0, 1, 3, 4]):
            for local_col, col in enumerate([0, 1, 3, 4]):
                local_k_b[row, col] += k_b[local_row, local_col]
                local_m[row, col] += m_b[local_row, local_col]
        for local_row, row in enumerate([2, 5]):
            for local_col, col in enumerate([2, 5]):
                local_k_t[row, col] += k_t[local_row, local_col]
                local_m[row, col] += m_t[local_row, local_col]

        for elem in range(n_elem):
            dofs = [
                _node_dof(rod_index, elem, 0, nodes_per_rod),
                _node_dof(rod_index, elem, 1, nodes_per_rod),
                _node_dof(rod_index, elem, 2, nodes_per_rod),
                _node_dof(rod_index, elem + 1, 0, nodes_per_rod),
                _node_dof(rod_index, elem + 1, 1, nodes_per_rod),
                _node_dof(rod_index, elem + 1, 2, nodes_per_rod),
            ]
            assemble(dofs, local_k_b, local_k_t, local_m)

    return bending_stiffness, torsion_stiffness, mass


def _assemble_full_matrices(
    *,
    mu: float,
    epsilon: float,
    eta: float,
    poisson: float,
    n_elements_per_rod: int,
) -> tuple[np.ndarray, np.ndarray]:
    bending_stiffness, torsion_stiffness, mass = _assemble_full_split_matrices(
        mu=mu,
        epsilon=epsilon,
        eta=eta,
        poisson=poisson,
        n_elements_per_rod=n_elements_per_rod,
    )
    return bending_stiffness + torsion_stiffness, mass


def _reduction_matrix(
    *,
    beta: float,
    n_elements_per_rod: int,
) -> tuple[np.ndarray, list[tuple[str, int, str]]]:
    beta_f = _finite_float("beta", beta)
    n_elem = _positive_int("n_elements_per_rod", n_elements_per_rod)
    nodes_per_rod = n_elem + 1
    full_ndof = 2 * nodes_per_rod * 3
    labels: list[tuple[str, int, str]] = []

    for rod_name in ("rod1", "rod2"):
        for node in range(1, n_elem):
            for component in ("z", "phi", "psi"):
                labels.append((rod_name, node, component))
    joint_z_col = len(labels)
    labels.append(("joint", n_elem, "z"))
    joint_psi_col = len(labels)
    labels.append(("joint", n_elem, "psi"))
    joint_phi_col = len(labels)
    labels.append(("joint", n_elem, "phi"))

    reduction = np.zeros((full_ndof, len(labels)), dtype=float)
    component_index = {"z": 0, "phi": 1, "psi": 2}

    col = 0
    for rod_index in (0, 1):
        for node in range(1, n_elem):
            for component in ("z", "phi", "psi"):
                row = _node_dof(rod_index, node, component_index[component], nodes_per_rod)
                reduction[row, col] = 1.0
                col += 1

    joint_node = n_elem
    # Rod 1 local joint variables equal the global reduced joint variables.
    reduction[_node_dof(0, joint_node, 0, nodes_per_rod), joint_z_col] = 1.0
    reduction[_node_dof(0, joint_node, 1, nodes_per_rod), joint_phi_col] = 1.0
    reduction[_node_dof(0, joint_node, 2, nodes_per_rod), joint_psi_col] = 1.0

    # Rod 2 local variables from theta_1 = theta_2 in the adopted local bases.
    cb = cos(beta_f)
    sb = sin(beta_f)
    reduction[_node_dof(1, joint_node, 0, nodes_per_rod), joint_z_col] = 1.0
    reduction[_node_dof(1, joint_node, 1, nodes_per_rod), joint_psi_col] = -sb
    reduction[_node_dof(1, joint_node, 1, nodes_per_rod), joint_phi_col] = -cb
    reduction[_node_dof(1, joint_node, 2, nodes_per_rod), joint_psi_col] = -cb
    reduction[_node_dof(1, joint_node, 2, nodes_per_rod), joint_phi_col] = sb

    return reduction, labels


def _assemble_out_of_plane_fem_1d_matrices_with_labels(
    beta: float,
    mu: float,
    epsilon: float,
    eta: float = 0.0,
    poisson: float = 0.3,
    n_elements_per_rod: int = 16,
) -> tuple[np.ndarray, np.ndarray, list[tuple[str, int, str]]]:
    full_k, full_m = _assemble_full_matrices(
        mu=mu,
        epsilon=epsilon,
        eta=eta,
        poisson=poisson,
        n_elements_per_rod=n_elements_per_rod,
    )
    reduction, labels = _reduction_matrix(beta=beta, n_elements_per_rod=n_elements_per_rod)
    reduced_k = reduction.T @ full_k @ reduction
    reduced_m = reduction.T @ full_m @ reduction
    return reduced_k, reduced_m, labels


def assemble_out_of_plane_fem_1d_matrices(
    beta: float,
    mu: float,
    epsilon: float,
    eta: float = 0.0,
    poisson: float = 0.3,
    n_elements_per_rod: int = 16,
) -> tuple[np.ndarray, np.ndarray]:
    """Return reduced 1D EB+torsion FEM matrices after clamps and joint MPCs."""

    stiffness, mass, _labels = _assemble_out_of_plane_fem_1d_matrices_with_labels(
        beta=beta,
        mu=mu,
        epsilon=epsilon,
        eta=eta,
        poisson=poisson,
        n_elements_per_rod=n_elements_per_rod,
    )
    warnings = out_of_plane_fem_1d_matrix_warnings(stiffness, mass)
    if warnings:
        raise ValueError("; ".join(warnings))
    return stiffness, mass


def assemble_out_of_plane_fem_1d_energy_matrices(
    beta: float,
    mu: float,
    epsilon: float,
    eta: float = 0.0,
    poisson: float = 0.3,
    n_elements_per_rod: int = 16,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return reduced bending stiffness, torsion stiffness, and total mass."""

    full_k_b, full_k_t, full_m = _assemble_full_split_matrices(
        mu=mu,
        epsilon=epsilon,
        eta=eta,
        poisson=poisson,
        n_elements_per_rod=n_elements_per_rod,
    )
    reduction, _labels = _reduction_matrix(beta=beta, n_elements_per_rod=n_elements_per_rod)
    k_b = reduction.T @ full_k_b @ reduction
    k_t = reduction.T @ full_k_t @ reduction
    mass = reduction.T @ full_m @ reduction
    warnings = out_of_plane_fem_1d_matrix_warnings(k_b + k_t, mass)
    if warnings:
        raise ValueError("; ".join(warnings))
    return k_b, k_t, mass


def assemble_out_of_plane_fem_1d_beta0_blocks(
    mu: float,
    epsilon: float,
    eta: float = 0.0,
    poisson: float = 0.3,
    n_elements_per_rod: int = 16,
) -> tuple[tuple[np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray]]:
    """Return beta=0 bending and torsion blocks for block-level checks."""

    stiffness, mass, labels = _assemble_out_of_plane_fem_1d_matrices_with_labels(
        beta=0.0,
        mu=mu,
        epsilon=epsilon,
        eta=eta,
        poisson=poisson,
        n_elements_per_rod=n_elements_per_rod,
    )
    bending = [index for index, label in enumerate(labels) if label[2] in {"z", "phi"}]
    torsion = [index for index, label in enumerate(labels) if label[2] == "psi"]
    k_b = stiffness[np.ix_(bending, bending)]
    m_b = mass[np.ix_(bending, bending)]
    k_t = stiffness[np.ix_(torsion, torsion)]
    m_t = mass[np.ix_(torsion, torsion)]
    return (k_b, m_b), (k_t, m_t)


def out_of_plane_fem_1d_matrix_warnings(stiffness: np.ndarray, mass: np.ndarray) -> list[str]:
    """Return symmetry and definiteness warnings for reduced FEM matrices."""

    stiffness_a = np.asarray(stiffness, dtype=float)
    mass_a = np.asarray(mass, dtype=float)
    warnings: list[str] = []
    if stiffness_a.shape != mass_a.shape or stiffness_a.ndim != 2 or stiffness_a.shape[0] != stiffness_a.shape[1]:
        return ["K and M must be square matrices with identical shape."]
    if not np.all(np.isfinite(stiffness_a)):
        warnings.append("K contains non-finite entries.")
    if not np.all(np.isfinite(mass_a)):
        warnings.append("M contains non-finite entries.")
    if not np.allclose(stiffness_a, stiffness_a.T, rtol=1e-11, atol=1e-12):
        warnings.append("K is not symmetric.")
    if not np.allclose(mass_a, mass_a.T, rtol=1e-11, atol=1e-12):
        warnings.append("M is not symmetric.")
    if warnings:
        return warnings

    try:
        mass_eigs = np.linalg.eigvalsh(mass_a)
    except np.linalg.LinAlgError as exc:
        warnings.append(f"M eigensolve failed: {exc}")
        return warnings
    if float(np.min(mass_eigs)) <= 0.0:
        warnings.append("M is not positive definite.")

    try:
        stiffness_eigs = np.linalg.eigvalsh(stiffness_a)
    except np.linalg.LinAlgError as exc:
        warnings.append(f"K eigensolve failed: {exc}")
        return warnings
    if float(np.min(stiffness_eigs)) <= -1e-10:
        warnings.append("K has a negative eigenvalue.")
    return warnings


def out_of_plane_fem_1d_energy_fractions(
    bending_stiffness: np.ndarray,
    torsion_stiffness: np.ndarray,
    eigenvectors: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return bending and torsion stiffness-energy fractions for modes."""

    k_b = np.asarray(bending_stiffness, dtype=float)
    k_t = np.asarray(torsion_stiffness, dtype=float)
    vectors = np.asarray(eigenvectors, dtype=float)
    bending: list[float] = []
    torsion: list[float] = []
    for mode_index in range(vectors.shape[1]):
        vector = vectors[:, mode_index]
        e_b = max(0.0, float(vector @ (k_b @ vector)))
        e_t = max(0.0, float(vector @ (k_t @ vector)))
        total = e_b + e_t
        if total <= 0.0:
            bending.append(float("nan"))
            torsion.append(float("nan"))
        else:
            bending.append(e_b / total)
            torsion.append(e_t / total)
    return np.asarray(bending, dtype=float), np.asarray(torsion, dtype=float)


def solve_out_of_plane_fem_1d_modes(
    beta: float,
    mu: float,
    epsilon: float,
    eta: float = 0.0,
    poisson: float = 0.3,
    n_elements_per_rod: int = 16,
    n_modes: int = 8,
) -> OutOfPlaneFem1DModeResult:
    """Return sorted positive FEM modes with Lambda conversion and energy fractions."""

    try:
        from scipy.linalg import eigh
    except ImportError as exc:
        raise ImportError("scipy.linalg.eigh is required for the 1D out-of-plane FEM solver.") from exc

    n_modes_i = _positive_int("n_modes", n_modes)
    epsilon_f = _positive_float("epsilon", epsilon)
    k_b, k_t, mass = assemble_out_of_plane_fem_1d_energy_matrices(
        beta=beta,
        mu=mu,
        epsilon=epsilon_f,
        eta=eta,
        poisson=poisson,
        n_elements_per_rod=n_elements_per_rod,
    )
    stiffness = k_b + k_t
    max_modes = stiffness.shape[0]
    if n_modes_i > max_modes:
        raise ValueError(f"n_modes={n_modes_i} exceeds reduced matrix size {max_modes}.")
    eigenvalues, eigenvectors = eigh(stiffness, mass, subset_by_index=[0, n_modes_i - 1])
    eigenvalues = np.asarray(eigenvalues, dtype=float)
    if not np.all(np.isfinite(eigenvalues)):
        raise ValueError("FEM eigenvalues must be finite.")
    if np.any(eigenvalues <= 0.0):
        raise ValueError("FEM eigenvalues must be positive.")
    lambdas = (eigenvalues / epsilon_f**2) ** 0.25
    if not np.all(np.isfinite(lambdas)) or np.any(lambdas <= 0.0):
        raise ValueError("Converted FEM Lambda values must be finite and positive.")
    bending, torsion = out_of_plane_fem_1d_energy_fractions(k_b, k_t, eigenvectors)
    return OutOfPlaneFem1DModeResult(
        lambdas=np.asarray(lambdas, dtype=float),
        eigenvalues=eigenvalues,
        eigenvectors=np.asarray(eigenvectors, dtype=float),
        bending_energy_fraction=bending,
        torsion_energy_fraction=torsion,
    )


def solve_out_of_plane_fem_1d_frequencies(
    beta: float,
    mu: float,
    epsilon: float,
    eta: float = 0.0,
    poisson: float = 0.3,
    n_elements_per_rod: int = 16,
    n_modes: int = 8,
) -> np.ndarray:
    """Return sorted positive FEM frequencies in the analytic ``Lambda`` scale."""

    return solve_out_of_plane_fem_1d_modes(
        beta=beta,
        mu=mu,
        epsilon=epsilon,
        eta=eta,
        poisson=poisson,
        n_elements_per_rod=n_elements_per_rod,
        n_modes=n_modes,
    ).lambdas


def first_beta0_eta0_torsion_fem_root(
    epsilon: float,
    poisson: float = 0.3,
    n_elements_per_rod: int = 32,
) -> float:
    """Return the first beta=0, eta=0 torsion-block FEM root in Lambda scale."""

    try:
        from scipy.linalg import eigh
    except ImportError as exc:
        raise ImportError("scipy.linalg.eigh is required for the 1D out-of-plane FEM solver.") from exc

    epsilon_f = _positive_float("epsilon", epsilon)
    (_k_b, _m_b), (k_t, m_t) = assemble_out_of_plane_fem_1d_beta0_blocks(
        mu=0.0,
        epsilon=epsilon_f,
        eta=0.0,
        poisson=poisson,
        n_elements_per_rod=n_elements_per_rod,
    )
    eigenvalue = float(eigh(k_t, m_t, subset_by_index=[0, 0], eigvals_only=True)[0])
    return float((eigenvalue / epsilon_f**2) ** 0.25)


__all__ = [
    "OutOfPlaneFem1DModeResult",
    "OutOfPlaneFem1DRodProperties",
    "assemble_out_of_plane_fem_1d_beta0_blocks",
    "assemble_out_of_plane_fem_1d_energy_matrices",
    "assemble_out_of_plane_fem_1d_matrices",
    "bending_element_matrices_phi",
    "first_beta0_eta0_torsion_fem_root",
    "out_of_plane_fem_1d_energy_fractions",
    "out_of_plane_fem_1d_matrix_warnings",
    "out_of_plane_fem_1d_rod_properties",
    "solve_out_of_plane_fem_1d_frequencies",
    "solve_out_of_plane_fem_1d_modes",
    "torsion_element_matrices",
]
