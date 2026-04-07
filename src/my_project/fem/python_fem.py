"""
МКЭ-верификация спектра связанных стержней
===========================================
Два стержня Эйлера–Бернулли, соединённых ЖЁСТКИМ СТЫКОМ под углом beta.
Оба внешних конца — заделка (CC).

Физическая модель:
    Плечо 1 (горизонтальное) --- жёсткий стык (угол beta) --- Плечо 2

    Стык передаёт ВСЕ усилия и моменты.
    В узле стыка общие 3 DOF: (u, v, theta).
    Плечо 2 повёрнуто на beta => осевая-изгибная связь через sin/cos(beta).

Тип элемента:
    Комбинированный: осевой стержень + балка Эйлера–Бернулли.
    6 DOF на элемент: [u1, v1, theta1, u2, v2, theta2]
    Согласованная матрица масс (consistent) для обеих частей.

Трансформация плеча 2:
    Каждый элемент плеча 2 поворачивается из локальной в глобальную СК:
        T = block_diag(R, R),  R = [[cos, -sin, 0], [sin, cos, 0], [0, 0, 1]]
        K_глоб = T^T @ K_лок @ T
        M_глоб = T^T @ M_лок @ T

Отслеживание мод:
    MAC (Modal Assurance Criterion) + венгерский алгоритм
    (scipy.optimize.linear_sum_assignment) для оптимального сопоставления.
    Фиксированная сетка (N_ELEM на плечо) => DOF = const => MAC всегда работает.

Параметры задачи:
    r      = 0.005 м          радиус кругового сечения
    E      = 2.1e11 Па        сталь
    rho    = 7800 кг/м³
    L_tot  = 2.0 м            L1 + L2
    ell    = L_tot / 2 = 1 м  характерная полудлина
    beta   = 15°               угол сопряжения
    eps    = r / (2*ell) = 0.0025   параметр осевой-изгибной связи

Безразмерная формулировка:
    EI = 1,  rhoA = 1,  ell = 1
    EA = 1 / eps^2 = 160000
    f~ = omega * ell^2 * sqrt(rhoA / EI)
    f [Гц] = f~ * scale,  scale = sqrt(EI_разм / rhoA_разм) / (2*pi*ell^2)

Параметр асимметрии:
    mu = (L2 - L1) / (L1 + L2)
    L1 = ell * (1 - mu),  L2 = ell * (1 + mu)

Выход:
    fem_spectrum.csv — частоты в Гц и безразмерные, референсы FP/FF

Автор: Нестерчук Г.А., СПбГУ, кафедра ТиПМ
"""

from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import linear_sum_assignment


# ============================================================
# Физические параметры
# ============================================================
r     = 0.005           # радиус сечения [м]
E     = 2.1e11          # модуль Юнга [Па]
rho   = 7800.0          # плотность [кг/м³]
L_tot = 2.0             # полная длина L1 + L2 [м]
ell   = L_tot / 2       # характерная полудлина [м]
BETA  = 15.0            # угол сопряжения [градусы]

A     = np.pi * r**2                    # площадь сечения
I     = np.pi * r**4 / 4                # момент инерции
EI    = E * I                           # изгибная жёсткость
EA    = E * A                           # осевая жёсткость
rhoA  = rho * A                         # погонная масса
eps   = np.sqrt(I / A) / ell            # = r / (2*ell), параметр связи
scale = np.sqrt(EI / rhoA) / (2 * np.pi * ell**2)  # Гц на единицу f~


# ============================================================
# Безразмерные параметры МКЭ
# ============================================================
EI_nd   = 1.0
rhoA_nd = 1.0
EA_nd   = 1.0 / eps**2   # = 160000
N_ELEM  = 80             # элементов на плечо (фиксировано для MAC)


# ============================================================
# Элементные матрицы: осевой стержень + балка ЭБ (6 DOF)
# Порядок DOF: [u1, v1, theta1, u2, v2, theta2]
# ============================================================
def elem_K(Le):
    """Матрица жёсткости элемента (6×6): осевая + изгибная части."""
    K = np.zeros((6, 6))

    # --- Осевая часть: DOF 0 (u1) и 3 (u2) ---
    ea = EA_nd / Le
    K[0, 0] += ea;  K[0, 3] -= ea
    K[3, 0] -= ea;  K[3, 3] += ea

    # --- Изгибная часть: DOF 1,2,4,5 = [v1, theta1, v2, theta2] ---
    c = EI_nd / Le**3
    Kb = c * np.array([
        [ 12,      6*Le,   -12,      6*Le  ],
        [ 6*Le,    4*Le**2, -6*Le,   2*Le**2],
        [-12,     -6*Le,    12,     -6*Le   ],
        [ 6*Le,    2*Le**2, -6*Le,   4*Le**2]])
    for i, di in enumerate([1, 2, 4, 5]):
        for j, dj in enumerate([1, 2, 4, 5]):
            K[di, dj] += Kb[i, j]

    return K


def elem_M(Le):
    """Согласованная матрица масс элемента (6×6): осевая + изгибная."""
    M = np.zeros((6, 6))

    # --- Осевая масса (consistent bar) ---
    ma = rhoA_nd * Le / 6
    M[0, 0] += 2*ma;  M[0, 3] += ma
    M[3, 0] += ma;    M[3, 3] += 2*ma

    # --- Изгибная масса (consistent beam) ---
    c = rhoA_nd * Le / 420
    Mb = c * np.array([
        [ 156,      22*Le,    54,     -13*Le  ],
        [ 22*Le,     4*Le**2, 13*Le,   -3*Le**2],
        [  54,      13*Le,   156,     -22*Le   ],
        [-13*Le,    -3*Le**2,-22*Le,    4*Le**2]])
    for i, di in enumerate([1, 2, 4, 5]):
        for j, dj in enumerate([1, 2, 4, 5]):
            M[di, dj] += Mb[i, j]

    return M


def rotation_matrix_6x6(beta_rad):
    """
    Матрица поворота для одного элемента (2 узла, 6 DOF).
    Переводит из локальной СК плеча 2 в глобальную.

    Для каждого узла: [u_лок, v_лок, theta_лок] -> [u_глоб, v_глоб, theta_глоб]
        R = [[cos, -sin, 0],
             [sin,  cos, 0],
             [0,    0,   1]]

    Для элемента: T = block_diag(R, R)
    """
    c, s = np.cos(beta_rad), np.sin(beta_rad)
    R = np.array([[c, -s, 0],
                  [s,  c, 0],
                  [0,  0, 1]])
    T = np.zeros((6, 6))
    T[:3, :3] = R
    T[3:, 3:] = R
    return T


# ============================================================
# Сборка глобальных матриц и решение задачи на собственные значения
# ============================================================
def fem_solve(mu, beta_deg=0.0, n_modes=14):
    """
    Решает обобщённую задачу на собственные значения K*phi = omega^2*M*phi
    для двух связанных стержней под углом beta.

    Параметры
    ---------
    mu : float
        Параметр асимметрии: mu = (L2 - L1) / (L1 + L2).
        L1 = ell*(1 - mu),  L2 = ell*(1 + mu).
    beta_deg : float
        Угол сопряжения в градусах.
    n_modes : int
        Число вычисляемых собственных значений.

    Возвращает
    ----------
    omega : ndarray, shape (n_modes,)
        Собственные частоты (безразмерные).
    vecs : ndarray, shape (n_free, n_modes)
        Собственные вектора на свободных DOF.

    Нумерация DOF (жёсткий стык, все 3 DOF общие):
        Плечо 1, узел i:      [3*i, 3*i+1, 3*i+2] = [u, v, theta]
        Стык (узел n1):        [3*n1, 3*n1+1, 3*n1+2] = общие [U, V, theta]
        Плечо 2, узел j>=1:   [3*(n1+j), 3*(n1+j)+1, 3*(n1+j)+2]
        Итого DOF = 3*(n1 + n2 + 1)
    """
    L1 = max(ell * (1 - mu), 1e-3)
    L2 = max(ell * (1 + mu), 1e-3)

    n1 = n2 = N_ELEM
    le1, le2 = L1 / n1, L2 / n2
    T6 = rotation_matrix_6x6(np.radians(beta_deg))

    ndof = 3 * (n1 + n2 + 1)
    K = np.zeros((ndof, ndof))
    M = np.zeros((ndof, ndof))

    def assemble(dof_map, Ke, Me):
        for r_idx, dr in enumerate(dof_map):
            for c_idx, dc in enumerate(dof_map):
                K[dr, dc] += Ke[r_idx, c_idx]
                M[dr, dc] += Me[r_idx, c_idx]

    # --- Плечо 1: горизонтальное, без поворота ---
    for i in range(n1):
        dofs = [3*i, 3*i+1, 3*i+2, 3*(i+1), 3*(i+1)+1, 3*(i+1)+2]
        assemble(dofs, elem_K(le1), elem_M(le1))

    # --- Плечо 2: повёрнуто на beta, каждый элемент трансформируется ---
    for j in range(n2):
        # Узел j плеча 2 в глобальной нумерации: узел (n1 + j)
        dofs = [3*(n1+j),   3*(n1+j)+1,   3*(n1+j)+2,
                3*(n1+j+1), 3*(n1+j+1)+1, 3*(n1+j+1)+2]
        Ke_glob = T6.T @ elem_K(le2) @ T6
        Me_glob = T6.T @ elem_M(le2) @ T6
        assemble(dofs, Ke_glob, Me_glob)

    # --- Граничные условия: заделки на обоих концах ---
    bc_left  = [0, 1, 2]                                  # плечо 1, узел 0
    bc_right = [3*(n1+n2), 3*(n1+n2)+1, 3*(n1+n2)+2]     # плечо 2, последний узел
    bc = set(bc_left + bc_right)
    free = sorted(set(range(ndof)) - bc)

    Kf = K[np.ix_(free, free)]
    Mf = M[np.ix_(free, free)]
    eigenvalues, eigenvectors = eigh(Kf, Mf, subset_by_index=[0, n_modes - 1])
    omega = np.sqrt(np.maximum(eigenvalues, 0))
    return omega, eigenvectors


# ============================================================
# Отслеживание мод через MAC (венгерский алгоритм)
# ============================================================
def mac_value(phi_a, phi_b):
    """MAC (Modal Assurance Criterion) между двумя формами колебаний."""
    numerator = abs(np.dot(phi_a, phi_b))**2
    denominator = np.dot(phi_a, phi_a) * np.dot(phi_b, phi_b)
    return numerator / denominator if denominator > 0 else 0.0


def track_modes(mu_array, beta_deg=0.0, n_track=8):
    """
    Вычисляет собственные частоты по диапазону mu с отслеживанием мод
    через MAC для гладких кривых через veering-области.

    Использует венгерский алгоритм (scipy.optimize.linear_sum_assignment)
    для оптимального сопоставления мод на каждом шаге.

    Параметры
    ---------
    mu_array : ndarray
        Значения параметра асимметрии.
    beta_deg : float
        Угол сопряжения в градусах.
    n_track : int
        Число отслеживаемых мод.

    Возвращает
    ----------
    result : ndarray, shape (len(mu_array), n_track)
        Отслеженные безразмерные собственные частоты.
    """
    n_mu = len(mu_array)
    result = np.zeros((n_mu, n_track))
    n_extra = 5  # запас кандидатов для устойчивости

    # --- Первый шаг: инициализация референсных мод ---
    omega_0, vecs_0 = fem_solve(mu_array[0], beta_deg)

    # Пропускаем околонулевые моды (ищем первую выше порога)
    skip = 0
    for k in range(len(omega_0)):
        if omega_0[k] > 5.0:
            skip = k
            break

    result[0] = omega_0[skip : skip + n_track]
    prev_vecs = vecs_0[:, skip : skip + n_track]

    # --- Проход по mu с MAC-отслеживанием ---
    for step in range(1, n_mu):
        omega_k, vecs_k = fem_solve(mu_array[step], beta_deg)
        candidates_omega = omega_k[skip:]
        candidates_vecs  = vecs_k[:, skip:]
        n_cand = min(candidates_vecs.shape[1], n_track + n_extra)

        # Матрица MAC: mac_mat[i, j] = MAC(пред_мода_i, канд_мода_j)
        mac_mat = np.array([
            [mac_value(prev_vecs[:, i], candidates_vecs[:, j])
             for j in range(n_cand)]
            for i in range(n_track)
        ])

        # Оптимальное назначение: минимизация суммарного (1 - MAC)
        _, col_assignment = linear_sum_assignment(1.0 - mac_mat)

        result[step] = candidates_omega[col_assignment]
        prev_vecs = candidates_vecs[:, col_assignment]

    return result


# ============================================================
# Референсные частоты одиночных балок
# ============================================================
# Характеристические волновые числа:
#   FP (заделка–шарнир): tan(lambda) = tanh(lambda)
#   FF (заделка–заделка): cos(lambda)*cosh(lambda) = 1
LAMBDA_FP = [3.9266, 7.0686, 10.2102, 13.3518, 16.4934]
LAMBDA_FF = [4.7300, 7.8532, 10.9956, 14.1372, 17.2788]


def reference_frequencies_hz(mu_array, lambdas, arm='right'):
    """
    Частоты одиночной балки (FP или FF) в зависимости от mu.

    arm='right': L = ell*(1 + mu)  (длинное плечо)
    arm='left':  L = ell*(1 - mu)  (короткое плечо)
    """
    return np.array([
        [(lam / ((1 + mu) if arm == 'right' else (1 - mu)))**2 * scale
         for lam in lambdas]
        for mu in mu_array
    ])


# ============================================================
# Основной расчёт
# ============================================================
if __name__ == "__main__":
    repo_root = Path(__file__).resolve().parents[3]
    results_dir = repo_root / "results"
    results_dir.mkdir(exist_ok=True)
    console_rule = "-" * 70


    mu_values = np.linspace(0.0, 0.9, 200)
    N_MODES_TRACK = 6

    # --- Спектр для beta=0 (справочно) и beta=15 (основной результат) ---
    print(f"Расчёт МКЭ (N_ELEM={N_ELEM}, DOF={3*(2*N_ELEM+1)})...")
    print(f"  beta = 0 deg  ...", end=" ", flush=True)
    spectrum_beta0  = track_modes(mu_values, beta_deg=0.0,  n_track=N_MODES_TRACK + 2)
    print("готово")

    print(f"  beta = {BETA} deg ...", end=" ", flush=True)
    spectrum_beta15 = track_modes(mu_values, beta_deg=BETA, n_track=N_MODES_TRACK + 2)
    print("готово")

    # --- Перевод в Гц ---
    hz_beta0  = spectrum_beta0[:, :N_MODES_TRACK]  * scale
    hz_beta15 = spectrum_beta15[:, :N_MODES_TRACK] * scale

    # --- Таблица при mu=0 ---
    print(f"\n{'='*70}")
    print(f"  Частоты при mu=0  (r={r*1e3:.0f} мм, L={L_tot} м, beta={BETA} deg)")
    print(f"{'='*70}")
    print(f"  {'Мода':>5}  {'f~ (b=0)':>10}  {'f~ (b=15)':>10}  "
          f"{'Гц (b=0)':>10}  {'Гц (b=15)':>10}")
    print(console_rule)
    for i in range(N_MODES_TRACK):
        print(f"  {i+1:>5}  {spectrum_beta0[0,i]:>10.3f}  {spectrum_beta15[0,i]:>10.3f}  "
              f"{hz_beta0[0,i]:>10.3f}  {hz_beta15[0,i]:>10.3f}")
    print(f"{'='*70}")
    print(f"  Масштаб: 1 ед. f~ = {scale:.5f} Гц")

    # --- Референсные кривые ---
    fp_right = reference_frequencies_hz(mu_values, LAMBDA_FP, 'right')
    ff_right = reference_frequencies_hz(mu_values, LAMBDA_FF, 'right')
    fp_left  = reference_frequencies_hz(mu_values, LAMBDA_FP, 'left')
    ff_left  = reference_frequencies_hz(mu_values, LAMBDA_FF, 'left')

    # --- Сохранение CSV ---
    header = (['mu']
              + [f'f{i+1}_Hz_beta0'  for i in range(N_MODES_TRACK)]
              + [f'f{i+1}_Hz_beta15' for i in range(N_MODES_TRACK)]
              + [f'f{i+1}_nd_beta0'  for i in range(N_MODES_TRACK)]
              + [f'f{i+1}_nd_beta15' for i in range(N_MODES_TRACK)]
              + [f'FP{i+1}_Hz_Lright' for i in range(5)]
              + [f'FF{i+1}_Hz_Lright' for i in range(5)]
              + [f'FP{i+1}_Hz_Lleft'  for i in range(5)]
              + [f'FF{i+1}_Hz_Lleft'  for i in range(5)])

    data = np.column_stack([
        mu_values,
        hz_beta0, hz_beta15,
        spectrum_beta0[:, :N_MODES_TRACK], spectrum_beta15[:, :N_MODES_TRACK],
        fp_right, ff_right, fp_left, ff_left,
    ])

    csv_path = results_dir / 'fem_spectrum.csv'
    np.savetxt(csv_path, data, delimiter=',',
               header=','.join(header), comments='', fmt='%.6f')
    print(f"\n  CSV: {csv_path}  ({len(mu_values)} строк x {data.shape[1]} столбцов)")
