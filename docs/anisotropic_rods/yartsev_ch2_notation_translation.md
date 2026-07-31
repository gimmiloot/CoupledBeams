# Yartsev Chapter-2 notation translation

This note records notation only. It does not add a physical derivation or
modify either existing model.

## State-component translation

| Chapter-2 notation | Old out-of-plane project notation | Exact relation |
| --- | --- | --- |
| `w_i` | `z_i` | `z_i^{old} = w_i` |
| `psi_i` | `varphi_i` | `varphi_i^{old} = -psi_i` |
| `Phi_i` | `psi_i^{old}` | `psi_i^{old} = Phi_i` |
| `Q_i` | `Q_i^{old}` | `Q_i^{old} = Q_i` |
| `M_i` | `M_i^{old}` | `M_i^{old} = -M_i` |
| `M_{T,i}` | `T_i^{old}` | `T_i^{old} = M_{T,i}` |

With the state orders

```text
y_book = [w, psi, Phi, Q, M, M_T]^T,
y_old  = [z, varphi, psi_old, Q_old, M_old, T_old]^T,
```

the translation is

```text
y_old = P y_book,
P = diag(1, -1, 1, 1, -1, 1).
```

For two arms, `P_2 = block_diag(P, P)`.

## Naming rules

- The project bases `e_z`, `t_i`, and `n_i` are unchanged.
- The new Chapter-2 coupled model uses book notation.
- The old isotropic out-of-plane baseline retains its own notation.
- Mixing the two meanings of `psi` without an explicit qualifier is
  prohibited.
- In Python, `psi_b` means the book variable `psi`, while `Phi_t` means the
  book variable `Phi`.

