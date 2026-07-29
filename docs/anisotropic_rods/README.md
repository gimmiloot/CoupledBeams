# Anisotropic Rods — Reserved Research Direction

```text
Status: planned and deferred
Scientific work: not started
Code implementation: not started
Primary source: electronic book to be supplied or identified by the user
Start condition: completion of the current repository refactoring
```

## Purpose

This directory is reserved for a future research direction. It does not contain
an approved physical model, and the term "anisotropic rods" has not yet been
made more specific. The existing isotropic model remains the unchanged
baseline. No current equation is to be generalized before the source book has
been studied and its theory reproduced.

## Deferred Until Refactoring Is Complete

- Reading the electronic book.
- Creating a source note.
- Reproducing the book's derivations and formulas.
- Checking its notation.
- Determining the type of anisotropy.
- Choosing a planar or spatial formulation.
- Choosing generalized variables.
- Formulating assumptions.
- Designing a constitutive API.
- Deriving a coupled determinant.
- Creating tests.
- FEM validation.
- Parameter studies.

## Intended Starting Sequence

1. Complete the current repository refactoring.
2. Obtain the exact book file and its bibliographic data from the user.
3. Register the source in the project's established literature structure.
4. Read the book without immediately changing the existing code.
5. Reproduce the key theory, derivations, and formulas separately from the
   baseline model.
6. Check dimensions, notation, and limiting cases.
7. Compare the book's model with the existing isotropic formulation.
8. Only then create:

   - a scientific research plan;
   - assumptions;
   - notation;
   - a verification contract;
   - a code architecture.

## Preservation Rule

> No existing isotropic equation, determinant entry, sign convention, unknown ordering, root solver or FEM model may be changed merely to make room for the future anisotropic direction.

## What This Directory Currently Contains

Only this reservation/status document.
