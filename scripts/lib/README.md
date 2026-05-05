# Internal helpers

This directory contains reusable helper modules that are not meant to be run directly.

- `tracked_bending_descendant_shapes.py` provides the shared tracked-state extraction, one-case normalization, one-case drawing, and output-path helpers used by both the single-shape and multi-panel tracked bending descendant commands.

Some historical helpers remain at root-level paths, especially `scripts/sweep_grid_policy.py`, because moving them would require broader import updates with no numerical benefit.
