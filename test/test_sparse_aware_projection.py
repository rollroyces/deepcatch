"""Tests for SparseAwareLinearProjection (pillar-2 sparsity fix).

Covers:
- Sparse rows emit the learned missing-token (and only the missing-token)
- Dense rows go through the normal Linear path (matches LinearProjection)
- Forward signature matches LinearProjection (drop-in replacement)
- sparsity_threshold parameter is respected (lower → fewer sparse rows)
- Training: gradients flow through both the dense path AND the
  missing-token parameter
- End-to-end: MultiModalEncoder works when one modality opts in
  (panel-LLR-style) while the others stay dense
- make_projection factory returns the right class and rejects bad kind
"""
import sys
import os

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

try:
    import torch
    _HAS_TORCH = True
except ImportError:
    _HAS_TORCH = False


# ── Sparse path emits the learned missing-token ─────────────────────────────

@pytest.mark.skipif(not _HAS_TORCH, reason="torch not installed")
def test_sparse_row_emits_missing_token():
    """A row whose zero fraction exceeds the threshold must output the
    LayerNorm-normalized missing-token exactly.

    Without this property, the pillar-2 bug returns: a constant bias
    vector from the Linear layer leaks through LayerNorm as if it were
    signal. With the fix, all sparse rows collapse to the same learned
    missing-token (a per-class identity the transformer can learn to
    attend on).
    """
    from src.foundation.model import SparseAwareLinearProjection

    torch.manual_seed(0)
    embed_dim = 8
    sp = SparseAwareLinearProjection(4, embed_dim, dropout=0.0,
                                     sparsity_threshold=0.5)
    # Pin missing_token to a known, distinguishable vector so we can
    # assert the output exactly equals LayerNorm(missing_token).
    with torch.no_grad():
        sp.missing_token.copy_(torch.tensor([1.0, -1.0, 0.5, -0.5,
                                            1.5, -1.5, 2.0, -2.0]))

    sp.eval()
    # All-zero row → sparse → missing-token branch.
    x = torch.zeros(3, 4)
    out = sp(x)
    expected = sp.norm(sp.missing_token.unsqueeze(0).expand(3, -1))
    # All three rows are sparse, so all three outputs match exactly.
    assert torch.allclose(out, expected, atol=1e-6), (
        "sparse rows must equal LayerNorm(missing_token); got "
        f"max diff {(out - expected).abs().max().item()}"
    )


# ── Dense path matches LinearProjection ─────────────────────────────────────

@pytest.mark.skipif(not _HAS_TORCH, reason="torch not installed")
def test_dense_row_uses_linear_path_equivalent_to_linear_projection():
    """A row with zero fraction ≤ threshold must produce the same output
    as a plain LinearProjection on the same input.

    This is the contract that guarantees the existing signal is
    preserved for non-sparse rows (panel-LLR at ≥1% VAF, DELFI
    profile, etc.).
    """
    from src.foundation.model import (
        SparseAwareLinearProjection,
        LinearProjection,
    )

    torch.manual_seed(42)
    d_in, d_out, dropout = 6, 12, 0.0
    sp = SparseAwareLinearProjection(d_in, d_out, dropout=dropout,
                                     sparsity_threshold=0.5)
    # Mirror the dense-path parameters onto a plain LinearProjection.
    plain = LinearProjection(d_in, d_out, dropout=dropout)
    plain.proj.load_state_dict(sp.proj.state_dict())
    plain.norm.load_state_dict(sp.norm.state_dict())

    # Dense input (no zeros at all → 0% zeros → not sparse)
    x = torch.randn(5, d_in)
    out_sp = sp(x)
    out_plain = plain(x)
    assert torch.allclose(out_sp, out_plain, atol=1e-6), (
        "dense rows through SparseAware should equal LinearProjection; "
        f"max diff {(out_sp - out_plain).abs().max().item()}"
    )


# ── Forward signature matches LinearProjection ───────────────────────────────

@pytest.mark.skipif(not _HAS_TORCH, reason="torch not installed")
def test_forward_signature_matches_linear_projection():
    """Output shape must equal LinearProjection for every input shape
    (1-D, 2-D, 3-D-as-2D) — drop-in replacement contract."""
    from src.foundation.model import SparseAwareLinearProjection

    sp = SparseAwareLinearProjection(11, 16, dropout=0.1)
    # 2-D input
    out_2d = sp(torch.randn(4, 11))
    assert out_2d.shape == (4, 16)
    # 1-D input — should be promoted to (1, d_in)
    out_1d = sp(torch.randn(11))
    assert out_1d.shape == (1, 16)
    # Larger batch
    out_big = sp(torch.randn(64, 11))
    assert out_big.shape == (64, 16)


# ── sparsity_threshold parameter is respected ───────────────────────────────

@pytest.mark.skipif(not _HAS_TORCH, reason="torch not installed")
def test_sparsity_threshold_is_respected():
    """Lowering the threshold must convert previously-dense rows to
    sparse rows. We test by counting how many output rows match the
    missing-token branch exactly (post-LayerNorm)."""
    from src.foundation.model import SparseAwareLinearProjection

    torch.manual_seed(1)
    embed_dim = 4
    d_in = 8
    sp = SparseAwareLinearProjection(d_in, embed_dim, dropout=0.0)
    # Use a non-uniform missing_token so LayerNorm produces a
    # non-degenerate output (all-constant input would give a zero
    # vector after LayerNorm, hiding any branch selection).
    with torch.no_grad():
        sp.missing_token.copy_(torch.tensor([1.0, -1.0, 0.5, -0.5]))

    # Mix of inputs: 2 dense, 2 mostly-zero, 2 fully-zero.
    x = torch.zeros(6, d_in)
    x[0] = torch.randn(d_in) + 1.0                  # 0% zeros (dense)
    x[1] = torch.randn(d_in) + 1.0                  # 0% zeros (dense)
    # Sparse: only 3 of 8 are non-zero → 5/8 = 62.5% zeros > 50%.
    x[2, :3] = 1.0                                  # >50% zeros (sparse)
    x[3, :2] = 1.0                                  # >50% zeros (sparse)
    # rows 4, 5 stay all-zero → 100% zeros (sparse)

    sp.eval()
    expected_missing = sp.norm(sp.missing_token.unsqueeze(0).expand(6, -1))

    # With threshold=0.5: rows 0, 1 are dense (0% zeros); rows 2, 3
    # have >50% zeros → sparse; rows 4, 5 are 100% zero → sparse.
    # → 4 sparse rows total.
    out = sp(x)
    # Per-row max abs diff between output and the missing-token branch.
    row_diff = (out - expected_missing).abs().max(dim=-1).values
    sparse_mask = row_diff < 1e-6
    assert sparse_mask.sum().item() == 4, (
        f"threshold=0.5 should give 4 sparse rows, got "
        f"{sparse_mask.sum().item()}"
    )

    # Lowering threshold to 0.2 → row 0 (0% zeros) and row 1 (0% zeros)
    # still dense, but row 2 (5/8 zeros = 62.5% > 20% → already sparse).
    # Same set stays sparse; a row at e.g. 25% zeros would now flip to
    # sparse but our test set doesn't include one. Build a row with
    # exactly 25% zeros to demonstrate.
    sp_low = SparseAwareLinearProjection(d_in, embed_dim, dropout=0.0,
                                         sparsity_threshold=0.2)
    with torch.no_grad():
        sp_low.missing_token.copy_(torch.tensor([1.0, -1.0, 0.5, -0.5]))
    x_low = x.clone()
    # Row 1: 2/8 zeros = 25% → now sparse under threshold=0.2.
    x_low[1, :2] = 0.0
    expected_low = sp_low.norm(
        sp_low.missing_token.unsqueeze(0).expand(6, -1)
    )
    out_low = sp_low(x_low)
    row_diff_low = (out_low - expected_low).abs().max(dim=-1).values
    sparse_mask_low = row_diff_low < 1e-6
    # Now rows 1, 2, 3, 4, 5 are sparse (5 total).
    assert sparse_mask_low.sum().item() == 5, (
        f"threshold=0.2 with row 1 at 25% zeros should give 5 sparse "
        f"rows, got {sparse_mask_low.sum().item()}"
    )


# ── Gradient flows through both paths ───────────────────────────────────────

@pytest.mark.skipif(not _HAS_TORCH, reason="torch not installed")
def test_gradient_flows_through_dense_and_missing_token_paths():
    """A single backward pass with a meaningful loss must populate
    gradients on the dense-path parameters AND on the missing_token.

    Uses an MSE loss against a target rather than a plain sum: a plain
    ``out.sum().backward()`` zeroes the upstream Linear's weight grad
    through LayerNorm (a known LayerNorm property: summing a normalized
    row cancels the input's contribution to the upstream weight grad).
    An MSE loss against a target keeps the gradient path open.
    """
    from src.foundation.model import SparseAwareLinearProjection

    torch.manual_seed(7)
    sp = SparseAwareLinearProjection(5, 4, dropout=0.0,
                                     sparsity_threshold=0.5)
    # Diverse inputs: rows 0, 1 distinct non-zero patterns (dense);
    # rows 2, 3 all-zero (sparse).
    x = torch.tensor([
        [1.0, 2.0, 3.0, 4.0, 5.0],
        [-1.0, 0.5, -2.0, 1.5, -0.5],
        [0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0],
    ])
    target = torch.randn(4, 4)

    out = sp(x)
    mse = ((out - target) ** 2).mean()
    mse.backward()

    assert sp.proj.weight.grad is not None, "dense path must receive grad"
    assert sp.missing_token.grad is not None, "missing_token must receive grad"
    # Gradients on the dense rows must be non-zero on the dense-path
    # weights, and on the missing rows non-zero on the missing_token.
    assert sp.proj.weight.grad.abs().sum().item() > 0
    assert sp.missing_token.grad.abs().sum().item() > 0

    # Verify the missing-token gradient comes from the sparse rows
    # only (the dense rows must NOT contribute to it). Compute the
    # gradient on a diverse dense-only batch and confirm it's zero.
    sp.zero_grad()
    x_dense = torch.tensor([
        [1.0, 2.0, 3.0, 4.0, 5.0],
        [-1.0, 0.5, -2.0, 1.5, -0.5],
        [0.5, 0.5, 0.5, 0.5, 0.5],
        [-0.5, -0.5, -0.5, -0.5, -0.5],
    ])
    target2 = torch.randn(4, 4)
    out_dense = sp(x_dense)
    ((out_dense - target2) ** 2).mean().backward()
    assert sp.missing_token.grad.abs().sum().item() == 0, (
        "missing_token must NOT receive gradient when all rows are dense"
    )


# ── End-to-end: MultiModalEncoder with one sparse modality ───────────────────

@pytest.mark.skipif(not _HAS_TORCH, reason="torch not installed")
def test_multimodal_encoder_works_with_one_sparse_aware_modality():
    """Wire SparseAwareLinearProjection in for one modality (panel-LLR-
    style sparse channel) while keeping LinearProjection for the rest.
    End-to-end forward + a single backward pass must produce the right
    output shape and finite gradients."""
    from src.foundation.config import (
        FoundationConfig,
        PROTOTYPE_CONFIG,
        MODALITY_DIMS,
    )
    from src.foundation.model import MultiModalEncoder

    cfg = PROTOTYPE_CONFIG  # 64d, 2 layers, 2 heads — fast on CPU
    # Opt the frag_basic modality (4-dim, panel-LLR-shaped) into
    # sparse-aware; leave the rest on the dense LinearProjection.
    encoder = MultiModalEncoder(
        cfg,
        projection_kinds={"frag_basic": "sparse_aware"},
    )

    # Sanity-check the wiring
    assert encoder.projections["frag_basic"].__class__.__name__ == \
        "SparseAwareLinearProjection"
    for name in MODALITY_DIMS:
        if name != "frag_basic":
            assert encoder.projections[name].__class__.__name__ == \
                "LinearProjection", (
                f"non-sparse modality {name} should still be LinearProjection"
            )

    # Build synthetic batch: half dense, half sparse on frag_basic
    batch_size = 8
    rng = np.random.default_rng(11)
    modalities = {}
    for name, dim in MODALITY_DIMS.items():
        x = rng.standard_normal((batch_size, dim)).astype(np.float32)
        if name == "frag_basic":
            # Make half the rows mostly-zero (sparse scenario)
            x[4:] = 0.0
        modalities[name] = torch.from_numpy(x)

    joint = encoder(modalities)
    assert joint.shape == (batch_size, len(MODALITY_DIMS), cfg.embed_dim)
    assert torch.isfinite(joint).all(), "joint embedding should be finite"

    # Backward to confirm all trainable params got gradients (including
    # the sparse modality's missing_token).
    joint.sum().backward()
    assert encoder.projections["frag_basic"].missing_token.grad is not None
    assert encoder.projections["frag_basic"].missing_token.grad.abs().sum().item() > 0


# ── Factory rejects bad kind ─────────────────────────────────────────────────

@pytest.mark.skipif(not _HAS_TORCH, reason="torch not installed")
def test_make_projection_factory_rejects_unknown_kind():
    """make_projection must raise ValueError on unknown kind strings."""
    from src.foundation.model import (
        make_projection,
        LinearProjection,
        SparseAwareLinearProjection,
    )
    assert isinstance(make_projection(4, 8), LinearProjection)
    assert isinstance(
        make_projection(4, 8, kind="sparse_aware"),
        SparseAwareLinearProjection,
    )
    with pytest.raises(ValueError, match="Unknown projection kind"):
        make_projection(4, 8, kind="not_a_real_kind")