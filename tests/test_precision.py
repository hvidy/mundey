"""Tests for locally scoped 64-bit JAX precision helpers in mundey._precision."""

from __future__ import annotations

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from mundey._precision import cast_tree, in_float64, run_in


def test_in_float64_promotes_inputs_and_restores_default_float32_state() -> None:
    """@in_float64 must run in float64 locally without leaking global x64 state."""
    assert jax.config.jax_enable_x64 is False

    @in_float64
    def accumulate(tree: dict[str, object], count: int) -> tuple[jax.Array, jax.Array, int, bool]:
        inside_flag = bool(jax.config.jax_enable_x64)
        arr = tree["x"]
        cplx = tree["z"]
        assert isinstance(arr, jax.Array)
        assert isinstance(cplx, jax.Array)
        assert arr.dtype == jnp.float64
        assert cplx.dtype == jnp.complex128
        # Verify precision-sensitive cancellation (1.0 + 1e-12 - 1.0 != 0 in float64)
        val = (jnp.ones(1) + 1e-12) - 1.0
        return arr + val, cplx, count, inside_flag

    payload = {
        "x": np.array([1.0, 2.0], dtype=np.float32),
        "z": jnp.asarray(np.array([1.0 + 2.0j], dtype=np.complex64)),
        "idx": np.array([3, 4], dtype=np.int32),
        "label": "output-A",
    }
    out_real, out_cplx, out_count, was_x64 = accumulate(payload, count=7)

    assert was_x64 is True
    assert out_real.dtype == jnp.float64
    assert out_cplx.dtype == jnp.complex128
    assert out_count == 7
    assert np.isclose(float(np.asarray(out_real)[0] - 1.0), 1e-12, rtol=1e-4)
    assert jax.config.jax_enable_x64 is False


def test_in_float64_restores_prior_state_on_exception_and_nested_contexts() -> None:
    """Prior JAX x64 state (False or True) must be restored even if an exception occurs."""

    @in_float64
    def boom(x: np.ndarray) -> jax.Array:
        assert x.dtype == jnp.float64
        raise RuntimeError("boom inside x64")

    assert jax.config.jax_enable_x64 is False
    with pytest.raises(RuntimeError, match="boom inside x64"):
        boom(np.ones(3, dtype=np.float32))
    assert jax.config.jax_enable_x64 is False

    with run_in("float64"):
        assert jax.config.jax_enable_x64 is True
        with pytest.raises(RuntimeError, match="boom inside x64"):
            boom(np.ones(3, dtype=np.float32))
        assert jax.config.jax_enable_x64 is True

        with run_in("float32"):
            assert jax.config.jax_enable_x64 is False
            cast_32 = cast_tree({"a": np.ones(2, dtype=np.float64)}, "float32")
            assert cast_32["a"].dtype == jnp.float32
        assert jax.config.jax_enable_x64 is True

    assert jax.config.jax_enable_x64 is False


def test_iterators_pass_through_unchanged_and_invalid_dtype_rejected() -> None:
    """Iterator arguments must not be eagerly consumed by @in_float64."""

    @in_float64
    def take_first(it) -> float:
        return float(next(it))

    gen = iter([1.25, 2.5, 3.75])
    assert take_first(gen) == 1.25
    assert next(gen) == 2.5

    with pytest.raises(ValueError, match="dtype must be"):
        with run_in("float16"):
            pass
