"""Locally scoped 64-bit JAX precision helpers for calibration entry points.

``mundey`` never enables ``jax_enable_x64`` globally on import. Numerical entry
points run inside :func:`run_in` (a locally scoped ``jax.enable_x64`` context)
after casting floating-point inputs via :func:`cast_tree`, restoring the caller's
prior JAX precision state on exit.
"""

from __future__ import annotations

import contextlib
import functools
from collections.abc import Callable, Iterator
from typing import Any, TypeVar

import jax
import jax.numpy as jnp
import numpy as np

F = TypeVar("F", bound=Callable[..., Any])

_DTYPES: dict[str, tuple[Any, Any]] = {
    "float32": (jnp.float32, jnp.complex64),
    "float64": (jnp.float64, jnp.complex128),
}


def _check_dtype(dtype: str) -> None:
    if dtype not in _DTYPES:
        raise ValueError(f"dtype must be 'float32' or 'float64', not {dtype!r}.")


@contextlib.contextmanager
def run_in(dtype: str = "float64") -> Iterator[None]:
    """Execute the enclosed block with JAX 64-bit mode enabled or disabled.

    Parameters
    ----------
    dtype : {"float32", "float64"}
        Target floating-point precision. The previous JAX ``enable_x64`` state
        is unconditionally restored on exit.
    """
    _check_dtype(dtype)
    with jax.enable_x64(dtype == "float64"):
        yield


def cast_tree(tree: Any, dtype: str = "float64") -> Any:
    """Cast floating and complex array leaves in a PyTree to ``dtype``.

    Integer, boolean, and non-array leaves pass through unchanged. Must be
    called inside ``run_in("float64")`` when ``dtype == "float64"`` so JAX
    retains 64-bit device buffers.
    """
    _check_dtype(dtype)
    real_dtype, complex_dtype = _DTYPES[dtype]

    def _cast_leaf(leaf: Any) -> Any:
        if not isinstance(leaf, (jax.Array, np.ndarray, np.generic)):
            return leaf
        if jnp.issubdtype(leaf.dtype, jnp.complexfloating):
            return jnp.asarray(leaf, dtype=complex_dtype)
        if jnp.issubdtype(leaf.dtype, jnp.floating):
            return jnp.asarray(leaf, dtype=real_dtype)
        return leaf

    return jax.tree_util.tree_map(_cast_leaf, tree)


def in_float64(fn: F) -> F:
    """Run ``fn`` inside locally scoped JAX ``float64`` with cast array inputs.

    Non-array arguments and iterators pass through unchanged, and the caller's
    prior ``jax_enable_x64`` configuration is restored when ``fn`` returns or
    raises.
    """

    def _cast_arg(arg: Any) -> Any:
        return arg if isinstance(arg, Iterator) else cast_tree(arg, "float64")

    @functools.wraps(fn)
    def wrapped(*args: Any, **kwargs: Any) -> Any:
        with run_in("float64"):
            cast_args = [_cast_arg(a) for a in args]
            cast_kwargs = {k: _cast_arg(v) for k, v in kwargs.items()}
            return fn(*cast_args, **cast_kwargs)

    return wrapped  # type: ignore[return-value]
