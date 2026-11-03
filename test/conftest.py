"""Shared pytest fixtures."""

import pytest

from kitaev_mf.infra.linalg_backend import TorchBackend


@pytest.fixture
def cpu_backend() -> TorchBackend:
    """Fixture providing a CPU-bound TorchBackend for unit tests."""
    return TorchBackend(device="cpu")
