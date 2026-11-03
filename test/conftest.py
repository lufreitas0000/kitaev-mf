"""Shared pytest fixtures."""
import pytest
from kitaev_mf.infra.linalg_backend import TorchBackend

@pytest.fixture
def cpu_backend() -> TorchBackend:
    return TorchBackend(device="cpu")
