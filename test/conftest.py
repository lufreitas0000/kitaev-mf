"""Shared pytest fixtures."""

import pytest

from src.infra.linalg_backend import TorchBackend


@pytest.fixture
def cpu_backend() -> TorchBackend:
    return TorchBackend(device="cpu")
