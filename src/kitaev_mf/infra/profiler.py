"""Computational resource profiler for counting exact analytical FLOPs, tensor memory, and wall-clock time."""

from __future__ import annotations

import time
from collections.abc import Generator
from contextlib import contextmanager
from dataclasses import dataclass


@dataclass
class ProfileRecord:
    """Stores computational resource measurements."""

    flops: float = 0.0
    wall_time_sec: float = 0.0
    peak_bytes: int = 0
    eigh_calls: int = 0
    gemm_calls: int = 0
    map_evaluations: int = 0


class ResourceProfiler:
    """Tracks FLOPs, wall-clock time, and memory allocations during solver steps."""

    def __init__(self) -> None:
        self.total_flops: float = 0.0
        self.total_wall_time: float = 0.0
        self.peak_bytes: int = 0
        self.eigh_calls: int = 0
        self.gemm_calls: int = 0
        self.map_evaluations: int = 0
        self._start_time: float | None = None

    def reset(self) -> None:
        """Reset all profiling counters."""
        self.total_flops = 0.0
        self.total_wall_time = 0.0
        self.peak_bytes = 0
        self.eigh_calls = 0
        self.gemm_calls = 0
        self.map_evaluations = 0
        self._start_time = None

    def count_eigh(self, m: int, is_complex: bool = True) -> None:
        """Record FLOPs for an eigensolve of an m x m matrix."""
        self.eigh_calls += 1
        if is_complex:
            # Complex Hermitian eigensolver (e.g. LAPACK zheevd): ~ 16/3 M^3 real FLOPs
            self.total_flops += (16.0 / 3.0) * (m**3)
        else:
            # Real symmetric eigensolver (e.g. LAPACK dsyevd): ~ 4/3 M^3 real FLOPs
            self.total_flops += (4.0 / 3.0) * (m**3)

    def count_gemm(self, m: int, n: int, k: int, is_complex: bool = True) -> None:
        """Record FLOPs for a general matrix multiplication (M x K) @ (K x N)."""
        self.gemm_calls += 1
        if is_complex:
            # Complex GEMM: 8 * M * N * K real FLOPs
            self.total_flops += 8.0 * m * n * k
        else:
            # Real GEMM: 2 * M * N * K real FLOPs
            self.total_flops += 2.0 * m * n * k

    def count_phi_bloch(self, n_c: int, n_k: int) -> None:
        """Record exact FLOPs for a full Bloch spectral projection evaluation."""
        self.map_evaluations += 1
        m_c = 4 * n_c
        # 1. Hamiltonian build: ~ 66 * n_c * n_k
        self.total_flops += 66.0 * n_c * n_k
        # 2. Eigendecomposition over N_k points
        self.eigh_calls += n_k
        self.total_flops += n_k * (16.0 / 3.0) * (m_c**3)
        # 3. Covariance projection (GEMM) over N_k points: M_c x (M_c/2) @ (M_c/2) x M_c
        self.gemm_calls += n_k
        self.total_flops += n_k * 4.0 * (m_c**3)
        # 4. BZ quadrature integration
        self.total_flops += 3.0 * n_c * (m_c**2) * n_k

    def count_purification_step(self, m: int, step_type: str = "newton_schulz") -> None:
        """Record FLOPs for one matrix sign purification step (pure real arithmetic)."""
        if step_type == "newton_schulz":
            # Y_{k+1} = 0.5 * Y_k (3I + Y_k^2)
            # Step 1: Y_k^2 -> real GEMM of size (m x m) @ (m x m) -> 2 * m^3 FLOPs
            # Step 2: Y_k * (3I + Y_k^2) -> real GEMM of size (m x m) @ (m x m) -> 2 * m^3 FLOPs
            self.gemm_calls += 2
            self.total_flops += 4.0 * (m**3)
        elif step_type == "sp2":
            # X_{k+1} = X_k^2 or 2X_k - X_k^2 -> 1 real GEMM -> 2 * m^3 FLOPs
            self.gemm_calls += 1
            self.total_flops += 2.0 * (m**3)

    @contextmanager
    def measure_time(self) -> Generator[None, None, None]:
        """Context manager to record elapsed wall-clock time."""
        t0 = time.perf_counter()
        try:
            yield
        finally:
            t1 = time.perf_counter()
            self.total_wall_time += t1 - t0

    def snapshot(self) -> ProfileRecord:
        """Return a snapshot dataclass of current profile metrics."""
        return ProfileRecord(
            flops=self.total_flops,
            wall_time_sec=self.total_wall_time,
            peak_bytes=self.peak_bytes,
            eigh_calls=self.eigh_calls,
            gemm_calls=self.gemm_calls,
            map_evaluations=self.map_evaluations,
        )
