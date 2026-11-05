SHELL := /bin/bash
VENV := .venv
PYTHON := $(if $(wildcard $(VENV)/bin/python3),$(VENV)/bin/python3,$(if $(wildcard $(VENV)/bin/python),$(VENV)/bin/python,python3))
PIP := $(VENV)/bin/pip
PIP_COMPILE := $(VENV)/bin/pip-compile
PIP_SYNC := $(VENV)/bin/pip-sync
PYTEST := $(VENV)/bin/pytest
RUFF := $(VENV)/bin/ruff
MYPY := $(VENV)/bin/mypy
SNAPSHOT_FILE := REPO_SNAPSHOT.md

.PHONY: help venv lock sync install test k006-quadrature lint format verify snapshot clean

help:
	@echo "kitaev-mf development targets:"
	@echo "  make venv     - Create .venv and install pip-tools"
	@echo "  make lock     - Compile requirements.in -> requirements.txt via pip-compile"
	@echo "  make sync     - Synchronize venv to exact locked requirements"
	@echo "  make install  - Editable install of kitaev-mf in venv"
	@echo "  make test     - Run pytest suite with coverage"
	@echo "  make k006-quadrature - Run the high-resolution K-006 convergence diagnostic"
	@echo "  make lint     - Run Ruff lint/format checks and Mypy"
	@echo "  make format   - Automatically fix lint, format code, and sanitize whitespace"
	@echo "  make verify   - Run tests, static checks, Git whitespace checks, and status"
	@echo "  make snapshot - Generate uncommitted $(SNAPSHOT_FILE) for AI context ingestion"
	@echo "  make clean    - Remove build, test caches, and temporary snapshots"

venv:
	@if [ ! -d "$(VENV)" ]; then \
		echo "Creating virtual environment..."; \
		python3 -m venv $(VENV); \
	fi
	@$(PIP) install --upgrade pip setuptools wheel pip-tools

lock: venv
	@echo "Compiling requirements.txt..."
	@$(PIP_COMPILE) --resolver=backtracking -o requirements.txt requirements.in
	@echo "Compiling requirements-dev.txt..."
	@$(PIP_COMPILE) --resolver=backtracking -o requirements-dev.txt requirements-dev.in

sync: venv
	@$(PIP_SYNC) requirements.txt requirements-dev.txt
	@$(PIP) install -e . --no-deps

install: venv lock sync

test:
	@$(PYTEST) --cov=kitaev_mf --cov-report=term-missing test/

k006-quadrature:
	@echo "[K-006] Recompute the midpoint sequence and empirical convergence diagnostics..."
	@$(PYTHON) scripts/k006_quadrature.py

format:
	@echo "[format] Ruff check & fix..."
	@$(RUFF) check --fix src test scripts
	@echo "[format] Ruff format..."
	@$(RUFF) format src test scripts
	@echo "[format] Sanitizing file whitespace and EOF lines..."
	@$(PYTHON) scripts/fix_whitespace.py

lint:
	@echo "[lint] Ruff rules across src/, test/, and scripts/..."
	@$(RUFF) check src test scripts
	@echo "[format] Check formatting across the complete Python source, test, and script tree..."
	@$(RUFF) format --check src test scripts
	@echo "[types] Mypy across src/, test/, and scripts/..."
	@$(MYPY) src test scripts

verify:
	@set -eu; \
	printf '\n[1/5] Unit tests with coverage\n'; \
	$(MAKE) --no-print-directory test; \
	printf '\n[2/5] K-006 high-resolution numerical diagnostic\n'; \
	$(MAKE) --no-print-directory k006-quadrature; \
	printf '\n[3/5] Ruff lint, formatting, and Mypy\n'; \
	$(MAKE) --no-print-directory lint; \
	printf '\n[4/5] Git whitespace checks (excluding read-only doc/reference/)\n'; \
	git diff --check main...HEAD -- . ':(exclude)doc/reference'; \
	git diff --cached --check -- . ':(exclude)doc/reference'; \
	git diff --check -- . ':(exclude)doc/reference'; \
	printf '\n[5/5] Working-tree status (informational; does not require a clean tree)\n'; \
	git status --short --branch; \
	printf '\nVerification completed: all preceding checks passed.\n'

snapshot:
	@$(PYTHON) scripts/make_snapshot.py

clean:
	@rm -rf build/ dist/ *.egg-info .pytest_cache/ .coverage .mypy_cache/ .ruff_cache/ $(SNAPSHOT_FILE) repo_snapshot.txt
	@find . -type d -name "__pycache__" -exec rm -rf {} +
