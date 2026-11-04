SHELL := /bin/bash
VENV := .venv
PYTHON := $(if $(wildcard $(VENV)/bin/python),$(VENV)/bin/python,python3)
PIP := $(VENV)/bin/pip
PIP_COMPILE := $(VENV)/bin/pip-compile
PIP_SYNC := $(VENV)/bin/pip-sync
PYTEST := $(VENV)/bin/pytest
RUFF := $(VENV)/bin/ruff
MYPY := $(VENV)/bin/mypy
SNAPSHOT_FILE := REPO_SNAPSHOT.md

.PHONY: help venv lock sync install test lint format verify snapshot clean

help:
	@echo "kitaev-mf development targets:"
	@echo "  make venv     - Create .venv and install pip-tools"
	@echo "  make lock     - Compile requirements.in -> requirements.txt via pip-compile"
	@echo "  make sync     - Synchronize venv to exact locked requirements"
	@echo "  make install  - Editable install of kitaev-mf in venv"
	@echo "  make test     - Run pytest suite with coverage"
	@echo "  make lint     - Run Ruff lint/format checks and Mypy"
	@echo "  make format   - Automatically fix lint and reformat all files via Ruff"
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

format:
	@$(RUFF) check --fix src test scripts
	@$(RUFF) format src test scripts

lint:
	@echo "[lint] Ruff rules across src/, test/, and scripts/..."
	@$(RUFF) check src test scripts
	@echo "[format] Check production code and tests touched in the active batch..."
	@$(RUFF) format --check src scripts test/unit/test_energetics.py test/unit/test_k006_quadrature.py
	@echo "[types] Mypy across src/, test/, and scripts/..."
	@$(MYPY) src test scripts

verify:
	@set -eu; \
	printf '\n[1/4] Unit tests with coverage\n'; \
	$(MAKE) --no-print-directory test; \
	printf '\n[2/4] Ruff lint, formatting, and Mypy\n'; \
	$(MAKE) --no-print-directory lint; \
	printf '\n[3/4] Git whitespace checks (committed branch diff and local changes)\n'; \
	git diff --check main...HEAD; \
	git diff --cached --check; \
	git diff --check; \
	printf '\n[4/4] Working-tree status (informational; does not require a clean tree)\n'; \
	git status --short --branch; \
	printf '\nVerification completed: all preceding checks passed.\n'

snapshot:
	@$(PYTHON) scripts/make_snapshot.py

clean:
	@rm -rf build/ dist/ *.egg-info .pytest_cache/ .coverage .mypy_cache/ .ruff_cache/ $(SNAPSHOT_FILE) repo_snapshot.txt
	@find . -type d -name "__pycache__" -exec rm -rf {} +
