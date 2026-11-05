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

.PHONY: help venv lock sync install test test-verbose benchmark k006-quadrature lint format verify snapshot clean

help:
	@echo "kitaev-mf targets:"
	@echo "  make verify        - Run concise complete gate (tests, static checks, whitespace)"
	@echo "  make test          - Run pytest suite concisely with coverage"
	@echo "  make test-verbose  - Run pytest suite with full per-test verbosity"
	@echo "  make benchmark     - Generate and export reproducible .dat benchmark tables"
	@echo "  make format        - Automatically fix lint, format code, and sanitize whitespace"
	@echo "  make lint          - Run Ruff check, format verification, and strict Mypy"
	@echo "  make snapshot      - Generate REPO_SNAPSHOT.md for agent ingestion"

venv:
	@if [ ! -d "$(VENV)" ]; then \
		echo "Creating virtual environment..."; \
		python3 -m venv $(VENV); \
	fi
	@$(PIP) install --upgrade pip setuptools wheel pip-tools

lock: venv
	@$(PIP_COMPILE) --resolver=backtracking -o requirements.txt requirements.in
	@$(PIP_COMPILE) --resolver=backtracking -o requirements-dev.txt requirements-dev.in

sync: venv
	@$(PIP_SYNC) requirements.txt requirements-dev.txt
	@$(PIP) install -e . --no-deps

install: venv lock sync

test:
	@$(PYTEST) -q --tb=short --cov=kitaev_mf --cov-report=term-missing:skip-covered test/

test-verbose:
	@$(PYTEST) -v --tb=short --cov=kitaev_mf --cov-report=term-missing test/

benchmark:
	@echo "[benchmark] Generating reproducible .dat benchmark files in data/processed/benchmarks/..."
	@$(PYTHON) scripts/generate_benchmarks.py --benchmark all

k006-quadrature:
	@$(PYTHON) scripts/k006_quadrature.py

format:
	@$(RUFF) check --fix src test scripts
	@$(RUFF) format src test scripts
	@$(PYTHON) scripts/fix_whitespace.py

lint:
	@$(RUFF) check src test scripts
	@$(RUFF) format --check src test scripts
	@$(MYPY) src test scripts

verify:
	@set -eu; \
	printf '[1/4] Running test suite... '; \
	$(PYTEST) -q --tb=short --cov=kitaev_mf --cov-report=term:skip-covered test/ > /dev/null 2>&1 \
		&& echo 'PASS' \
		|| { echo 'FAILED'; $(PYTEST) --tb=short test/; exit 1; }; \
	printf '[2/4] Static analysis (Ruff & Mypy)... '; \
	if $(RUFF) check src test scripts > /dev/null 2>&1 \
		&& $(RUFF) format --check src test scripts > /dev/null 2>&1 \
		&& $(MYPY) src test scripts > /dev/null 2>&1; then \
		echo 'PASS'; \
	else \
		echo 'FAILED'; \
		$(MAKE) lint; \
		exit 1; \
	fi; \
	printf '[3/4] Whitespace & Git hygiene... '; \
	if git diff --check -- . ':(exclude)doc/reference' > /dev/null 2>&1; then \
		echo 'PASS'; \
	else \
		echo 'FAILED'; \
		git diff --check -- . ':(exclude)doc/reference'; \
		exit 1; \
	fi; \
	printf '[4/4] Coverage summary:\n'; \
	$(PYTEST) -q --cov=kitaev_mf --cov-report=term:skip-covered test/ | grep -E "(TOTAL|TOTAL.*%|===)" || true; \
	printf '\nVerification completed successfully.\n'

snapshot:
	@$(PYTHON) scripts/make_snapshot.py

clean:
	@rm -rf build/ dist/ *.egg-info .pytest_cache/ .coverage .mypy_cache/ .ruff_cache/ $(SNAPSHOT_FILE) repo_snapshot.txt
	@find . -type d -name "__pycache__" -exec rm -rf {} +
