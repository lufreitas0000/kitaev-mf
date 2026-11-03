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

.PHONY: help venv lock sync install test lint format snapshot clean

help:
	@echo "kitaev-mf development targets:"
	@echo "  make venv     - Create .venv and install pip-tools"
	@echo "  make lock     - Compile requirements.in -> requirements.txt via pip-compile"
	@echo "  make sync     - Synchronize venv to exact locked requirements"
	@echo "  make install  - Editable install of kitaev-mf in venv"
	@echo "  make test     - Run pytest suite with coverage"
	@echo "  make lint     - Run ruff format/lint checks and mypy"
	@echo "  make format   - Automatically fix lint and reformat all files via ruff"
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
	@$(RUFF) check src test scripts
	@$(RUFF) format --check src test scripts
	@$(MYPY) src test scripts

snapshot:
	@$(PYTHON) scripts/make_snapshot.py

clean:
	@rm -rf build/ dist/ *.egg-info .pytest_cache/ .coverage .mypy_cache/ .ruff_cache/ $(SNAPSHOT_FILE) repo_snapshot.txt
	@find . -type d -name "__pycache__" -exec rm -rf {} +
