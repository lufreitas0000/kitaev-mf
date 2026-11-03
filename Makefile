SHELL := /bin/bash
VENV := .venv
PYTHON := $(VENV)/bin/python
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
	@$(RUFF) check --fix src test
	@$(RUFF) format src test

lint:
	@$(RUFF) check src test
	@$(RUFF) format --check src test
	@$(MYPY) src test

snapshot:
	@echo "Generating $(SNAPSHOT_FILE)..."
	@echo "# Repository Context Snapshot: kitaev-mf" > $(SNAPSHOT_FILE)
	@echo "" >> $(SNAPSHOT_FILE) @echo "- **Generated:** $$(date -u '+%Y-%m-%d %H:%M:%S UTC')" >> $(SNAPSHOT_FILE) @echo "- **Git Branch:** $$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo 'unknown')" >> $(SNAPSHOT_FILE) @echo "- **Git Commit:** $$(git rev-parse --short HEAD 2>/dev/null \vert{}\vert{} echo 'unknown')" >> $(SNAPSHOT_FILE)
	@echo "" >> $(SNAPSHOT_FILE)
	@echo "## 1. Directory Tree" >> $(SNAPSHOT_FILE)
	@echo '```text' >> $(SNAPSHOT_FILE)
	@find . -maxdepth 4 \
		-not -path '*/.*' \
		-not -path './.venv*' \
		-not -path './doc/reference*' \
		-not -path '*__pycache__*' \
		-not -path './*.egg-info*' \
		| sort >> $(SNAPSHOT_FILE)
	@echo '```' >> $(SNAPSHOT_FILE)
	@echo "" >> $(SNAPSHOT_FILE)
	@echo "## 2. Governance, State & Documentation" >> $(SNAPSHOT_FILE) @for f in AGENTS.md PROJECT_STATE.yaml pyproject.toml Makefile .agents/rules/*.md doc/adr/*.md doc/rnd/*.md; do \ 	if [ -f "$$f" ]; then \ 		echo "### File: \`$$f\`" >> $(SNAPSHOT_FILE); \
			echo '```markdown' >> $(SNAPSHOT_FILE); \
			cat "$$f" >> $(SNAPSHOT_FILE); \
			echo "" >> $(SNAPSHOT_FILE); \
			echo '```' >> $(SNAPSHOT_FILE); \
			echo "" >> $(SNAPSHOT_FILE); \
		fi; \
	done
	@echo "## 3. Source Code & Tests" >> $(SNAPSHOT_FILE) @for f in $$(find src test -name "*.py" \vert{} sort); do \ 	echo "### File: \`$$f\`" >> $(SNAPSHOT_FILE); \
		echo '```python' >> $(SNAPSHOT_FILE); \ 	cat "$$f" >> $(SNAPSHOT_FILE); \
		echo "" >> $(SNAPSHOT_FILE); \
		echo '```' >> $(SNAPSHOT_FILE); \
		echo "" >> $(SNAPSHOT_FILE); \
	done
	@echo "Snapshot written to $(SNAPSHOT_FILE) ($$(wc -l < $(SNAPSHOT_FILE)) lines)."

clean:
	@rm -rf build/ dist/ *.egg-info .pytest_cache/ .coverage .mypy_cache/ .ruff_cache/ $(SNAPSHOT_FILE) repo_snapshot.txt
	@find . -type d -name "__pycache__" -exec rm -rf {} +
