SHELL := /bin/bash
MAIN := main
TEX_ENGINE := pdflatex

# Terminal ANSI Palette
CLR_RED := \033[1;31m
CLR_GRN := \033[1;32m
CLR_YEL := \033[1;33m
CLR_BLU := \033[1;34m
CLR_CYN := \033[1;36m
CLR_RST := \033[0m

.PHONY: all pdf clean

all: pdf

pdf:
	@printf "$(CLR_BLU)==> [1/2] Pass 1: Typesetting LaTeX source and benchmark tables...$(CLR_RST)\n"
	@START_TIME=$$(date +%s); \
	if ! $(TEX_ENGINE) -halt-on-error -file-line-error -interaction=nonstopmode $(MAIN).tex > build.log 2>&1; then \
		printf "$(CLR_RED)==> [FAIL] LaTeX compilation failed at Pass 1! Error details:$(CLR_RST)\n"; \
		awk '/:[0-9]+: |^!/ {print "$(CLR_RED)" $$0 "$(CLR_RST)"; next} \
		     /^(l\.[0-9]+|<argument>)/ {print "$(CLR_YEL)" $$0 "$(CLR_RST)"; next} \
		     {print}' build.log | grep -E -A 3 -B 1 "(\.tex:[0-9]+:|!)" | head -n 35; \
		exit 1; \
	fi; \
	printf "$(CLR_BLU)==> [2/2] Pass 2: Resolving cross-references and Table of Contents...$(CLR_RST)\n"; \
	if ! $(TEX_ENGINE) -halt-on-error -file-line-error -interaction=nonstopmode $(MAIN).tex >> build.log 2>&1; then \
		printf "$(CLR_RED)==> [FAIL] LaTeX compilation failed at Pass 2! Error details:$(CLR_RST)\n"; \
		awk '/:[0-9]+: |^!/ {print "$(CLR_RED)" $$0 "$(CLR_RST)"; next} \
		     /^(l\.[0-9]+|<argument>)/ {print "$(CLR_YEL)" $$0 "$(CLR_RST)"; next} \
		     {print}' build.log | grep -E -A 3 -B 1 "(\.tex:[0-9]+:|!)" | head -n 35; \
		exit 1; \
	fi; \
	END_TIME=$$(date +%s); \
	PAGES=$$(pdfinfo $(MAIN).pdf 2>/dev/null | awk '/Pages:/ {print $$2}'); \
	PAGES=$${PAGES:-$$(grep -c '\[[0-9]\+\]' build.log || echo '?')}; \
	ELAPSED=$$((END_TIME - START_TIME)); \
	printf "$(CLR_GRN)==> [SUCCESS] Generated doc/labbook/$(MAIN).pdf ($$PAGES pages, $${ELAPSED}s)$(CLR_RST)\n"

clean:
	@printf "$(CLR_YEL)==> Cleaning labbook build artifacts...$(CLR_RST)\n"
	@rm -f *.aux *.log *.out *.toc *.fls *.fdb_latexmk *.synctex.gz build.log
	@rm -f chapters/*.aux
