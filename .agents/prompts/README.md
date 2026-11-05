# Portable Multi-Agent Preference Kit

This directory contains your standardized, cross-platform system prompts and role headers so that every AI assistant—**Gemini Web, ChatGPT Web, Claude Web, Qwen Web, Qwen Local (Ollama), DeepSeek Web, Google Antigravity, and Codex Desktop**—follows the exact same scientific, coding, editorial, and LaTeX conventions.

---

## 1. Architecture of the Kit (3 Layers)

| Layer | File(s) | Purpose | Target Token Budget |
| :--- | :--- | :--- | :---: |
| **Layer 1: Universal Profile** | [`00_universal_profile.md`](./00_universal_profile.md) | Permanent system instructions defining your mathematical rigor, Three-Front R&D layout, LaTeX design style, software architecture, and version-control rules. | ~950 tokens |
| **Layer 2: Role Task Headers** | [`roles/role_rnd_coder.md`](./roles/role_rnd_coder.md)<br/>[`roles/role_labbook_synthesizer.md`](./roles/role_labbook_synthesizer.md)<br/>[`roles/role_manuscript_editor.md`](./roles/role_manuscript_editor.md)<br/>[`roles/role_adversarial_referee.md`](./roles/role_adversarial_referee.md) | Task-specific operational contracts prepended at the start of a chat session depending on which R&D front or audit gate you are running. | ~300 tokens each |
| **Layer 3: Repo Context Bundle** | [`../templates/make_snapshot.py`](../templates/make_snapshot.py) | Deterministic script that packages `PROJECT_STATE.yaml`, `LOG_DISCOVERY_AND_PIVOTS.md`, active `doc/rnd/` notes, and code into a single `REPO_SNAPSHOT.md` file. | Configurable (`compact`, `active`, `full`) |

---

## 2. How to Use Across Your AI Fleet

### A. Web Chat UIs (ChatGPT Web, Claude Web, Gemini Web, Qwen Web, DeepSeek Web)
1. **One-Time Setup:** Paste the contents of [`00_universal_profile.md`](./00_universal_profile.md) into:
   - **ChatGPT Web:** *Project Instructions* (or *Custom Instructions*)
   - **Claude Web:** *Project Instructions* (or a custom *Style*)
   - **Gemini Web:** *Gem Instructions* (or saved info)
   - **Qwen Web / DeepSeek Web:** Pinned system prompt or initial turn
2. **Per-Session Dispatch:** Generate a snapshot with your desired role header embedded:
   ```bash
   python3 scripts/make_snapshot.py --mode active --role rnd
   ```
   Upload or paste the generated `REPO_SNAPSHOT.md` and state your concrete task (`T1`, `T2`, etc.).

### B. Local LLMs (Qwen on Ollama)
1. Create an Ollama `Modelfile` using [`00_universal_profile.md`](./00_universal_profile.md) as the `SYSTEM` block.
2. Generate a compact snapshot tailored for local context windows:
   ```bash
   python3 scripts/make_snapshot.py --mode compact --role rnd --output /tmp/qwen_context.md
   ```

### C. Desktop & CLI Agents (Google Antigravity, Codex Desktop)
- Desktop agents automatically read `AGENTS.md`, `.agents/rules/`, and `PROJECT_STATE.yaml` from the repository root. You can also reference `@prompts/roles/role_adversarial_referee.md` when invoking an independent review subagent.
