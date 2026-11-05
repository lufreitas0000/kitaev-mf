#!/usr/bin/env python3
"""Deterministic CLI tool for managing tasks in PROJECT_STATE.yaml.

Usage:
    python3 scripts/manage_project_state.py transition-task --done T1 --next T2
    python3 scripts/manage_project_state.py set-task-status --id T2 --status completed
    python3 scripts/manage_project_state.py add-task --id T5 --title "..." --files "f1,f2" --contracts "c1,c2" --validation "pytest ..."
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

import yaml

STATE_FILE = Path("PROJECT_STATE.yaml")
ALLOWED_STATUSES = {
    "planned",
    "in_progress",
    "blocked",
    "awaiting_local_validation",
    "completed",
}


def load_state() -> dict[str, Any]:
    if not STATE_FILE.is_file():
        print(f"Error: {STATE_FILE} not found.", file=sys.stderr)
        sys.exit(1)
    with STATE_FILE.open("r", encoding="utf-8") as f:
        data = yaml.safe_load(f)
    if not isinstance(data, dict):
        print(f"Error: {STATE_FILE} is not a valid YAML dictionary.", file=sys.stderr)
        sys.exit(1)
    return data


def save_state(data: dict[str, Any]) -> None:
    # Validate minimal schema consistency
    active_m = data.get("active_milestone", {})
    if data.get("current_milestone") != active_m.get("id"):
        print(
            f"Error: current_milestone ({data.get('current_milestone')}) != "
            f"active_milestone.id ({active_m.get('id')})",
            file=sys.stderr,
        )
        sys.exit(1)

    # Dump formatted YAML
    content = yaml.dump(data, sort_keys=False, allow_unicode=True)

    # Ensure single trailing newline
    content = content.rstrip() + "\n"
    STATE_FILE.write_text(content, encoding="utf-8")
    print(f"Successfully updated {STATE_FILE}")


def find_task(active_tasks: list[dict[str, Any]], task_id: str) -> dict[str, Any] | None:
    for task in active_tasks:
        if task.get("id") == task_id:
            return task
    return None


def cmd_set_task_status(args: argparse.Namespace) -> None:
    data = load_state()
    active_tasks: list[dict[str, Any]] = data.get("active_milestone", {}).get("active_tasks", [])
    task = find_task(active_tasks, args.id)
    if task is None:
        print(f"Error: Task {args.id} not found in active_milestone.active_tasks.", file=sys.stderr)
        sys.exit(1)

    if args.status not in ALLOWED_STATUSES:
        print(
            f"Error: Invalid status '{args.status}'. Allowed: {ALLOWED_STATUSES}", file=sys.stderr
        )
        sys.exit(1)

    old_status = task.get("status")
    task["status"] = args.status
    print(f"Task {args.id}: status changed from '{old_status}' to '{args.status}'")
    save_state(data)


def cmd_transition_task(args: argparse.Namespace) -> None:
    data = load_state()
    active_tasks: list[dict[str, Any]] = data.get("active_milestone", {}).get("active_tasks", [])

    task_done = find_task(active_tasks, args.done)
    if task_done is None:
        print(f"Error: Task {args.done} not found.", file=sys.stderr)
        sys.exit(1)

    task_next = find_task(active_tasks, args.next)
    if task_next is None:
        print(f"Error: Task {args.next} not found.", file=sys.stderr)
        sys.exit(1)

    task_done["status"] = "completed"
    task_next["status"] = "in_progress"
    print(f"Task {args.done} marked as 'completed'")
    print(f"Task {args.next} marked as 'in_progress'")
    save_state(data)


def cmd_add_task(args: argparse.Namespace) -> None:
    data = load_state()
    active_tasks: list[dict[str, Any]] = data.get("active_milestone", {}).get("active_tasks", [])

    if find_task(active_tasks, args.id) is not None:
        print(f"Error: Task {args.id} already exists.", file=sys.stderr)
        sys.exit(1)

    files = [f.strip() for f in args.files.split(",") if f.strip()]
    contracts = [c.strip() for c in args.contracts.split(",") if c.strip()]

    new_task: dict[str, Any] = {
        "id": args.id,
        "title": args.title,
        "status": args.status,
        "target_files": files,
        "contracts": contracts,
        "validation": args.validation,
    }

    active_tasks.append(new_task)
    print(f"Added task {args.id}: {args.title}")
    save_state(data)


def main() -> None:
    parser = argparse.ArgumentParser(description="Deterministic manager for PROJECT_STATE.yaml")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # set-task-status
    p_set = subparsers.add_parser("set-task-status", help="Update task status")
    p_set.add_argument("--id", required=True, help="Task ID (e.g. T1)")
    p_set.add_argument("--status", required=True, choices=list(ALLOWED_STATUSES), help="New status")
    p_set.set_defaults(func=cmd_set_task_status)

    # transition-task
    p_trans = subparsers.add_parser(
        "transition-task", help="Mark one task completed and start the next"
    )
    p_trans.add_argument("--done", required=True, help="Task ID to mark completed (e.g. T1)")
    p_trans.add_argument("--next", required=True, help="Task ID to mark in_progress (e.g. T2)")
    p_trans.set_defaults(func=cmd_transition_task)

    # add-task
    p_add = subparsers.add_parser("add-task", help="Add a new task to active milestone")
    p_add.add_argument("--id", required=True, help="Task ID (e.g. T5)")
    p_add.add_argument("--title", required=True, help="Concise task title")
    p_add.add_argument("--status", default="planned", choices=list(ALLOWED_STATUSES))
    p_add.add_argument("--files", required=True, help="Comma-separated target file paths")
    p_add.add_argument("--contracts", required=True, help="Comma-separated contracts")
    p_add.add_argument("--validation", required=True, help="Verification command")
    p_add.set_defaults(func=cmd_add_task)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
