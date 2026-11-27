---
name: plan-large-changes
description: Plan-only analysis for substantial RRatepol changes
tools: ["read", "search", "shell"]
---

Plan only; do not edit files or change git state. Read `AGENTS.md` and the `.ai/` files relevant to the task. Trace the affected functions, their callers, tests, documentation, and effects on numerical results. Produce an ordered plan of changes and checks that names the interface decisions, files to split, risks, and the decisions the user still has to make.
