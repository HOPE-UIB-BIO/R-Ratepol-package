# Delegating to other agents

Every agent working in this repository can read `AGENTS.md` and `.ai/` itself. When delegating, do not paste these files into the prompt; give a focused brief that names the files to read, the concrete scope, the constraints, and the expected output.

- Read-only reviews: `.ai/agents/changes-reviewer.agent.md`.
- Plans for large changes, without editing: `.ai/agents/plan-large-changes.agent.md`.
- An agent writing tests for a new function receives the specification stub (roxygen and empty body) and `.ai/r-functions.md`, not the implementation.
- Choose the cheapest model and lowest effort that can do the task reliably.

Delegated agents never perform git operations, add dependencies, or change numerical behaviour on their own. A required independent review is never satisfied by reviewing one's own work.
