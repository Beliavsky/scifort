# Claude Code instructions

`AGENTS.md` is the authoritative repository instruction file. Read and follow
it before making any change.

Additional Claude Code requirements:

- Inspect `CODE_PROVENANCE.md` before translating or adapting any algorithm.
- Ask for no license assumptions. If a source license cannot be verified, do
  not use that implementation.
- Keep edits scoped and review the complete diff before finishing.
- Run `fpm test` after implementation changes. Do not describe unrun tests as
  passing.
- Preserve public procedure names and `bind(c)` symbol names unless the task
  explicitly authorizes an API or ABI break.
- Do not generate Python, R, MATLAB, or Octave wrappers that bypass the common
  C ABI merely because a compiler-specific shortcut is easier.
