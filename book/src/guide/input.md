# TOML input reference

An input file is TOML with these tables: `[beam]`, `[materials]` (optional),
`[target]`, `[physics]`, `[stopping]` (optional), `[run]` and `[tally]`
(optional). The physics behind each `[physics]` choice is in the
[physics manual](../models/overview.md): [potentials](../models/potentials.md)
and [screening lengths](../models/screening-lengths.md),
[electronic stopping](../models/stopping.md), the
[free-path conventions and weak collisions](../models/bca.md), and
[user stopping tables](../models/stopping-user-tables.md). The
[first run](first-run.md) walks through a complete file, and the
repository's
[`examples/`](https://github.com/2AMLogic/lindhard/tree/main/examples)
directory has more.

This page is the input section of
[`docs/cli.md`](https://github.com/2AMLogic/lindhard/blob/main/docs/cli.md),
included here so there is one copy.

{{#include ../../../docs/cli.md:input}}
