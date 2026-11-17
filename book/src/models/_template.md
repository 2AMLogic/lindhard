<!--
Template for a physics-manual page. Copy it to book/src/models/<model>.md,
add the page to book/src/SUMMARY.md, and fill in every section. This file is
not part of the book and does not count for the coverage check.

Rules (CONTRIBUTING.md, CLAUDE.md):
- Every equation and every constant cites a published source, listed in this
  page's own References section.
- Every number that enters the code as data has a row in
  docs/data-provenance.md; say here how far it has been verified.
- Tables from any source (SRIM, ICRU, NIST and others) are admissible as
  cited data; label their origin, and say when a model is fitted to a code's
  output. No code from Tier B or Tier C programs. Nothing that is not already
  published physics.
- If the model is selected by an enum, name every variant in its Rust
  spelling, `Enum::Variant`, and, if the TOML input exposes it, in its TOML
  spelling, `key = "kebab-case"`. Add the enum to ENUMS in
  validation/check_manual_coverage.py; CI then fails if a variant is missing.
- Math: inline \\( ... \\), display \\[ ... \\]; a backslash before
  punctuation must be doubled (\\, for a thin space).
-->

# Model name

Code: `lindhard/src/<path>.rs` (`TypeName`).

## Model

What the model describes, in one paragraph. Then the equations:

\\[ \text{equation} \\]

with every symbol defined, and the source of each equation and constant.

## Selecting it

| TOML | Rust |
|---|---|
| `key = "variant-name"` | `EnumName::VariantName` |

Defaults, and what is library-only.

## Assumptions

- The approximations the model makes.

## Validity

The range of energies, materials or geometries over which the model is meant
to be used, and what happens outside it (refusal, warning, silent use).

## Verification status

What has been checked against what (primary paper, secondary source, test,
validation level), and what has not. Link the `docs/data-provenance.md` row.

## References

- Author, Journal volume, page (year), doi.
