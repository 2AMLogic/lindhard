# book

Source of the lindhard book: the physics manual (one page per model) and the
user guide. It is published at <https://2amlogic.github.io/lindhard/> on every
merge to `main` (`.github/workflows/book.yml`).

Build it locally with [mdBook](https://rust-lang.github.io/mdBook/):

```sh
mdbook build book      # renders into book/book/
mdbook serve book      # live preview
```

CI also checks that every model enum variant is named in the manual
(`validation/check_manual_coverage.py`) and that the user guide's walkthrough
quotes the output of an actual run (`validation/book_walkthrough.py`).
