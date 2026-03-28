# Human-AI Vibe Coding Edge List

Non-aggregated edge list of all consecutive action transitions across
429 sessions. Each row is one transition from one action to the next
within a session.

## Usage

``` r
human_ai_edges
```

## Format

A data frame with 18,918 rows and 14 columns:

- from:

  Character. Source action (code level).

- to:

  Character. Target action (code level).

- weight:

  Integer. Always 1 (non-aggregated).

- session_id:

  Character. Unique session hash.

- session:

  Integer. Numeric session ID (1..429).

- project:

  Character. Project identifier.

- order:

  Integer. Position of the transition within the session.

- timepoint:

  Character. ISO 8601 timestamp of the source action.

- from_actor:

  Character. Actor of the source action.

- to_actor:

  Character. Actor of the target action.

- from_category:

  Character. Category of the source action.

- to_category:

  Character. Category of the target action.

- from_superclass:

  Character. Superclass of the source action.

- to_superclass:

  Character. Superclass of the target action.

## Source

Saqr, M. (2026). Human-AI vibe coding interaction study.

## Examples

``` r
# \donttest{
# Filter to Human -> AI transitions only
handoffs <- human_ai_edges[
  human_ai_edges$from_actor == "Human" &
  human_ai_edges$to_actor == "AI", ]
# }
```
