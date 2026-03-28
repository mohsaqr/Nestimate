# Human-AI Vibe Coding Interaction Data

Coded interaction sequences from 429 human-AI pair programming sessions
across 34 projects. Three coding granularities: code (32 states),
category (17 states), and superclass (6 states).

## Usage

``` r
human_ai

human_ai_cat

human_ai_super

human_detailed

human_cat

human_super

ai_detailed

ai_cat

ai_super

human_wide

ai_wide
```

## Format

Long-format data frames with columns:

- id:

  Integer. Turn index within the session.

- project:

  Character. Project identifier (Project_1 .. Project_34).

- session_id:

  Character. Unique session hash.

- timestamp:

  Character. ISO 8601 timestamp.

- session_date:

  Character. Date of the session (YYYY-MM-DD).

- actor:

  Character. `"Human"` or `"AI"`.

- code:

  Character. Fine-grained action code (32 states).

- category:

  Character. Mid-level category (17 states).

- superclass:

  Character. High-level superclass (6 states).

An object of class `data.frame` with 19347 rows and 9 columns.

An object of class `data.frame` with 19347 rows and 9 columns.

An object of class `data.frame` with 19347 rows and 9 columns.

An object of class `data.frame` with 10796 rows and 9 columns.

An object of class `data.frame` with 10796 rows and 9 columns.

An object of class `data.frame` with 10796 rows and 9 columns.

An object of class `data.frame` with 8551 rows and 9 columns.

An object of class `data.frame` with 8551 rows and 9 columns.

An object of class `data.frame` with 8551 rows and 9 columns.

An object of class `data.frame` with 429 rows and 164 columns.

An object of class `data.frame` with 428 rows and 138 columns.

## Source

Saqr, M. (2026). Human-AI vibe coding interaction study.

## Details

Nine long-format datasets are provided, filtered by actor and named by
granularity level:

|                  |           |                       |
|------------------|-----------|-----------------------|
| **Dataset**      | **Actor** | **Granularity**       |
| `human_ai`       | Both      | code (32 states)      |
| `human_ai_cat`   | Both      | category (17 states)  |
| `human_ai_super` | Both      | superclass (6 states) |
| `human_detailed` | Human     | code (32 states)      |
| `human_cat`      | Human     | category (17 states)  |
| `human_super`    | Human     | superclass (6 states) |
| `ai_detailed`    | AI        | code (32 states)      |
| `ai_cat`         | AI        | category (17 states)  |
| `ai_super`       | AI        | superclass (6 states) |

Two wide-format datasets at category level (rows = sessions, columns =
T1, T2, ...):

|              |                                       |
|--------------|---------------------------------------|
| `human_wide` | Human actions in wide sequence format |
| `ai_wide`    | AI actions in wide sequence format    |

## Examples

``` r
# \donttest{
# Build a transition network from human category sequences
net <- build_network(human_wide, method = "relative")

# Use the edge list directly
head(human_ai_edges)
#>            from            to weight   session_id session   project order
#> 1       Context        Direct      1 0086cabebd15       1 Project_7     1
#> 2        Direct Specification      1 0086cabebd15       1 Project_7     2
#> 3 Specification     Interrupt      1 0086cabebd15       1 Project_7     3
#> 4     Interrupt      Delegate      1 0086cabebd15       1 Project_7     4
#> 5      Delegate          Plan      1 0086cabebd15       1 Project_7     5
#> 6          Plan  Verification      1 0086cabebd15       1 Project_7     6
#>                  timepoint from_actor to_actor from_category to_category
#> 1 2026-03-05T11:32:52.057Z      Human    Human       Specify     Command
#> 2 2026-03-05T11:32:52.057Z      Human    Human       Command     Specify
#> 3 2026-03-05T11:32:52.057Z      Human    Human       Specify   Interrupt
#> 4 2026-03-05T11:32:52.068Z      Human       AI     Interrupt    Delegate
#> 5 2026-03-05T11:32:56.945Z         AI       AI      Delegate        Plan
#> 6 2026-03-05T11:32:56.945Z         AI    Human          Plan      Verify
#>   from_superclass to_superclass
#> 1       Directive     Directive
#> 2       Directive     Directive
#> 3       Directive Metacognitive
#> 4   Metacognitive        Action
#> 5          Action        Repair
#> 6          Repair    Evaluative
# }
```
