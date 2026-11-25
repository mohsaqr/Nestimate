# Coded turns of human-AI vibe coding sessions

`human_long` and `ai_long` are a coded dataset of the turns of 429
human-AI pair-programming ("vibe coding") sessions across 34 projects,
held between 5 February and 8 March 2026: `human_long` has the codes of
the human turns and `ai_long` those of the AI turns (428 of the sessions
have AI turns). Both are in long format, one row per code; they share
their columns and their sessions, and interleave on `order_in_session`.

## Usage

``` r
human_long

ai_long
```

## Format

Two data frames with the same 9 columns:

- message_id:

  Integer identifier of the message (turn); the codes of one message
  share it, and no identifier occurs in both data frames.

- project:

  Character project identifier, `Project_1` to `Project_34`.

- session_id:

  Character session identifier, a 12-character hash.

- timestamp:

  Integer Unix time (seconds) of midnight, UTC+2, on `session_date`. It
  has day resolution, so the rows of one day share it;
  `build_network(time = "timestamp")` keeps tied rows in row order,
  which follows `order_in_session`. With the default `time_threshold`
  (900 seconds), each change of day also starts a new sequence, so the
  84 sessions that run over more than one day are split (526 sequences
  instead of 429); to keep every session whole, use
  `time = "order_in_session"`.

- session_date:

  Character date of the row (`YYYY-MM-DD`); 84 sessions run over more
  than one day.

- code:

  Character interaction code. The two data frames use disjoint code
  sets: 9 labels in `human_long` (Command, Correct, Frustrate, Inquire,
  Interrupt, Refine, Request, Specify, Verify) and 8 in `ai_long` (Ask,
  Delegate, Execute, Explain, Investigate, Plan, Repair, Report), 17 in
  total.

- cluster:

  Character high-level cluster of the code. `human_long` uses Directive,
  Evaluative and Metacognitive; `ai_long` uses Action, Communication and
  Repair.

- code_order:

  Integer order of the code within its message, starting at 1.

- order_in_session:

  Integer position of the row in its session, counting the rows of both
  data frames together.

|              |                               |
|--------------|-------------------------------|
| `human_long` | Human turns only, 10,796 rows |
| `ai_long`    | AI turns only, 8,551 rows     |

An object of class `data.frame` with 10796 rows and 9 columns.

An object of class `data.frame` with 8551 rows and 9 columns.

## Source

Saqr, M. (2026). Human-AI vibe coding interaction study.
[saqr.me](https://saqr.me/blog/2026/human-ai-interaction-cograph/)

## Examples

``` r
net <- build_network(human_long, method = "tna",
                     action = "code", actor = "session_id",
                     time = "timestamp")
#> Metadata aggregated per session: ties resolved by first occurrence in 'cluster' (20 sessions)
net
#> Transition Network (relative probabilities) [directed]
#>   Weights: [0.018, 0.620]  |  mean: 0.111
#> 
#>   Strongest edges (10 of 81):
#>    from      to      weight
#>    Request   Specify 0.620 
#>    Refine    Specify 0.457 
#>    Command   Specify 0.280 
#>    Correct   Specify 0.278 
#>    Specify   Command 0.263 
#>    Interrupt Command 0.259 
#>    Command   Command 0.235 
#>    Verify    Command 0.224 
#>    Inquire   Command 0.196 
#>    Inquire   Inquire 0.188  
#> 
#>   extract_edges(include_self = TRUE) lists every edge; extract_transition_matrix() gives the weight matrix;
#>   extract_initial_probs() the initial probabilities.
```
