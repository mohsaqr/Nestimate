# Human-AI Vibe Coding Interaction Data (Long Format)

Coded interaction sequences from 429 human-AI pair programming sessions
across 34 projects, in long format.

## Usage

``` r
human_long

ai_long
```

## Format

Data frames in long format with 9 columns:

- message_id:

  Integer. Turn index.

- project:

  Character. Project identifier (Project_1 .. Project_34).

- session_id:

  Character. Unique session hash.

- timestamp:

  Integer. Unix timestamp for ordering.

- session_date:

  Character. Date of the session (YYYY-MM-DD).

- code:

  Character. Interaction code (17 action labels).

- cluster:

  Character. High-level cluster: Action, Communication, Directive,
  Evaluative, Metacognitive, or Repair.

- code_order:

  Integer. Order of the code within the session.

- order_in_session:

  Integer. Absolute turn order within the session.

|              |                               |
|--------------|-------------------------------|
| `human_long` | Human turns only, 10,796 rows |
| `ai_long`    | AI turns only, 8,551 rows     |

An object of class `data.frame` with 10796 rows and 9 columns.

An object of class `data.frame` with 8551 rows and 9 columns.

## Source

Saqr, M. (2026). Human-AI vibe coding interaction study.
<https://saqr.me/blog/2026/human-ai-interaction-cograph/>

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
#>   Weight matrix:
#>             Command Correct Frustrate Inquire Interrupt Refine Request Specify
#>   Command     0.235   0.088     0.055   0.066     0.035  0.038   0.155   0.280
#>   Correct     0.093   0.091     0.138   0.057     0.054  0.112   0.120   0.278
#>   Frustrate   0.103   0.115     0.176   0.075     0.047  0.171   0.111   0.125
#>   Inquire     0.196   0.126     0.098   0.188     0.094  0.062   0.084   0.106
#>   Interrupt   0.259   0.081     0.094   0.123     0.102  0.080   0.069   0.122
#>   Refine      0.058   0.075     0.072   0.047     0.042  0.086   0.146   0.457
#>   Request     0.096   0.019     0.044   0.067     0.051  0.033   0.039   0.620
#>   Specify     0.263   0.058     0.081   0.072     0.167  0.074   0.070   0.172
#>   Verify      0.224   0.079     0.164   0.118     0.043  0.097   0.116   0.083
#>             Verify
#>   Command    0.048
#>   Correct    0.056
#>   Frustrate  0.077
#>   Inquire    0.044
#>   Interrupt  0.069
#>   Refine     0.018
#>   Request    0.032
#>   Specify    0.043
#>   Verify     0.077 
#> 
#>   Initial probabilities:
#>   Specify       0.679  ████████████████████████████████████████
#>   Command       0.175  ██████████
#>   Request       0.055  ███
#>   Interrupt     0.030  ██
#>   Correct       0.019  █
#>   Refine        0.019  █
#>   Frustrate     0.013  █
#>   Inquire       0.008  
#>   Verify        0.002  
```
