# ChatGPT Self-Regulated Learning Scale Scores

Scores on five self-regulated learning (SRL) scales for 1,000 responses
that ChatGPT generated to a validated SRL questionnaire. The data are
part of a larger set comparing the responses of seven large language
models with human norms.

## Usage

``` r
chatgpt_srl
```

## Format

A data frame with 1,000 rows and 5 columns, one row per response; each
column is the numeric mean of one scale's items (13 items for `CSU`, 4
for `TA` and 9 for each of the others), with observed values from 2 to
6.5 and no missing values:

- CSU:

  Cognitive Strategy Use.

- IV:

  Intrinsic Value.

- SE:

  Self-Efficacy.

- SR:

  Self-Regulation.

- TA:

  Test Anxiety.

## Source

Vogelsmeier, L.V.D.E., Oliveira, E., Misiejuk, K., Lopez-Pernas, S., &
Saqr, M. (2025). Delving into the psychology of Machines: Exploring the
structure of self-regulated learning via LLM-generated survey responses.
*Computers in Human Behavior*, 173, 108769.
[doi:10.1016/j.chb.2025.108769](https://doi.org/10.1016/j.chb.2025.108769)

## Examples

``` r
net <- build_network(chatgpt_srl, method = "glasso",
                     params = list(gamma = 0.5))
net
#> Partial Correlation Network (EBICglasso) [undirected]
#>   Sample size: 1000
#> 
#>   Weight matrix:
#>       CSU IV SE SR TA
#>   CSU   0  0  0  0  0
#>   IV    0  0  0  0  0
#>   SE    0  0  0  0  0
#>   SR    0  0  0  0  0
#>   TA    0  0  0  0  0 
#> 
#>   Predictability (R²):
#>   CSU           0.000  
#>   IV            0.000  
#>   SE            0.000  
#>   SR            0.000  
#>   TA            0.000  
#> 
#>   Gamma: 0.50  |  Lambda: 0.0517
```
