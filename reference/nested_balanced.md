# Nested, balanced, and spatially well-spread stratified samples

Draws a sequence of **nested** samples per stratum that are both
probability-balanced (via BalancedSampling) and well-spread in space.
The largest sample is drawn first; each subsequent (smaller) sample is
drawn *from the previous sample only*, but balanced against the
auxiliary variables at the level of the whole population, while spatial
spreading is computed based on the coordinated of the included
population units only.

## Usage

``` r
nested_balanced(
  samplingFrame,
  n_seq,
  id_col = "ID",
  stratum_col = "stratum",
  easting_col = "Easting",
  northing_col = "Northing",
  area_col = "area",
  xbal_formula = ~1,
  exclude_offset = 1e+06,
  return_dataframe = FALSE,
  out_name = "mysample",
  safetyNumber = 1,
  auto_n = TRUE,
  min_step = 10
)
```

## Arguments

- samplingFrame:

  A `data.frame` (your population of sampling units, e.g. `mypop`).

- n_seq:

  Numeric vector of desired **per-stratum** sample sizes in decreasing
  order, e.g. `c(100, 80, 60, 40, 20)`. The first element is the largest
  (top-level) sample; each subsequent element is a nested subsample of
  the previous. If auto_n is TRUE (which is the default), then only two
  values are needed; a start and a stop values (i.e. c(100, 10)).

- id_col:

  Name of the ID column in `samplingFrame`. Default `"ID"`. Needs to be
  unique to each population unit.

- stratum_col:

  Name of the stratum column in `samplingFrame`. Default `"stratum"`.

- easting_col:

  Name of the easting (x) coordinate column. Default `"Easting"`.

- northing_col:

  Name of the northing (y) coordinate column. Default `"Northing"`.

- area_col:

  Name of the area/size measure column used to set inclusion
  probabilities within stratum (your `area`). Default `"area"`. If you
  are not using sampling with probabilities proportional to size, then
  `area_col` might just be the population or strata size *N*.

- xbal_formula:

  A right-hand-side formula for the balancing variables, evaluated on
  `samplingFrame`, e.g. `~ aux1 + aux2 - 1`. Default `~ 1` (no balancing
  auxiliaries beyond spreading/strata).

- exclude_offset:

  Numeric offset added to the max of each coordinate to push excluded
  units *far away* for the spreading step. Default `1e6`.

- return_dataframe:

  Logical; if `TRUE`, also returns a filtered sampling frame, with all
  columns, instead of just the populatio units ID's. Default `FALSE`.

- safetyNumber:

  Numeric; constant which is subtracted from difference between the
  current n (for a given sample) and the remaining area, in order to
  keep probabilities above 1. Default is 1, and the parameter normally
  don't need to be altered.

- auto_n:

  Logical; if the function should try to find the smallest possible step
  (but still larger than 'min_step') that still produces inclusion
  probabilities \< 1. Default is TRUE, and then the n_seq can be set to
  just two numbers, a start and a stop value. If FALSE, then the n is
  given ny 'n_seq'.

- min_step:

  Numeric;, the minimal reduction in n between subsequent samples.
  Default is 10.

- mysample:

  Character; text to be prepended to the sample names, followed by the
  sample size. Default to 'mysample' which will name the output like
  'mysample_n' where n is the sample size.

## Value

A named list of lists. For each `n` you get an element named
`mysample{n}` (e.g. `wetlands_100`, `wetlands_80`, ...), each containing
a list of population ID's, and a list of inclusion probabilities, which
are the accumulated probabilities. This means that if the probability of
a unit in the initial sample is 0.1, and in the second sample it was
0.9, then the accumulated probaili is 0.09. If `return_dataframe = TRUE`
a `data.frame` with the rows of `samplingFrame` that were selected.

## Details

For each requested sample size `n` (per stratum), the function sets
inclusion probabilities \\\pi_i = n \cdot a_i / \sum\_{h} a_i\\ within
each stratum using the `area_col` values \\a_i\\. For the first draw
this sum is taken over the entire stratum; for subsequent draws it is
taken over the *previous* sample only, so that samples are nested.
Spatial spreading uses `Xspread` (Easting/Northing). Units not available
in the current step are assigned coordinates far outside the study area
to avoid influencing the spread.

Requires BalancedSampling.

## Examples

``` r
if (FALSE) { # \dontrun{
# Suppose mypop has: ID, stratum, Easting, Northing, area2, aux1, aux2
# We want nested samples of sizes 100, 80, 60, 40, 20 per stratum:

data(mypop)
out <- nested_balanced(
  data = mypop,
  n_seq = c(100, 80, 60, 40, 20),
  id_col = "ID",
  stratum_col = "stratum",
  easting_col = "Easting",
  northing_col = "Northing",
  area_col = "area2",
  xbal_formula = ~ aux1 + aux2 - 1,
  auto_n = FALSE,
  mysample = "sample"
)

# Access the largest and a nested subset:
sample_100 <- out$mysample_100
sample_60  <- out$mysample_60
} # }
```
