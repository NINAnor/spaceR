#' Nested, balanced, and spatially well-spread stratified samples
#'
#' Draws a sequence of **nested** samples per stratum that are both
#' probability-balanced (via \pkg{BalancedSampling}) and well-spread in space.
#' The largest sample is drawn first; each subsequent (smaller) sample is drawn
#' *from the previous sample only*, but balanced against the auxiliary variables at
#' the level of the whole population, while spatial spreading is computed based on 
#' the coordinated of the included population units only. 
#'
#' @param samplingFrame A `data.frame` (your population of sampling units, e.g. `mypop`).
#' @param n_seq Numeric vector of desired **per-stratum** sample sizes in
#'   decreasing order, e.g. `c(100, 80, 60, 40, 20)`. The first element is the
#'   largest (top-level) sample; each subsequent element is a nested subsample
#'   of the previous. If auto_n is TRUE (which is the default), then only two values are needed; a
#' start and a stop values (i.e. c(100, 10)).
#' @param id_col Name of the ID column in `samplingFrame`. Default `"ID"`. Needs to be unique to each population unit. 
#' @param stratum_col Name of the stratum column in `samplingFrame`. Default `"stratum"`.
#' @param easting_col Name of the easting (x) coordinate column. Default `"Easting"`.
#' @param northing_col Name of the northing (y) coordinate column. Default `"Northing"`.
#' @param area_col Name of the area/size measure column used to set inclusion
#'   probabilities within stratum (your `area`). Default `"area"`. If you are not using 
#'   sampling with probabilities proportional to size, then `area_col` might just be 
#'   the population or strata size _N_.
#' @param xbal_formula A right-hand-side formula for the balancing variables,
#'   evaluated on `samplingFrame`, e.g. `~ aux1 + aux2 - 1`. Default `~ 1` (no balancing
#'   auxiliaries beyond spreading/strata).
#' @param exclude_offset Numeric offset added to the max of each coordinate to
#'   push excluded units *far away* for the spreading step. Default `1e6`.
#' @param return_dataframe Logical; if `TRUE`, also returns a filtered sampling frame, with all columns, 
#'   instead of just the populatio units ID's. Default `FALSE`.
#' @param mysample Character; text to be prepended to the sample names, followed by the sample size. 
#' Default to 'mysample' which will name the output like 'mysample_n' where n is the sample size.
#' @param safetyNumber Numeric; constant which is subtracted from difference between the current n
#' (for a given sample) and the remaining area, in order to keep probabilities above 1. Default is 1, 
#' and the parameter normally don't need to be altered.
#' @param auto_n Logical; if the function should try to find the smallest possible step (but still larger than 
#' 'min_step') that still produces inclusion probabilities < 1. Default is TRUE, and then the n_seq can be
#' set to just two numbers, a start and a stop value. If FALSE, then the n is given ny 'n_seq'. 
#' @param min_step Numeric;, the minimal reduction in n between subsequent samples. Default is 10.
#'
#' @details
#' For each requested sample size `n` (per stratum), the function sets inclusion
#' probabilities \eqn{\pi_i = n \cdot a_i / \sum_{h} a_i} within each stratum
#' using the `area_col` values \eqn{a_i}. For the first draw this sum is taken
#' over the entire stratum; for subsequent draws it is taken over the *previous*
#' sample only, so that samples are nested. Spatial spreading uses \code{Xspread}
#' (Easting/Northing). Units not available in the current step are assigned
#' coordinates far outside the study area to avoid influencing the spread.
#'
#' Requires \pkg{BalancedSampling}.
#'
#' @return A named list of lists. For each `n` you get an element named
#'   `mysample{n}` (e.g. `wetlands_100`, `wetlands_80`, ...), each containing a list of population ID's,
#'   and a list of inclusion probabilities, which are the accumulated probabilities. This means that if the 
#' probability of a unit in the initial sample is 0.1, and in the second sample it was 0.9, then the accumulated 
#' probaili is 0.09. If `return_dataframe = TRUE` a `data.frame` with the rows
#'   of `samplingFrame` that were selected. 
#'
#' @examples
#' \dontrun{
#' # Suppose mypop has: ID, stratum, Easting, Northing, area2, aux1, aux2
#' # We want nested samples of sizes 100, 80, 60, 40, 20 per stratum:
#' 
#' data(mypop)
#' out <- nested_balanced(
#'   data = mypop,
#'   n_seq = c(100, 80, 60, 40, 20),
#'   id_col = "ID",
#'   stratum_col = "stratum",
#'   easting_col = "Easting",
#'   northing_col = "Northing",
#'   area_col = "area2",
#'   xbal_formula = ~ aux1 + aux2 - 1,
#'   auto_n = FALSE,
#'   mysample = "sample"
#' )
#'
#' # Access the largest and a nested subset:
#' sample_100 <- out$mysample_100
#' sample_60  <- out$mysample_60
#' }
#'
#' @importFrom stats model.matrix
#' @export
nested_balanced <- function(
    samplingFrame,
    n_seq,
    id_col = "ID",
    stratum_col = "stratum",
    easting_col = "Easting",
    northing_col = "Northing",
    area_col = "area",
    xbal_formula = ~ 1,
    exclude_offset = 1e6,
    return_dataframe = FALSE,
    out_name = "mysample",
    safetyNumber = 1,
    auto_n = TRUE,
    min_step = 10
) {
  if (!requireNamespace("BalancedSampling", quietly = TRUE)) {
    stop("Package 'BalancedSampling' is required but not installed.")
  }
  
  if (!is.data.frame(samplingFrame)) stop("'samplingFrame' must be a data.frame.")
  if (length(n_seq) < 1) stop("'n_seq' must have at least one value.")
  if (is.unsorted(rev(n_seq))) {
    warning("It's recommended that 'n_seq' be in decreasing order (largest first).")
  }
  
  # Pull columns
  id        <- samplingFrame[[id_col]]
  stratum   <- samplingFrame[[stratum_col]]
  easting   <- samplingFrame[[easting_col]]
  northing  <- samplingFrame[[northing_col]]
  area      <- samplingFrame[[area_col]]
  
  if (anyNA(id)) stop("ID column contains NA.")
  if (anyDuplicated(id) != 0) stop("The ID columns contains duplicates.")
  if (anyNA(stratum)) stop("Stratum column contains NA.")
  if (anyNA(easting) || anyNA(northing)) stop("Coordinate columns contain NA.")
  if (anyNA(area)) stop("Area column contains NA.")
  if (any(area < 0)) stop("Area column has negative values.")
  
  # Balancing variables (Xbal)
  Xbal <- stats::model.matrix(
    stats::as.formula(paste("~", as.character(xbal_formula)[2])),
    data = samplingFrame
  )
  
  # Helper: compute per-stratum totals for a vector (e.g. area) but
  # only over a subset 'keep_idx' (logical).
  strat_totals <- function(x, keep_idx) {
    # sum by stratum considering only units where keep_idx is TRUE
    tapply(ifelse(keep_idx, x, 0), stratum, sum)
  }
  # Helper: get the per strata n
  get_n <- function(x) {
      x |>
        group_by(stratum) |>
        summarise(n = n()) |> 
        pull(n)
    }

  
  # Precompute a "far away" coordinate to nullify excluded units' influence
  far_x <- max(easting, na.rm = TRUE) + exclude_offset
  far_y <- max(northing, na.rm = TRUE) + exclude_offset
  
  # Prepare a container for outputs
  out <- list()
  
  # State: which units are eligible at current step
  # For the first (largest) sample: all units are eligible
  eligible <- rep(TRUE, nrow(samplingFrame))
  
  # set prob_last as all 1's for the first iteration of the loop
  prob_last <- rep(1, nrow(samplingFrame))
  
  # make seq longer (denser) of using auto_n
  if(auto_n) {
    n_seq2 <- seq(max(n_seq), min(n_seq), -1)
  } 
  
  # Iterate over requested nested sizes
  for (n_now in n_seq2) {
    # Remaining area per stratum among eligible units
    rem_area_by_stratum <- strat_totals(area, keep_idx = eligible)

    # Map each unit's stratum to the remaining area value
    rem_area_unit <- rem_area_by_stratum[as.character(stratum)]
    if (any(rem_area_unit <= 0)) {
      stop("Some strata have zero remaining area among eligible units; cannot compute probabilities.")
    }
    
    
    if(n_now == n_seq2[1]) {
      n_ <- n_now
    } else {
      # Auto n - find the highest possible next n that results in probabilities <1
      max_n <- n_last-ceiling(max(n_last - rem_area_by_stratum)) - safetyNumber
      if(n_last-max_n < min_step) {max_n <- n_last-min_step}
      n_ <- ifelse(auto_n, max_n, n_now)
    }


    # Inclusion probabilities for this draw (0 for ineligible)
    prob <- ifelse(eligible, n_ * area / rem_area_unit, 0)
    if (any(prob >= 1)) {
      stop("Some population units have >1 probability of being selected.")
    }
    
    # Build Xspread: real coords for eligible units; far away for others
    Xspread <- cbind(
      ifelse(eligible, easting, far_x),
      ifelse(eligible, northing, far_y)
    )
    
    # Draw sample indices (row positions in 'samplingFrame')
    sel_idx <- BalancedSampling::lcubestratified(
      prob         = prob,
      Xspread      = Xspread,
      Xbal         = Xbal,
      integerStrat = stratum
    )
    
    # Keep only selected rows, and append probabilities
    current_sample <- samplingFrame[sel_idx, , drop = FALSE] |>
      cbind(current_pi = prob[sel_idx],
            accumulated_pi = prob[sel_idx]*prob_last[sel_idx])
        
    n_string <- get_n(current_sample)

    if(length(unique(n_string)) != 1) {
      warning("The n differs amongst strata in one of the samples.")
    }

    if(unique(n_string)[1] != n_) {
      warning("The n for at least one of the samples differs from the prescribed n.")
    }

    if (anyDuplicated(current_sample[[id_col]]) != 0) stop("The ID columns on one of the samples duplicates.")
      
    # Save result
    name_now <- paste(out_name, n_, sep = "_")
    if (return_dataframe) {
      out[[name_now]] <- current_sample
      
    } else {
      out[[name_now]] <- list(
        ID = current_sample[[id_col]], 
        prob = current_sample$accumulated_pi
      )
    }
    
    # For the next (smaller) sample, only the current selection remains eligible
    # Reset 'eligible' to FALSE and flag selected rows TRUE
    eligible <- rep(FALSE, nrow(samplingFrame))
    eligible[sel_idx] <- TRUE

    if(n_ < min(n_seq)) break

    n_last <- n_
    prob_last <- prob
    
  }
  
  return(out)
}
