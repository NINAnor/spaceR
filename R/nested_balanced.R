#' Nested, balanced, and spatially well-spread stratified samples
#'
#' Draws a sequence of **nested** samples per stratum that are both
#' probability-balanced (via \pkg{BalancedSampling}) and well-spread in space.
#' The first sample can also use probabilities proportional to size (pps).
#' The largest sample is drawn first; each subsequent (smaller) sample is drawn
#' *from the previous sample only*, but balanced against the auxiliary variables at
#' the level of the whole population, while spatial spreading is computed based on 
#' the coordinated of the included population units only. 
#'
#' @param samplingFrame A `data.frame` (your population of sampling units, e.g. `mypop`).
#' @param n_seq Numeric vector of desired **per-stratum** sample sizes in
#'   decreasing order, e.g. `c(100, 80, 60, 40, 20)`. The first element is the
#'   largest (top-level) sample; each subsequent element is a nested subsample
#'   of the previous. 
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
#'   instead of just the population units ID's and the probabilities. Default `FALSE`.
#' @param out_name Character; text to be prepended to the sample names, followed by the sample size. 
#' Default to 'mysample' which will name the output like 'mysample_n' where n is the sample size.
#'
#' @details
#' For each requested sample size `n` in n_seq (per stratum), the function returns
#' a balanced, well-spread and stratified sample. All samples are nested,
#' so that samples with small n is a (balanced) subset of samples with greater n.
#' The largest sample (the first iteration) can use pps, and sets inclusion
#' probabilities \eqn{\pi_i = n \cdot a_i / \sum_{h} a_i} within each stratum
#' using the `area_col` values \eqn{a_i}. Subsequent samples use equal probabilities.
#' Spatial spreading uses \code{Xspread} (Easting/Northing) and is based on the position of 
#' the remaining population units only. The balancing is, however, always done against the 
#' entire population.
#'
#' Requires \pkg{BalancedSampling}.
#'
#' @return A named list of lists. For each `n` you get an element named
#'   `mysample{n}` (e.g. `wetlands_100`, `wetlands_80`, ...), each containing a list of population ID's,
#'   and a list of inclusion probabilities, which are the accumulated probabilities and the initial probabilities. 
#' This means that if the 
#' probability of a unit in the initial sample is 0.1, and in the second sample it was 0.9, then the accumulated 
#' probaili is 0.09. If `return_dataframe = TRUE` a `data.frame` with the rows
#' of `samplingFrame` that were selected. 
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
#'   out_name = "mysample"
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
    quiet = FALSE
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
  
  # Calculate the probabilities for the first iteration, using pps
    # Area per stratum among eligible units
      rem_area_by_stratum <- strat_totals(area, keep_idx = eligible)
      
    # Map each unit's stratum to the remaining area value
      rem_area_unit <- rem_area_by_stratum[as.character(stratum)]
        if (any(rem_area_unit <= 0)) {
          stop("Some strata have zero area; cannot compute probabilities.")
        }
      
      prob_init <- ifelse(eligible, n_seq[1] * area / rem_area_unit, 0)
      prob_last <- rep(1, nrow(samplingFrame))
  # Iterate over requested nested sizes
  for (n_now in n_seq) {
    
    if (!quiet) {
      message(paste("Sampling with n = ", n_now))
    }

    # Inclusion probabilities for this draw (0 for ineligible)
    if(n_now == n_seq[1]) {
      prob <- prob_init
    } else {       
      # Remaining n per stratum among eligible units
      rem_n_by_stratum <- tapply(as.numeric(eligible), stratum, sum)      
      # Map each unit's stratum to the remaining n
      rem_n_unit <- rem_n_by_stratum[as.character(stratum)]
        if (any(rem_n_unit <= 0)) {
          stop("Some strata have zero remaining population units; cannot compute probabilities.")
        }
      prob <- ifelse(eligible, n_now / rem_n_unit, 0)
    }

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
      cbind(
        initial_pi = prob_init[sel_idx],
        current_pi = prob[sel_idx],
        previous_pi = prob_last[sel_idx],
        previous_accumulated_pi = if(n_now == n_seq[1]) {
          prob_init[sel_idx] } else {
          pacc_pi[sel_idx]
          },
        accumulated_pi = if(n_now == n_seq[1]) {
          prob_init[sel_idx]
        } else {
        pacc_pi[sel_idx]*prob[sel_idx]   
        }
      )
        
    n_string <- get_n(current_sample)

    if(length(unique(n_string)) != 1) {
      warning("The n differs amongst strata in one of the samples.")
    }

    if(unique(n_string)[1] != n_now) {
      warning("The n for at least one of the samples differs from the prescribed n.")
    }

    if(anyDuplicated(current_sample[[id_col]]) != 0) stop("The ID columns on one of the samples duplicates.")
      
    # Save result
    name_now <- paste(out_name, n_now, sep = "_")
    if (return_dataframe) {
      out[[name_now]] <- current_sample
      
    } else {
      out[[name_now]] <- list(
        ID = current_sample[[id_col]], 
        initial_prob = current_sample$initial_pi,
        prob = current_sample$accumulated_pi
      )
    }
    
    # For the next (smaller) sample, only the current selection remains eligible
    # Reset 'eligible' to FALSE and flag selected rows TRUE
    eligible <- rep(FALSE, nrow(samplingFrame))
    eligible[sel_idx] <- TRUE


    pacc_pi <- prob_last*prob
    prob_last <- prob
        }
  
  return(out)
}
