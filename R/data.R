#' BCG Vaccine Meta-Analysis Example Data
#'
#' A classic meta-analysis of 13 trials examining the effectiveness of
#' the BCG vaccine against tuberculosis. Suitable for demonstrating
#' MAFI publication bias analysis.
#'
#' @format A data frame with 13 rows and 5 variables:
#' \describe{
#'   \item{study}{Study author(s)}
#'   \item{year}{Publication year}
#'   \item{yi}{Log risk ratio effect size}
#'   \item{vi}{Variance of effect size}
#'   \item{sei}{Standard error of effect size}
#' }
#'
#' @source Originally from: Colditz et al. (1994). Efficacy of BCG vaccine
#'   in the prevention of tuberculosis. JAMA, 271, 698-702.
#'
#' @examples
#' data(mafi_example_bcg)
#' result <- mafi(mafi_example_bcg$yi, mafi_example_bcg$vi)
#' print(result)
#'
"mafi_example_bcg"
