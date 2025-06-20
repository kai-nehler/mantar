#' Dummy data sets for illustration purposes in the mantar package
#'
#' These two \strong{simulated data sets} are provided for illustration purposes.
#' They are based on a sparse psychological network structure with a single underlying construct.
#' The column names represent core properties of neuroticism but are purely made up to make the example more illustrative.
#' \itemize{
#'   \item \strong{mantar_dummy_full}: A complete data set without missing values.
#'   \item \strong{mantar_dummy_mis}: A version with approximately 30% missing values per column.
#' }
#'
#' @format
#' \describe{
#'   \item{Both data frames}{8 columns; rows: 400 (`mantar_dummy_full`) and 600 (`mantar_dummy_mis`)}
#'   \item{Columns}{\describe{
#'     \item{Emotional reactivity}{Tending to feel emotions strongly in response to life events.}
#'     \item{Tendency toward worry}{Being more likely to feel concerned or uneasy.}
#'     \item{Sensitivity to stress}{Feeling more stressed in challenging or uncertain situations.}
#'     \item{Self-awareness}{Being conscious of one’s own feelings and how they shift.}
#'     \item{Moodiness}{Experiencing occasional changes in mood.}
#'     \item{Cautiousness}{Being careful and thinking ahead about possible negative outcomes.}
#'     \item{Thoughtfulness about future challenges}{Reflecting on what might go wrong and preparing for it.}
#'     \item{Responsiveness to criticism}{Being affected by others’ feedback or disapproval.}
#'   }}
#' }
#'
#' @examples
#' # Load the data sets
#' data(mantar_dummy_full)
#' data(mantar_dummy_mis)
#'
#' # View the first few rows of each data set
#' head(mantar_dummy_full)
#' head(mantar_dummy_mis)
#'
#' @keywords data sets
"mantar_dummy_full"

#' @rdname mantar_dummy_full
"mantar_dummy_mis"
