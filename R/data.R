#' Dummy data sets for illustration purposes in the mantar package
#'
#' These \strong{simulated data sets} are provided for illustration purposes.
#' They are based on a sparse psychological network structure with a single underlying construct.
#' The column names represent core properties of neuroticism but are purely made up to make the example more illustrative.
#' \itemize{
#'   \item \strong{mantar_dummy_full}: A complete data set of continuous variables without missing values.
#'   \item \strong{mantar_dummy_mis}: A data set of continuous variables with approximately 30% missing values in each column.
#'   \item \strong{mantar_dummy_full_cat}: A complete data set where variables are ordered categorical without missing values.
#' }
#'
#' @format
#' \describe{
#'   \item{Both data frames}{8 columns; rows: 400 (`mantar_dummy_full` and `mantar_dummy_full_cat`) and 600 (`mantar_dummy_mis`)}
#'   \item{Columns}{\describe{
#'     \item{EmoReactivity}{Tending to feel emotions strongly in response to life events.}
#'     \item{TendWorry}{Being more likely to feel concerned or uneasy.}
#'     \item{StressSens}{Feeling more stressed in challenging or uncertain situations.}
#'     \item{SelfAware}{Being conscious of one’s own feelings and how they shift.}
#'     \item{Moodiness}{Experiencing occasional changes in mood.}
#'     \item{Cautious}{Being careful and thinking ahead about possible negative outcomes.}
#'     \item{ThoughtFuture}{Reflecting on what might go wrong and preparing for it.}
#'     \item{RespCriticism}{Being affected by others’ feedback or disapproval.}
#'   }}
#' }
#'
#' @examples
#' # Load the data sets
#' data(mantar_dummy_full)
#' data(mantar_dummy_mis)
#' data(mantar_dummy_full_cat)
#'
#' # View the first few rows of each data set
#' head(mantar_dummy_full)
#' head(mantar_dummy_mis)
#' head(mantar_dummy_full_cat)
#'
#' @keywords data sets
"mantar_dummy_full"

#' @rdname mantar_dummy_full
"mantar_dummy_mis"

#' @rdname mantar_dummy_full
"mantar_dummy_full_cat"
