#' @title Dummy data sets for illustration purposes in the mantar package
#'
#' @description These simulated data sets are provided for illustration purposes.
#' They are based on a sparse psychological network structure with a single
#' underlying construct. The column names represent core properties of
#' neuroticism but are purely made up to make the example more illustrative.
#'
#' @details
#' The following data sets are available:
#' \itemize{
#'   \item \code{mantar_dummy_full_cont}: A complete data set of continuous
#'   variables without missing values.
#'   \item \code{mantar_dummy_mis_cont}: A data set of continuous variables with
#'   approximately 30% missing values in each column.
#'   \item \code{mantar_dummy_full_cat}: A complete data set where variables are
#'   ordered categorical without missing values.
#'   \item \code{mantar_dummy_mis_cat}: A data set where variables are ordered
#'   categorical with approximately 25% missing values in each column.
#'   \item \code{mantar_dummy_full_mix}: A complete data set with a mix of
#'   continuous and ordered categorical variables without missing values.
#'   \item \code{mantar_dummy_mis_mix}: A data set with a mix of continuous and
#'   ordered categorical variables with approximately 25% missing values in
#'   each column.
#' }
#'
#' @format
#' All data sets are data frames with 400 rows and 8 columns. The columns are:
#' \describe{
#'   \item{EmoReactivity}{Tending to feel emotions strongly in response to life events.}
#'   \item{TendWorry}{Being more likely to feel concerned or uneasy.}
#'   \item{StressSens}{Feeling more stressed in challenging or uncertain situations.}
#'   \item{SelfAware}{Being conscious of one’s own feelings and how they shift.}
#'   \item{Moodiness}{Experiencing occasional changes in mood.}
#'   \item{Cautious}{Being careful and thinking ahead about possible negative outcomes.}
#'   \item{ThoughtFuture}{Reflecting on what might go wrong and preparing for it.}
#'   \item{RespCriticism}{Being affected by others’ feedback or disapproval.}
#' }
#'
#' @name mantar_dummy_data
#' @docType data
#' @aliases mantar_dummy_full_cont mantar_dummy_mis_cont
#'   mantar_dummy_full_cat mantar_dummy_mis_cat
#'   mantar_dummy_full_mix mantar_dummy_mis_mix
#'
#' @examples
#' # Load selected data set
#' data(mantar_dummy_full_cont)
#' data(mantar_dummy_mis_cont)
#'
#' # View the first few rows of selected data sets
#' head(mantar_dummy_full_cont)
#' head(mantar_dummy_mis_cont)
#'
NULL
