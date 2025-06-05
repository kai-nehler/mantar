calculate_sample_size <- function(data, n_calc = c("individual", "average", "max", "total")){

  n_calc <- tolower(n_calc[1])
  n_calc <- match.arg(n_calc, c("individual", "average", "max", "total"))

  if (!anyNA(data)) {
    ns <- rep(nrow(data), ncol(data))
    message("No missing values in data. Sample size for each variable is equal to the number of rows in the data.")
  } else {
    n_calc <- match.arg(n_calc)
    if (n_calc == "average"){
      ns <- rep(mean(colSums(!is.na(data))), ncol(data))
      } else if (n_calc == "individual"){
        ns <- colSums(!is.na(data))
        } else if (n_calc == "max"){
          ns <- rep(max(colSums(!is.na(data))), ncol(data))
          } else if (n_calc == "total"){
            ns <- rep(nrow(data), ncol(data))
            }
  }
  return(ns)
}




reg_ic_calc <- function(resid_var, n, n_preds, k){

  k <- parse(text = k) |> eval()

  # calculate log-likelihood
  LL <- -(n/2) * (log(resid_var) + log(2*pi) + 1)

  # calculate information criterion; add 1 to the for the intercept
  IC <- -2 * LL + k * (n_preds + 1)

  return(IC)
}

