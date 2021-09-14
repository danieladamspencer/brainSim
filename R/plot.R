#' S3 method: use \code{\link[ciftiTools]{view_xifti_surface}} to plot a \code{"simCifti"} object
#'
#' This function also plots the design matrix, if specified
#'
#' @param x An object of class \code{simCifti}
#' @param type One of \code{coef}, \code{ar}, \code{simulated}, \code{design}.
#'   For \code{coef}, the simulated coefficients are plotted. For \code{ar}, the
#'    autoregressive errors are plotted. \code{simulated} plots the simulated
#'    BOLD time series, and \code{design} will plot the design matrix
#'    coefficients.
#' @param subject If "simulated" is specified for \code{type}, which subject's
#'   BOLD time series should be plotted?
#' @param session If "simulated" is specified for \code{type}, which session's
#'   BOLD time series should be plotted?
#' @param run If "simulated" is specified for \code{type}, which run's
#'   BOLD time series should be plotted?
#' @param zlim Overrides the \code{zlim} argument for
#'  \code{\link[ciftiTools]{view_xifti_surface}}. Default: \code{c(-1, 1)}.
#' @param ...  Additional arguments to \code{\link[ciftiTools]{view_xifti_surface}}
#'
#' @method plot simCifti
#'
# @importFrom ciftiTools view_xifti_surface
#' @importFrom graphics legend lines par
#'
#' @export
#'
plot.simCifti <- function(x, type = c("coef","ar","simulated","design"), subject = NULL, session = NULL, run = NULL, zlim = c(-1,1), ...) {
  # Make sure ciftiTools is installed
  if (!requireNamespace("ciftiTools", quietly = TRUE)) {
    stop("This function requires the `ciftiTools` package. Please install it.")
  }
  # Check the type of plot
  type = match.arg(type, c("coef","ar","simulated","design"))
  if(type %in% c("coef","ar","simulated")) type <- paste0(type,"_cifti")
  # Plot the object
  # >> For coef, ar, or simulated
  if(type %in% c("coef","ar","simulated")) {
    # Set subject, session, and run
    if(is.null(subject)) subject <- 1; if(is.null(session)) session <- 1
    if(is.null(run)) run <- 1
    sub_sess_run <- paste0("subj",subject,"_sess",session,"_run",run)
    ciftiTools::view_xifti_surface(xifti = x[[type]][[sub_sess_run]],zlim = zlim)
  }
  # >> For design
  if(type == "design") {
    n_tasks <- ncol(x$design)
    n_time <- nrow(x$design)
    y_lims <- round(c(min(x$design),max(x$design)),2)
    if(n_tasks == 1) {
      plot(x$design[,1], type = 'l', ylim = y_lims, ylab = "HRF Convolved Task", xlab = "Time")
    }
    if(n_tasks > 1) {
      par(mar=c(5.1, 4.1, 1.1, 8.1), xpd=TRUE)
      plot(x$design[,1], type = 'l', ylim = y_lims, ylab = "HRF Convolved Task", xlab = "Time")
      for(k in seq(2,n_tasks)){
        lines(x$design[,k], lty = k)
      }
      legend('topright', inset = c(-0.2,0),legend = paste("Task",seq(n_tasks)),lty = seq(n_tasks))
    }
  }
}
