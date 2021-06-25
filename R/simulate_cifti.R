#' Simulate cifti data for task fMRI
#'
#' @param wb_path Path to the connectome workbench (required)
#' @param ntasks The number of tasks with which to simulate response time series
#' @param ntime The length of the response time series
#' @param resamp_res The resolution at which to resample the simuluated cifti
#'   data (default \code{resamp_res = 5000}).
#' @param max_amplitude The maximum amplitude of the coefficients in percent
#'   signal change (default \code{max_amplitude = 1})
#' @param onsets A list of times for the onsets of tasks, each element
#'   representing a task. Overrides \code{ntasks} for determining the number of
#'   tasks if both are given.
#' @param durations A list of durations for each task, each element representing
#'   a task. The list should have the same length as \code{onsets}, and each
#'   list element should have the same length as the respective element within
#'   \code{onsets}.
#' @param TR The repetition time, in seconds. Must be a whole number.
#' @param ar_error A vector of length \code{p} of the AR(p) coefficients for the
#'   error term. The default is \code{NULL}, which results in no autoregressive
#'   error. These coefficient values will be simulated to vary across the
#'   spatial field in order to more accurately reflect most cifti data
#' @param surfL Left surface to use for the simulated cifti. If none are given, the
#'   default surfaces from \code{ciftiTools} are used.
#' @param surfR Right surface to use for the simulated cifti. If none are given, the
#'   default surfaces from \code{ciftiTools} are used.
#'
#' @return A list of objects containing the simulated data and the relevant
#'   information used to create the data. The simulated response data are found
#'   within the \code{cifti} object. The simulated amplitude fields are found
#'   as cifti objects within the \code{coefficients} object. The simulated
#'   design matrix is found in the \code{design} object.
#' @export
#'
#' @importFrom neuRosim specifydesign
#' @importFrom ciftiTools ciftiTools.setOption demo_files read_cifti concat_xifti resample_cifti
#' @importFrom stats arima.sim
#'
#' @examples
#' \dontrun{
#' simulate_cifti()
#' }
simulate_cifti <-
  function(wb_path,
           ntasks = 2,
           ntime = 300,
           resamp_res = NULL,
           max_amplitude = 1,
           onsets = NULL,
           durations = NULL,
           TR = 1,
           ar_error = NULL,
           surfL = NULL,
           surfR = NULL) {
    # This is necessary for the cifti functions
    ciftiTools::ciftiTools.setOption('wb_path',wb_path)
    # Checks on the inputs
    if(!is.null(onsets)) {
      ntasks = length(onsets)
      if(ntime < max(unlist(onsets))) ntime <- max(unlist(onsets))
    }
    if((is.null(onsets) + is.null(durations)) == 1) stop("Only onsets or durations were specified, but not both. Please give both onsets or durations, or give neither.")
    # Make onsets and durations, if they are not supplied
    if(is.null(onsets)) {
      block_period <- round(ntime / 5)
      onsets <- sapply(seq(ntasks), function(task_n) {
        seq((block_period / ntasks)*task_n,
            ntime - (block_period / ntasks)*(ntasks - task_n),
            length.out = 5) - block_period / ntasks
      }, simplify = FALSE)
      durations <- sapply(onsets, function(onset_n) {
        TR
      }, simplify = FALSE)
    }
    # Create the design matrix from the onsets and the durations
    design <-
      neuRosim::specifydesign(
        onsets = onsets,
        durations = durations,
        totaltime = ntime,
        TR = TR,
        effectsize = lapply(1:length(onsets),function(x) 1),
        accuracy = 0.1,
        conv = "double-gamma",
        cond.names = NULL,
        param = NULL
      )
    # Make a cifti for the coefficients
    cifti_files <- ciftiTools::demo_files()
    coef_cifti <-
      ciftiTools::read_cifti(
        cifti_fname = cifti_files$cifti[[1]],
        surfL_fname = cifti_files$surf[[1]],
        surfR_fname = cifti_files$surf[[2]]
      )
    smooth_coef_cifti <- sapply(seq(ntasks), function(n) {
      spatial_effects_cifti(
        cifti_obj = coef_cifti,
        centers_lambda = 0.2,
        smooth_FWHM = 20,
        max_amplitude = max_amplitude
      )
    }, simplify = FALSE)
    smooth_coef_cifti <- ciftiTools::concat_xifti(xifti_list = smooth_coef_cifti)
    # Make ciftis for the AR coefficients
    if(is.null(ar_error)) ar_error <- 0
    ar_coefs <- sapply(ar_error, function(p) {
      spatial_effects_cifti(
        cifti_obj = coef_cifti,
        centers_lambda = 0.2,
        smooth_FWHM = 100,
        max_amplitude = p
      )
    },simplify = F)
    ar_coefs <- ciftiTools::concat_xifti(xifti_list = ar_coefs)
    # Simulate the error term
    cifti_error <- ar_coefs
    cifti_error$data$cortex_left <-
      t(apply(ar_coefs$data$cortex_left, 1, function(cl_v)
        arima.sim(model = list(ar = cl_v), n = ntime)))
    cifti_error$data$cortex_right <-
      t(apply(ar_coefs$data$cortex_right, 1, function(cl_v)
        arima.sim(model = list(ar = cl_v), n = ntime)))
    final_cifti <- cifti_error
    final_cifti$data$cortex_left <- final_cifti$data$cortex_left +
      tcrossprod(smooth_coef_cifti$data$cortex_left,design)
    final_cifti$data$cortex_right <- final_cifti$data$cortex_right +
      tcrossprod(smooth_coef_cifti$data$cortex_right,design)
    final_cifti <- final_cifti + 250
    if(!is.null(resamp_res)) {
      if(resamp_res < 1000) message("Resampling to a resolution below 1000 may oversmooth and deliver undesirable results.")
      final_cifti <- ciftiTools::resample_cifti(x = final_cifti,surfL_original_fname = surfL,surfR_original_fname = surfR,resamp_res = resamp_res)
    }
    return(
      list(
        simulated_cifti = final_cifti,
        coef_cifti = smooth_coef_cifti,
        error_cifti = cifti_error,
        design = design
      )
    )
  }

#' Generate spatial effects over areas for a cifti
#'
#' @param cifti_obj a \code{xifti} object
#' @param centers_lambda the parameter controlling how many centers the activations will have
#' @param smooth_FWHM The full-width half-maximum smoothing value, in mm
#' @param max_amplitude The maximum value taken by the spatial effect
#'
#' @importFrom ciftiTools smooth_cifti
#' @importFrom stats rpois
#'
#' @return a \code{xifti} object
#' @keywords internal
spatial_effects_cifti <- function(cifti_obj, centers_lambda, smooth_FWHM, max_amplitude) {
  if(!"xifti" %in% class(cifti_obj)) stop("The cifti_obj must have class 'xifti'.")
  voxL <- nrow(cifti_obj$data$cortex_left)
  voxR <- nrow(cifti_obj$data$cortex_right)
  left_centers <- sample(x = 1:voxL,size = max(rpois(1,centers_lambda),1))
  right_centers <- sample(x = 1:voxR,size = max(rpois(1,centers_lambda),1))
  left_binary <- rep(0,voxL)
  left_binary[left_centers] <- 1
  right_binary <- rep(0,voxR)
  right_binary[right_centers] <- 1
  cifti_out <- cifti_obj
  cifti_out$data$cortex_left <- as.matrix(left_binary)
  cifti_out$data$cortex_right <- as.matrix(right_binary)
  # cifti_out <- ciftiTools:::fix_xifti(cifti_out) # Don't think this is necessary
  smooth_cifti <- ciftiTools::smooth_cifti(cifti_out, surf_FWHM = smooth_FWHM)
  smooth_cifti$data$cortex_left <- apply(smooth_cifti$data$cortex_left, 2, function(x) max_amplitude * x / max(x))
  smooth_cifti$data$cortex_right <- apply(smooth_cifti$data$cortex_right, 2, function(x) max_amplitude * x / max(x))
  return(smooth_cifti)
}
