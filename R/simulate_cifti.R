#' Generate spatial effects over areas for a cifti with cortical surface data
#'
#' @param cifti_obj a \code{xifti} object
#' @param centers_lambda the parameter controlling how many centers the activations will have
#' @param smooth_FWHM The full-width half-maximum smoothing value, in mm
#' @param max_amplitude The maximum value taken by the spatial effect
#' @param n_tasks The number of tasks to generate
#' @param n_subjects The number of subjects for which to generate data
#' @param n_sessions The number of sessions per subject to generate
#' @param n_runs The number of runs per session to generate
#' @param subjects_var The amount of average variance in coefficient location
#'   from subject to subject
#' @param sessions_var The amount of average variance in coefficient location
#'   from session to session within subject
#' @param runs_var The amount of average variance in coefficient location from
#'   run to run within session
#'
#' @importFrom ciftiTools smooth_cifti
#' @importFrom stats rpois
#'
#' @return a list of \code{xifti} objects nested by subject, session, and run
#' @keywords internal
spatial_effects_cifti_cs <- function(cifti_obj,
                                       centers_lambda,
                                       smooth_FWHM,
                                       max_amplitude,
                                       n_tasks,
                                       n_subjects = 1,
                                       n_sessions = 1,
                                       n_runs = 1,
                                       subjects_var = 4,
                                       sessions_var = 2,
                                       runs_var = 1) {
  if(!"xifti" %in% class(cifti_obj)) stop("The cifti_obj must have class 'xifti'.")
  hems <- c('left','right')
  n_vox <- sapply(hems, function(hem) nrow(cifti_obj$data[[paste0("cortex_",hem)]]), simplify = F)
  if(length(max_amplitude) == 1) max_amplitude <- rep(max_amplitude,2)

  all_ciftis <- sapply(paste("Task",seq(n_tasks)), function(h) {
    h_num <- as.numeric(sub("Task ","",h))
    if(length(max_amplitude) == 1) h_num <- 1
    # Simulate overall centers for the tasks
    hem_centers <-
      sapply(n_vox, function(hem_vox) {
        if(!is.null(hem_vox))
          return(sample(x = 1:hem_vox, size = max(rpois(1, centers_lambda), 1)
          ))
        if(is.null(hem_vox)) return(NULL)
      }, simplify = F)
    sapply(paste("Subject",seq(n_subjects)), function(i) {
      # Now jitter for subjects
      hem_centers_i <-
        sapply(hem_centers, function(hc) {
          if(!is.null(hc)) return(hc + rpois(1, subjects_var))
          if(is.null(hc)) return(NULL)
        }, simplify = F)
      sapply(paste("Session",seq(n_sessions)), function(j) {
        # Jitter a little less for sessions
        hem_centers_ij <-
          sapply(hem_centers_i, function(hc) {
            if(!is.null(hc)) return(hc + rpois(1, subjects_var))
            if(is.null(hc)) return(NULL)
          }, simplify = F)
        sapply(paste("Run", seq(n_runs)), function(k) {
          # And jitter just a little bit between runs
          binary_act_ijk <-
            mapply(function(hc,nvox) {
              if(is.null(hc)) return(NULL)
              hem_center_ijk <- hc + rpois(1, subjects_var)
              hem_binary <- rep(0,nvox)
              hem_binary[hem_center_ijk] <- 1
              return(as.matrix(hem_binary))
            }, hc = hem_centers_ij, nvox = n_vox, SIMPLIFY = F)
          cifti_out <- cifti_obj
          for(hem_num in 1:2) {
            if(!is.null(n_vox[[hem_num]]))
              cifti_out$data[[hem_num]] <- binary_act_ijk[[hem_num]]
          }
          # Smooth out the signal
          smooth_cifti <- ciftiTools::smooth_cifti(cifti_out, surf_FWHM = smooth_FWHM)
          # Make sure the amplitude matches the maximum amplitude as input
          for(hem_num in 1:2) {
            if(!is.null(n_vox[[hem_num]])){
              smooth_cifti$data[[hem_num]] <-
                apply(smooth_cifti$data[[hem_num]], 2, function(x)
                  max_amplitude[hem_num] * x / max(x))
            }
          }
          return(smooth_cifti)
        }, simplify = FALSE)
      }, simplify = FALSE)
    }, simplify = FALSE)
  }, simplify = FALSE)
  # Bring all of the tasks together into one cifti for each subject-session-run
  coef_ciftis <- Reduce(function(x,y) {
    mapply(function(xx,yy) {
      mapply(function(xxx,yyy) {
        mapply(function(xxxx,yyyy) {
          out <- xxxx
          if(!is.null(n_vox[[1]])) {
            out$data$cortex_left <- cbind(xxxx$data$cortex_left, yyyy$data$cortex_left)
          }
          if(!is.null(n_vox[[2]])) {
            out$data$cortex_right <- cbind(xxxx$data$cortex_right, yyyy$data$cortex_right)
          }
          return(out)
        }, xxxx = xxx, yyyy = yyy, SIMPLIFY = FALSE)
      }, xxx = xx, yyy = yy, SIMPLIFY = FALSE)
    },xx = x, yy = y, SIMPLIFY = FALSE)
  }, all_ciftis)
  return(coef_ciftis)
}

#' Generate spatial effects over areas for a cifti with subcortical data
#'
#' @param cifti_obj a \code{xifti} object
#' @param centers_lambda the parameter controlling how many centers the activations will have
#' @param smooth_FWHM The full-width half-maximum smoothing value, in mm
#' @param max_amplitude The maximum value taken by the spatial effect
#' @param n_tasks The number of tasks to generate
#' @param n_subjects The number of subjects for which to generate data
#' @param n_sessions The number of sessions per subject to generate
#' @param n_runs The number of runs per session to generate
#' @param subjects_var The amount of average variance in coefficient location
#'   from subject to subject
#' @param sessions_var The amount of average variance in coefficient location
#'   from session to session within subject
#' @param runs_var The amount of average variance in coefficient location from
#'   run to run within session
#'
#' @importFrom ciftiTools smooth_cifti
#' @importFrom stats rpois
#'
#' @return a list of \code{xifti} objects nested by subject, session, and run
#' @keywords internal
spatial_effects_cifti_sub <- function(cifti_obj,
                                     centers_lambda,
                                     smooth_FWHM,
                                     max_amplitude,
                                     n_tasks,
                                     n_subjects = 1,
                                     n_sessions = 1,
                                     n_runs = 1,
                                     subjects_var = 4,
                                     sessions_var = 2,
                                     runs_var = 1) {
  if(!"xifti" %in% class(cifti_obj)) stop("The cifti_obj must have class 'xifti'.")
  n_vox <- nrow(cifti_obj$data$subcort)

  all_ciftis <- sapply(paste("Task",seq(n_tasks)), function(h) {
    h_num <- as.numeric(sub("Task ","",h))
    if(length(max_amplitude) == 1) h_num <- 1
    # Simulate overall centers for the tasks
    hem_centers <-
      sapply(n_vox, function(hem_vox) {
        if(!is.null(hem_vox))
          return(sample(x = 1:hem_vox, size = max(rpois(1, centers_lambda), 1)
          ))
        if(is.null(hem_vox)) return(NULL)
      }, simplify = F)
    sapply(paste("Subject",seq(n_subjects)), function(i) {
      # Now jitter for subjects
      hem_centers_i <-
        sapply(hem_centers, function(hc) {
          if(!is.null(hc)) return(hc + rpois(1, subjects_var))
          if(is.null(hc)) return(NULL)
        }, simplify = F)
      sapply(paste("Session",seq(n_sessions)), function(j) {
        # Jitter a little less for sessions
        hem_centers_ij <-
          sapply(hem_centers_i, function(hc) {
            if(!is.null(hc)) return(hc + rpois(1, subjects_var))
            if(is.null(hc)) return(NULL)
          }, simplify = F)
        sapply(paste("Run", seq(n_runs)), function(k) {
          # And jitter just a little bit between runs
          binary_act_ijk <-
            mapply(function(hc,nvox) {
              if(is.null(hc)) return(NULL)
              hem_center_ijk <- hc + rpois(1, subjects_var)
              hem_binary <- rep(0,nvox)
              hem_binary[hem_center_ijk] <- 1
              return(as.matrix(hem_binary))
            }, hc = hem_centers_ij, nvox = n_vox, SIMPLIFY = F)
          cifti_out <- cifti_obj
          cifti_out$data$subcort <- binary_act_ijk[[1]]
          # Smooth out the signal
          smooth_cifti <- ciftiTools::smooth_cifti(cifti_out, surf_FWHM = smooth_FWHM)
          # Make sure the amplitude matches the maximum amplitude as input
          smooth_cifti$data$subcort <- apply(smooth_cifti$data$subcort,2,
                                             function(x) {
                                               max_amplitude * x / max(x)
                                             })
          return(smooth_cifti)
        }, simplify = FALSE)
      }, simplify = FALSE)
    }, simplify = FALSE)
  }, simplify = FALSE)
  # Bring all of the tasks together into one cifti for each subject-session-run
  coef_ciftis <- Reduce(function(x,y) {
    mapply(function(xx,yy) {
      mapply(function(xxx,yyy) {
        mapply(function(xxxx,yyyy) {
          out <- xxxx
          out$data$subcort <- cbind(xxxx$data$subcort, yyyy$data$subcort)
          # This next line is necessary in order to make a valid subcortical
          # xifti object.
          out$meta$cifti$names <- as.character(seq(ncol(out$data$subcort)))
          return(out)
        }, xxxx = xxx, yyyy = yyy, SIMPLIFY = FALSE)
      }, xxx = xx, yyy = yy, SIMPLIFY = FALSE)
    },xx = x, yy = y, SIMPLIFY = FALSE)
  }, all_ciftis)
  return(coef_ciftis)
}

#' Simulate cifti data for task fMRI for multiple subjects, sessions, and runs
#'
#' @param wb_path Path to the connectome workbench (required)
#' @param brainstructure Which structure(s) should be generated? One of "left",
#'   "right", "subcortical", or "both". The "both" option will return cortical surface simulations for the left hemisphere. It is currently not possible to simulate cortical surface and subcortical data in the same simuulation. Default is "both".
#' @param n_subjects The number of subjects for which data should be generated
#' @param n_sessions The number of sessions of data per subject to be generated
#' @param n_runs The number of runs per session to be generated
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
#' @param subject_var The amount of variance in the location of activations
#'   between subjects.
#' @param session_var The amount of variance in the location of activations
#'   between sessions.
#' @param run_var The amount of variance in the location of activations
#'   between runs.
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
#' @importFrom ciftiTools ciftiTools.setOption ciftiTools.files read_cifti merge_xifti resample_cifti remove_xifti
#' @importFrom stats arima.sim rnorm
#'
#' @examples
#' \dontrun{
#' simulate_cifti()
#' }
simulate_cifti_multiple <-
  function(wb_path,
           brainstructure = "both",
           n_subjects = 1,
           n_sessions = 1,
           n_runs = 1,
           ntasks = 2,
           ntime = 300,
           resamp_res = NULL,
           max_amplitude = 1,
           onsets = NULL,
           durations = NULL,
           TR = 1,
           subject_var = NULL,
           session_var = NULL,
           run_var = NULL,
           ar_error = NULL,
           surfL = NULL,
           surfR = NULL) {
    # This is necessary for the cifti functions
    ciftiTools::ciftiTools.setOption('wb_path',wb_path)
    # Check brainstructure entry
    brainstructure <- match.arg(brainstructure,c("left","right","subcortical", "both"))
    do_left <- brainstructure %in% c('left','both')
    do_right <- brainstructure %in% c('right','both')
    do_sub <- brainstructure == "subcortical"
    # Checks on the inputs
    if(!is.null(onsets)) {
      ntasks = length(onsets)
      if(ntime < max(unlist(onsets))) ntime <- max(unlist(onsets))
    }
    if((is.null(onsets) + is.null(durations)) == 1) stop("Only onsets or durations were specified, but not both. Please give both onsets or durations, or give neither.")
    if(do_sub & !is.null(resamp_res)) message("Subcortical data simulation is not compatible with resampling. No resampling will be done.\n")
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
    cifti_files <- ciftiTools::ciftiTools.files()
    if(do_left | do_right) {
      if(brainstructure == "both") brainstructure <- c("left","right")
      template_cifti <-
        ciftiTools::read_cifti(
          cifti_fname = cifti_files$cifti[[1]],
          surfL_fname = cifti_files$surf[[1]],
          surfR_fname = cifti_files$surf[[2]],
          brainstructures = brainstructure,
          resamp_res = resamp_res
        )
      spatial_effects_cifti <- spatial_effects_cifti_cs
    }
    if(do_sub) {
      template_cifti <-
        ciftiTools::read_cifti(
          cifti_fname = cifti_files$cifti[[4]],
          surfL_fname = NULL,
          surfR_fname = NULL,
          brainstructures = brainstructure,
          resamp_res = NULL
        )
      spatial_effects_cifti <- spatial_effects_cifti_sub
    }
    true_coef_cifti <-
      spatial_effects_cifti(
        cifti_obj = template_cifti,
        # How many distinct areas of activations is
        # simulated from a Poisson(centers_lambda)
        centers_lambda = 0.2,
        smooth_FWHM = 20, # in mm
        max_amplitude = max_amplitude,
        n_tasks = ntasks,
        n_subjects = n_subjects,
        n_sessions = n_sessions,
        n_runs = n_runs
      )
    # Make ciftis for the AR coefficients
    if(is.null(ar_error)) ar_error <- 0
    ar_coefs <- spatial_effects_cifti(
      cifti_obj = template_cifti,
      centers_lambda = 0.2,
      smooth_FWHM = 100,
      max_amplitude = ar_error,
      n_tasks = length(ar_error),
      n_subjects = n_subjects,
      n_sessions = n_sessions,
      n_runs = n_runs
    )
    # Combine information from the AR coefficients and the coefficients
    # to produce the final simulated data
    final_output <- mapply(function(coef_i,ar_i) {
      mapply(function(coef_ij,ar_ij) {
        mapply(function(coef_ijk, ar_ijk) {
          cifti_error <- ar_ijk
          if(do_left) {
            cifti_error$data$cortex_left <-
              t(apply(ar_ijk$data$cortex_left, 1, function(cl_v)
                if(cl_v == 0) {
                  return(rnorm(ntime))
                } else {
                  return(arima.sim(model = list(ar = cl_v), n = ntime))
                }))
          }
          if(do_right) {
            cifti_error$data$cortex_right <-
              t(apply(ar_ijk$data$cortex_right, 1, function(cl_v)
                if(cl_v == 0) {
                  return(rnorm(ntime))
                } else {
                  return(arima.sim(model = list(ar = cl_v), n = ntime))
                }))
          }
          if(do_sub) {
            cifti_error$data$subcort <-
              t(apply(ar_ijk$data$subcort, 1, function(cl_v)
                if(cl_v == 0) {
                  return(rnorm(ntime))
                } else {
                  return(arima.sim(model = list(ar = cl_v), n = ntime))
                }))
          }
          final_cifti <- cifti_error
          if(do_left) {
            final_cifti$data$cortex_left <- final_cifti$data$cortex_left +
              tcrossprod(coef_ijk$data$cortex_left,design)
          }
          if(do_right) {
            final_cifti$data$cortex_right <- final_cifti$data$cortex_right +
              tcrossprod(coef_ijk$data$cortex_right,design)
          }
          if(do_sub) {
            final_cifti$data$subcort <- final_cifti$data$subcort +
              tcrossprod(coef_ijk$data$subcort,design)
            # This next line is necessary in order to make a valid subcortical
            # xifti object.
            final_cifti$meta$cifti$names <- as.character(seq(ntime))
          }
          final_cifti <- final_cifti + 250 # This is done so that preprocessing does
                                            # not artifically inflate values in locations
                                            # with means close to zero.
          return(final_cifti)
        }, coef_ijk = coef_ij, ar_ijk = ar_ij, SIMPLIFY = FALSE)
      }, coef_ij = coef_i, ar_ij = ar_i, SIMPLIFY = FALSE)
    }, coef_i = true_coef_cifti, ar_i = ar_coefs, SIMPLIFY = FALSE)
    all_data_combos <- expand.grid(seq(n_subjects),seq(n_sessions),seq(n_runs))
    final_return <- apply(all_data_combos,1,function(x) {
      final_output[[x[1]]][[x[2]]][[x[3]]]
    })
    coef_return <- apply(all_data_combos,1,function(x) {
      true_coef_cifti[[x[1]]][[x[2]]][[x[3]]]
    })
    ar_return <- apply(all_data_combos,1,function(x) {
      ar_coefs[[x[1]]][[x[2]]][[x[3]]]
    })
    cifti_names <- apply(all_data_combos,1,function(x) paste0("subj",x[1],"_sess",x[2],"_run",x[3]))
    names(final_return) <- names(coef_return) <- names(ar_return) <-  cifti_names
    output <- list(
      simulated_cifti = final_return,
      coef_cifti = coef_return,
      ar_cifti = ar_return,
      design = design
    )
    class(output) <- "simCifti"
    return(output)
  }
