###########################################################################################
#' simOneSpecimen
#' Simulates a set of 1D landmarks for a single specimen for a given number of landmarks 
#' that have a given number of modules. Values are drawn from a normal distribution with 
#' a random mean from 1-50 for each module.
#' distribution 
#' @name simOneSpecimen
#' @param nmodule The number of modules required.
#' @param nlandmark The number of landmarks per module.
#' @param sd The standard deviation of a normal distribution to draw points from
#' @export

simOneSpecimen <- function(nmodule, nlandmark, sd) {
  res <- matrix(ncol = 3, nrow = (nlandmark * nmodule))
  landmarks <- matrix(
    c(seq(from = 1, to = (nlandmark * nmodule), by = nlandmark),
    seq(from = 1, to = (nlandmark * nmodule), by = nlandmark) + nlandmark - 1),
    nrow = nmodule
    )
  means <- runif(nmodule, min = 1, max = 50)

  for (i in 1:nrow(landmarks)) {
    res[landmarks[i, 1]:landmarks[i, 2], ] <- rnorm(length(res[landmarks[i, 1]:landmarks[i, 2], ]),
      means[i], sd = sd)
  }
  rownames(res) <- as.character(1:nrow(res))
  return(abs(res))
}

###########################################################################################
#' simXYZLandmarks
#' Simulates XYZ landmarks for some number of species/specimens. Invokes simOneSpecimen
#' to generat XYZ for each specimen n times. Generates modules of equal size.
#' @name simXYZLandmarks
#' @param nmodule The number of modules required
#' @param nlandmark The number of landmarks per module.
#' @param sd Standard deviation for each module.
#' @export

simXYZLandmarks <- function(nmodule, nlandmark, nspecimen, sd) {
  xyz <- array(dim = c((nlandmark * nmodule), 3, nspecimen),
    dimnames = list(as.character(1:(nlandmark * nmodule)), c("x", "y", "z"), paste0("specimen_", 1:nspecimen)))

  for (i in 1:dim(xyz)[3]) {
    xyz[,,i] <- simOneSpecimen(nmodule, nlandmark, sd)
  }

  return(abs(xyz))
}

###########################################################################################
#' generateModels
#' Generates a set of models to test in EMMLi for a simulated dataset to add to the known, 
#' true, model for testing purposes. Based on modules of equal size, but can be modified
#' post-hoc.
#' @name generateModels
#' @param nmodule Number of modules in the real model.
#' @param nlandmark Number of landmarks in each module in the real model.
#' @param nmods Number of models in the final models table.
#' @param random Logical - if TRUE then a random reshuffling of the true model and a totally random model are added in.
#' @export

generateModels <- function(nmodule, nlandmark, nmods, random = TRUE) {
  # Function to find factors.
  factorise <- function(x) {
    x <- as.integer(x)
    div <- seq_len(abs(x))
    factors <- div[x %% div == 0L]
    factors <- list(neg = -factors, pos = factors)
    return(factors)
  }

  # First generate true model. 
  module_size <- nlandmark
  len <- nmodule * module_size
  models <- matrix(nrow = len, ncol = (nmods + 2))
  models[ , 1] <- c(1:len)
  models[ , 2] <- unlist(lapply(1:nmodule, function(x) rep(x, module_size)))
  
  # Other models
  fact <- factorise(len)$pos
  # exclude the "no modules" model and the true model
  fact <- fact[-length(fact)]
  fact <- fact[!fact == len/nmodule]
  for (i in 1:nmods) {
    idx <- length(fact) - (i - 1)
    models[ , (i + 2)] <- unlist(lapply(1:(len / fact[idx]), function(x) rep(x, fact[idx])))
  }

  models <- models[order(models[ , 1]), ]
  if (random == TRUE) {
    colnames(models) <- c("landmark", "truemod", paste0("badmod_", 1:nmods))
    random_true <- sample(models[ , 2])
    random_complex <- sample(models[ , ncol(models)])
    models <- cbind(models, random_true, random_complex)
  } else {
    colnames(models) <- c("landmark", "truemod", paste0("badmod_", 1:nmods))  
  }
  
  models <- as.data.frame(models)
  models$landmark <- as.factor(models$landmark)
  return(as.data.frame(models))
}

 ######  ##     ## ########   ######     ###    ##     ## ########  ##       #### ##    ##  ######   
##    ## ##     ## ##     ## ##    ##   ## ##   ###   ### ##     ## ##        ##  ###   ## ##    ##  
##       ##     ## ##     ## ##        ##   ##  #### #### ##     ## ##        ##  ####  ## ##        
 ######  ##     ## ########   ######  ##     ## ## ### ## ########  ##        ##  ## ## ## ##   #### 
      ## ##     ## ##     ##       ## ######### ##     ## ##        ##        ##  ##  #### ##    ##  
##    ## ##     ## ##     ## ##    ## ##     ## ##     ## ##        ##        ##  ##   ### ##    ##  
 ######   #######  ########   ######  ##     ## ##     ## ##        ######## #### ##    ##  ######  

###########################################################################################
#' subsampleLandmarks
#' This function randomly removes landmarks from a dataset, and adjusts the corresponding 
#' models object to match the newly subsampled landmark data.
#' @name subsampleLandmarks
#' @param landmarks A dataframe of landmarks
#' @param fraction The decimal fraction to subsample down to. e.g. 0.2 will return 20% of the original landmarks.
#' @param models The models for the landmarks
#' @param min_landmarks The minimum number of landmarks for a module in any model. When subsampling goes below this threshold landmarks are randomly seleclted from the original data and added back in in order to make up to min_landmarks. This means that sometimes the exact fraction desired is not returned.
#' @export

subsampleLandmarks <- function(landmarks, fraction, models, min_landmark) {
  sampleSpecimen <- function(lms, kps) {
    return(lms[kps, ])
  }

  x <- lapply(models[2:length(models)], table)
  x <- sapply(x, length)
  kps <- sort(sample(nrow(landmarks[,,1]), round(nrow(landmarks[,,1]) * fraction, 0)))
  m <- models[ , 2:ncol(models)]
  for (i in 1:ncol(m)) {
    mod_original <- unique(m[ , i])
    ss_mod_count <- table(m[kps, i])

    if (length(ss_mod_count) != length(mod_original)) {
      missing <- mod_original[!mod_original %in% names(ss_mod_count)]
      for (j in missing) {
        kps <- c(kps, sample(which(m[ , i] == j), min_landmark))
      }
      kps <- sort(unique(kps))
    }

    if (any(ss_mod_count < min_landmark)) {
      short <- names(ss_mod_count)[ss_mod_count < min_landmark]
      for (j in short) {
        kps <- c(kps, sample(which(m[ , i] == j), min_landmark))
      }
      kps <- sort(unique(kps))
    }
  }

  res <- array(dim = c(length(kps), 3, dim(landmarks)[3]))
  for (i in 1:dim(landmarks)[3]) {
    res[,,i] <- sampleSpecimen(landmarks[,,i], kps)
  }
  
  new_models <- models[kps, ]
  return(list(landmarks = res, models = new_models, 
    true_subsample = round(dim(res)[1] / dim(landmarks)[1], 2)))
}

###########################################################################################
#' getCorrs
#' This function takes a fitted EMMLi model and the models used to fit it, and then returns
#' all the correaltions for the best fitting model
#' @name getCorrs
#' @param emm A fitted EMMLi model (output of EMMLi)
#' @param models The models object
#' @param corr The original correlation matrix.
#' @export

getCorrs <- function(emm, models, corr) {
  best_mod <- rownames(emm$results)[which(emm$res[ , "dAICc"] == 0)]
  best_mod <- strsplit(best_mod, " ")[[1]][1]
  best_mod <- paste(head(strsplit(best_mod, "\\.")[[1]], n = -2), collapse = ".")

  symmet = corr
  symmet[upper.tri(symmet)] = t(symmet)[upper.tri(symmet)]

  model <- paste0("models$", best_mod)
  lms <- array(eval(parse(text = model)))
  modNF <- stats::na.omit(cbind(1:nrow(models), lms))
  w <- unique(modNF[, 2])
  
  all_modules <- list()
  modules <- list()
  btw_mod = list()
  betweenModules = list()
  withinModules = list()
  unintegrated = list()
  betweenFloat = list()

  for(i in seq(length(w))){
    # identify landmarks within a class
    fg <- modNF[modNF[, 2] == w[i], ] 
    
    # coefficients between identified landmarks.
    l <- corr[fg[, 1], fg[, 1]] 
    modules[[i]] <- (as.array(l[!is.na(l)]))
  }
  names(modules) <- paste("Module", w)
  
  if (length(w) > 1) {
    # make combinations of modules for between module.
    cb <- utils::combn(w, 2) 
    for (i in seq(dim(cb)[2])){
      fg1 <- modNF[modNF[, 2] == cb[1, i], ]
      fg2 <- modNF[modNF[, 2] == cb[2, i], ]
      
      # setdiff(A,B) - present in A but not in B
      btw <- symmet[as.integer(setdiff(fg2[, 1], fg1[, 1])), as.integer(setdiff(fg1[, 1], fg2[, 1]))] 
      btw_mod[[i]] <- btw[!is.na(btw)]
    }
    names(btw_mod) <- paste(cb[1, ], "to", cb[2, ])
    betweenModules['betweenModules'] = list(as.vector(rle(unlist(btw_mod))$values))
    withinModules['withinModules'] = list(as.vector(rle(unlist(modules))$values))
  }

  all_modules[[model]] = c(modules, btw_mod, withinModules, betweenModules)
  
  return(all_modules)
}

###########################################################################################
#' subSampleEMMLi
#' This takes some landmarks, and then will subsample them down to one or more subsampling
#' fractions, calculate the correaltion matrix (using dotcorr only) and then fit EMMLi,
#' returning the result of all the subsampling (if multiple fractions). If a single fraction
#' the nsim argument is required with a number of times to repeat subsampling at that level.
#' @name subSampleEmmli
#' @param landmarks Landmarks dataframe
#' @param fractions Either a single subsampling fraction (in which case nsim is required) or a vector of fractions.
#' @param models The models to test.
#' @param min_landmark The minimum number of landmarks to subsample to
#' @param aic_cut This is the threshold of dAICc below which two models are considered to be not different. When this occurs multiple models can be returned as the best.
#' @param return_corr Logical - if TRUE then the full correlations of within and between modules are returned for the best models.
#' @param summarise Logical - if TRUE then the subsampling results are summarised and retuned.
#' @param nsim If a single subsampling fraction, this is the number of times that subsampling fraction is repeated.
#' @export

subSampleEMMLi <- function(landmarks, fractions, models, min_landmark, aic_cut = 7, 
  return_corr = TRUE, summarise = TRUE, nsim = NULL) {
  # Much like before, except this one takes actual data and actual models rather than 
  # simulation parameters.
  require(paleomorph)
  require(tibble)
  #require(EMMLi)
  
  if (!is.null(nsim)) {
    if (length(fractions) > 1) {
      stop("Only one subsampling fraction is allowed when multiple subsampling simulations are run.")  
    }
    fractions <- rep(fractions, nsim)    
  }

  fitEmmli <- function(landmarks, models, fraction, min_landmark, aic_cut) {

    dat <- subsampleLandmarks(landmarks = landmarks, fraction = fraction, models = models, 
      min_landmark = min_landmark)
    c <- dotcorr(dat$landmarks)

    emm <- EMMLi(corr = c, mod = dat$models, N_sample = dim(landmarks)[3], all_rhos = TRUE)
    emmli <- as.data.frame(emm$results)
    best_models <- emmli[emmli$dAICc <= aic_cut, ]
    best_names <- rownames(best_models)
    best_models <- do.call(rbind, best_models)
    colnames(best_models) <- best_names
    names(emm$all_rhos) <- sapply(names(emm$all_rhos), function(x) strsplit(x, "\\$")[[1]][2])
    names(emm$all_rhos) <- unlist(strsplit(names(emm$all_rhos), "1$"))
    rhos <- emm$all_rhos[names(emm$all_rhos) %in% best_names]
    all_data <- list(landmarks = dat$landmarks, models = dat$models)

    if (return_corr) {
      best_corr <- getCorrs(emm, dat$models, c)
      res <- list(best_models = best_models, rhos_best = rhos, best_corr = best_corr, 
        all_data = all_data, true_subsample = dat$true_subsample) 
    } else {
      res <- list(best_models = best_models, rhos_best = rhos, all_data = all_data,
        true_subsample = dat$true_subsample)
    }
    # print(ncol(best_models))
    # if (ncol(best_models) > 1) {
    #   print("multiple best")
    # }
    return(res)
  }

  res <- lapply(fractions, function(x) fitEmmli(fraction = x, landmarks = landmarks, 
    models = models, min_landmark = min_landmark, aic_cut = aic_cut))
  
  names(res) <- fractions

  if (summarise) {
    summ <- summariseResults(res)
    ret <- list(results = res, summary = summ)
  } else {
    ret <- res
  }
  return(ret)
}

########  ##     ## ##    ## ##        #######          ######## ##     ## ##     ## ##       #### 
##     ## ##     ##  ##  ##  ##       ##     ##         ##       ###   ### ###   ### ##        ##  
##     ## ##     ##   ####   ##       ##     ##         ##       #### #### #### #### ##        ##  
########  #########    ##    ##       ##     ## ####### ######   ## ### ## ## ### ## ##        ##  
##        ##     ##    ##    ##       ##     ##         ##       ##     ## ##     ## ##        ##  
##        ##     ##    ##    ##       ##     ##         ##       ##     ## ##     ## ##        ##  
##        ##     ##    ##    ########  #######          ######## ##     ## ##     ## ######## #### 

###########################################################################################
#' phyloEmmli
#' Takes landmarks and a phylogeny and then corrects the landmarks for the phylogeny according
#' to one of two methods, and then either returns the corrected landmarks, or calculates the 
#' correlation matrix (using dotcorr) and fits EMMLi. Species missing from data or tree are
#' automatically dropped.
#' @name phyloEmmli
#' @param landmarks The landmarks, with species names as rownames.
#' @param phylo A phylogeny describing the relationship between species in landmarks.
#' @param method Either "pgls" or "ic". If PGLS then corrected data are calculated as the residuals of a phylogenetic least squares regression against 1, if IC then independent contrasts.
#' @param EMMLi Logical - if TRUE EMMLi is fit and the results returned.
#' @param ... Extra arguments required for EMMLi (at minimum models, and N_sample)
#' @export

phyloEmmli <- function(landmarks, phylo, method = "pgls", EMMLi = FALSE, ...) {
  require(paleomorph)
  require(geomorph)
  require(caper)
  require(ape)
  # first check that the three is a tree.

  if (class(phylo) != "phylo") {
    stop("Tree must be an object of class 'phylo'.")
  }

  # now check that the data has rownames that are the species.
  if (is.null(rownames(landmarks))) {
    stop("Landmarks must have species names as rownames.")
  }

  # check if mod and N_sample are provided if EMMLi is TRUE.
  if (EMMLi) {
    if (!exists(mod)) {
      stop("mod must be provided to fit EMMLi to phylo landmarks.")
    }
  }

  # now check that the species on the tree match the species 
  sp_lm <- rownames(landmarks)
  sp_tr <- phylo$tip.label

  if (sum(sp_lm %in% sp_tr) != nrow(landmarks)) {
    missing <- sum(!sp_lm %in% sp_tr)
    if (missing == nrow(landmarks)) {
      stop("No species in dataset found on tree.")
    }
    print("Dropping", missing, "species from dataset - not in tree.")
    landmarks <- landmarks[sp_lm %in% sp_tr, ]
  }

  if (sum(sp_tr %in% sp_lm) != length(phylo$tip.label)) {
    missing <- sum(!sp_tr %in% sp_lm)
    if (missing == length(phylo$tip.label)) {
      stop("No species on tree found in dataset.")
    }
    print("Dropping", missing, "species from tree - not in dataset.")
    xtips <- phylo$tip.label[!sp_tr %in% sp_lm]
    phylo <- drop.tip(phylo, xtips)
  }

  # Now, depending on method, replace the real data with phylogenetically
  # corrected data.

  if (method == "pgls") {
    lms <- colnames(landmarks)
    landmarks$names <- rownames(landmarks)
    comp_data <- comparative.data(phylo, as.data.frame(landmarks), names = names)
    x <- lapply(lms, function(x) pgls(formula(paste(x, "~", 1)), comp_data)$phyres)
    phy_landmarks <- do.call(cbind, x)
    rownames(phy_landmarks) <- rownames(landmarks)
  } else if (method == "ic") {
    # match the order of the dataset to the order of the tree.
    landmarks <- landmarks[match(phylo$tip.label, rownames(landmarks)), ]
    phy_landmarks <- apply(landmarks, 2, function(x) pic(x, phylo))
  }

  # Fit or don't fit EMMLi.
  if (!EMMLi) {
    res <- phy_landmarks
  } else if (EMMLi) {
    arr <- arrayspecs(phy_landmarks, ncol(phy_landmarks), 3)
    corr <- dotcorr(arr)
    N_sample <- dim(landmarks)[3]
    emm <- EMMLi(corr = corr, mod = mod, N_sample = N_sample)
    res <- list(EMMLi = emm, phy_landmarks = phy_landmarks)
  }

  return(res)
}

 ######  #### ##     ## ########  ##       #### ########  ######     ###    ######## ####  #######  ##    ## 
##    ##  ##  ###   ### ##     ## ##        ##  ##       ##    ##   ## ##      ##     ##  ##     ## ###   ## 
##        ##  #### #### ##     ## ##        ##  ##       ##        ##   ##     ##     ##  ##     ## ####  ## 
 ######   ##  ## ### ## ########  ##        ##  ######   ##       ##     ##    ##     ##  ##     ## ## ## ## 
      ##  ##  ##     ## ##        ##        ##  ##       ##       #########    ##     ##  ##     ## ##  #### 
##    ##  ##  ##     ## ##        ##        ##  ##       ##    ## ##     ##    ##     ##  ##     ## ##   ### 
 ######  #### ##     ## ##        ######## #### ##        ######  ##     ##    ##    ####  #######  ##    ## 

###########################################################################################
#' simplifyEMMLi
#' Attempts to simplify a fitted EMMLi model, looking for a simpler model that fits the
#' data better by merging modules. Pairs of modules are identified as candidates for merging
#' if their between-module rho is higher than either of the within-module rhos by more than
#' 2 * SD of the modules combined. This repeats until the model does not change or improve.
#' Alternatively, pairs of modules can be offered as candidates for merging - in this 
#' instance there is not exploration and just the suggested pairs are tested. Pairs are 
#' offered as a list, where each element is a vector of two module numbers.
#' @name simplifyEMMLi
#' @param fitted_emmli Fitted EMMLi model output.
#' @param corr_matrix The correlation matrix that EMMLi was fitted to.
#' @param models The original models that EMMLi tested.
#' @param candidates A list of pairs of modules to test merging (each element is a vector of length 2 with module numbers).
#' @param N_sample The sample size for the original EMMLi fit.
#' @param correction if "normal" uses the normal EMMLi calculation for K, if "new" then uses the experimental adjustment to AICc (adding nmodules - 1 to K)

identifyCandidates <- function(fitted_emmli, corr_matrix, models) {
  fit_mods <- fitted_emmli$results
  best_mod <- names(fitted_emmli$rho)
  start_mod <- data.frame(data.point = models$Data.point,
    original_best = models[ , grep(strsplit(best_mod, "\\.")[[1]][1], colnames(models))])
  rhos <- fitted_emmli$rho[[1]]

  # Seperate between from within.
  within_rhos <- rhos[ , grep("^Module", colnames(rhos))]
  between_rhos <- rhos[ , grep("to", colnames(rhos))]

  all_corrs <- getCorrs(fitted_emmli, models, corr_matrix)[[1]] 
  between_rhos <- between_rhos[,order(between_rhos[2,], decreasing = TRUE)]

  # Select pairs that have a between rho greater than either of the withings by
  # a factor of 2*SD of the combined within-correlations of both modules pooled.

  pairs <- colnames(between_rhos)
  candidates <- list()
  for (i in seq_along(pairs)) {
    modules <- strsplit(pairs[1], " to ")[[1]]
    mod_1_rho <- rhos["MaxL_p" , grep(paste("Module", modules[1]), colnames(rhos))]
    mod_2_rho <- rhos["MaxL_p" , grep(paste("Module", modules[2]), colnames(rhos))]
    between_rho <- rhos["MaxL_p", grep(paste(modules[1], "to", modules[2]), colnames(rhos))]
    combined_mods <- c(
      all_corrs[[paste("Module", modules[1])]], 
      all_corrs[[paste("Module", modules[2])]]
      )
    sd_comb <- sd(combined_mods)
    if (between_rho - mod_1_rho > 2 * sd_comb) {
      candidates[[length(candidates) + 1]] <- pairs[i]
    } else if (between_rho - mod_2_rho > 2 * sd_comb) {
      candidates[[length(candidates) + 1]] <- pairs[i]
    }
  }

  candidates <- lapply(candidates, function(x) as.numeric(strsplit(x, " to ")[[1]]))
  # return candidates, and the starting model (the best model from the fitted EMMLi)
  return(list(candidates = candidates, start_mod = start_mod))
}

simplifyEMMLi <- function(fitted_emmli, corr_matrix, models, candidates = NULL,
  correction = "normal", N_sample) {
  # First of all look at the best-fitting EMMLi model, and find pairs
  # of modules that have similar between to within module correlations.

  if (is.null(candidates)) {
    x <- identifyCandidates(fitted_emmli, corr_matrix, models)
    candidates <- x$candidates
    start_mod <- start_mod
    if (length(candidates) == 0) {
      stop("No candidate pairs found.")
    }
  } else {
    # Use user-supplied candidates and identify starting model.
    candidates <- candidates
    best_mod <- names(fitted_emmli$rho)
    start_mod <- data.frame(data.point = models$Data.point,
      original_best = models[ , grep(strsplit(best_mod, "\\.")[[1]][1], colnames(models))])
  }

  repeat {
    # If there are no candidates, break the loop - means the model doesn't simplify any further
    if (length(candidates) == 0) {
      break
    }
    
    new_emms <- list()
    # test all candidate pairs.
    for (i in seq_along(candidates)) {      
      test_mods <- start_mod
      # Combine the two modules in a new_model
      test_mods$new_model <- test_mods$original_best
      test_mods$new_model[test_mods$new_model == candidates[[i]][1]] <- candidates[[i]][2]
      # fit EMMLi
      new_emms[[i]] <- EMMLi(corr = corr_matrix, mod = test_mods, N_sample = N_sample, correction = correction)

      # Now if the best model is "original" go to the next i, otherwise break this loop.
      # if (!grepl("original_best", rownames(emm$results)[emm$results[ , "dAICc"] == 0])) {
      #   print("Better model found.")
      #   models <- test_mods
      #   new_emm <- emm
      #   break
      # }
    
    }
    best_liks <- rep(NA, length(candidates))
    for (i in seq_along(new_emms)) {
      # If the original model isn't the best record the likelihood of it.
      if (!grepl("original_best", 
        rownames(new_emms[[1]]$results)[new_emms[[1]]$results[ , "dAICc"] == 0])) {
        best_liks[i] <- new_emms[[i]]$results[new_emms[[i]]$results[ , "dAICc"] == 0, "MaxL"]
      }
    }

    # If best_liks is all NA then none of the merges helped, and break the loop.
    # Else take the one with the best likelihood and calculate new canidate pairs.
    # If there are no candidates the loop breaks, if there are, it repeats.
    if (all(is.na(best_liks))) {
      break
    } else {
      new_emm <- new_emms[which.max(best_liks)]
      x <- identifyCandidates(fitted_emmli = new_emm, models = test_mods, 
        corr = corr_matrix)  
      candidates <- x$candidates
      start_mod <- x$start_mod
    }
  }
  # If new_emm exists it will be the simplest version of the model, so return it
  # and the test_mods. This will be the most recent test mods.
  # If it doesn't say so, and break.
  if (exists("new_emm")) {
    return (list(emmli = new_emm, model = test_mods))
  } else {
    stop("No simplification possible.")
  }
}

 ######  ##     ## ##     ## ##     ##    ###    ########  ##    ## 
##    ## ##     ## ###   ### ###   ###   ## ##   ##     ##  ##  ##  
##       ##     ## #### #### #### ####  ##   ##  ##     ##   ####   
 ######  ##     ## ## ### ## ## ### ## ##     ## ########     ##    
      ## ##     ## ##     ## ##     ## ######### ##   ##      ##    
##    ## ##     ## ##     ## ##     ## ##     ## ##    ##     ##    
 ######   #######  ##     ## ##     ## ##     ## ##     ##    ##  

###########################################################################################
#' summariseResults
#' This will summarise the output of subSampleEmmli. It returns the rhos of the best models
#' across the set of subsampled analyses. 
#' @name summariseResults
#' @param res The output of subSampleEmmli
#' @export

summariseResults <- function(res) {

  getBestMods <- function(rs) {
    bestMods <- vector(mode = "list", length = length(rs))
    for (i in seq_along(rs)) {
      bestMods[[i]] <- cbind(t(rs[[i]][[1]]), 
        subsample = as.numeric(names(rs)[i]), 
        true_subsample = as.numeric(rs[[i]]$true_subsample))
    }
    bestMods <- do.call(rbind, bestMods)
    bestMods <- cbind(bestMods, subsample = as.numeric(names(rs)))
    return(bestMods)
  }

  getRhos <- function(rs) {
    rho_names <- unlist(lapply(rs, function(x) names(x[[2]])))
    allRhos <- vector(mode = "list", length = length(unique(rho_names)))
    names(allRhos) <- unique(rho_names)
    for (i in seq_along(unique(rho_names))) {
      allRhos[[i]] <- vector(mode = "list", length = sum(rho_names == names(allRhos)[[i]]))
    }

    ref_names <- vector(mode = "list", length = length(allRhos))
    for (i in seq_along(rs)) {
      c_rs <- rs[[i]]$rhos_best
      for (j in seq_along(c_rs)) {
        slt <- which(names(allRhos) == names(c_rs)[j])
        first_null <- which(sapply(allRhos[[slt]], is.null))[[1]]
        xx <- cbind(c_rs[[j]], 
          subsample = as.numeric(names(rs)[i]), 
          true_subsample = as.numeric(rs[[i]]$true_subsample))
        bets <- grep("to", colnames(xx))
        ordered <- sapply(bets, function(x)
          paste(sort(as.numeric(strsplit(colnames(xx)[x], " to ")[[1]])), collapse = " to ")
        )
        # Sort out the right ordering for column names here...
        colnames(xx)[bets] <- ordered
        xx <- xx[ , order(colnames(xx), decreasing = TRUE)]
        xx <- xx[ , c(4:ncol(xx), 1:3)]
        allRhos[[slt]][[first_null]] <- xx
      }
    }

    for (i in seq_along(allRhos)) {
      allRhos[[i]] <- do.call(rbind, allRhos[[i]])
      allRhos[[i]] <- allRhos[[i]][rownames(allRhos[[i]]) != "MaxL", ]
    }
    return(allRhos)    
  }

  bestMods <- getBestMods(res)
  bestRhos <- getRhos(res)
  return(list(bestModels = bestMods, bestRho = bestRhos))
}

###########################################################################################
#' compareModules
#' This takes a correlation matrix, a model definition, and then two module numbers to 
#' compare. It plots a figure of three boxplots - the first two are the correltions within 
#' each of the two modules, and the third is the between-module correlations. These boxes
#' are coloured such that matching colours are not significantly different according to 
#' a Tukey HSD test. The results of the anova and tukey HSD test are also returned.
#' @name compareModules
#' @param corr A correlation matrix
#' @param model Either a vector of numbers describing a model of modules, or a 2 column dataframe with the first bein landmark names and the second being the module definitions.
#' @param test_modules A vector of two module numbers to compare.
#' @param plot Logical - if TRUE the plot is drawn.
#' @export

compareModules <- function(corr, model, test_modules, plot = TRUE) {
  if (plot) {
    require(ggplot2)
    require(multcompView)
  }

  # If landmarks supplied
  if (length(dim(corr)) == 3) {
    require(paleomorph)
    corr <- dotcorr(corr)
  }

  # If model has landmark names
  if (ncol(model) == 2) {
    lms <- array(model[ , 2])
  # If model is just a vector
  } else if (is.vector(model)) {
    lms <- model
  } else if (ncol(model) > 2) {
    stop("Model must either be a vector of model definitions or a data frame with the first column of landmark names and the second of module definitions.")
  }

  symmet = corr
  symmet[upper.tri(symmet)] = t(symmet)[upper.tri(symmet)]

  modNF <- stats::na.omit(cbind(1:nrow(model), lms))
  w <- unique(modNF[, 2])
  w <- w[w %in% test_modules]
  all_modules <- list()
  modules <- list()
  btw_mod = list()
  betweenModules = list()
  withinModules = list()
  unintegrated = list()
  betweenFloat = list()
  for(i in seq(length(w))){
    fg <- modNF[modNF[, 2] == w[i], ] 
    l <- corr[fg[, 1], fg[, 1]] 
    modules[[i]] <- (as.array(l[!is.na(l)]))
  }

  names(modules) <- paste("Module", w)
  if (length(w) > 1) {
    cb <- utils::combn(w, 2) 
    for (i in seq(dim(cb)[2])){
      fg1 <- modNF[modNF[, 2] == cb[1, i], ]
      fg2 <- modNF[modNF[, 2] == cb[2, i], ]
      btw <- symmet[
      as.integer(setdiff(fg2[, 1], fg1[, 1])), 
      as.integer(setdiff(fg1[, 1], fg2[, 1]))] 
      btw_mod[[i]] <- btw[!is.na(btw)]
    }
    names(btw_mod) <- paste(cb[1, ], "to", cb[2, ])
    betweenModules['betweenModules'] = list(as.vector(rle(unlist(btw_mod))$values))
    withinModules['withinModules'] = list(as.vector(rle(unlist(modules))$values))
  }
  all_modules = c(modules, btw_mod, withinModules, betweenModules)

  # Prepare data.
  td <- c(all_modules[[1]], all_modules[[2]], all_modules[[3]])
  group <- as.factor(c(
    rep(names(all_modules)[1], length(all_modules[[1]])), 
    rep(names(all_modules)[2], length(all_modules[[2]])), 
    rep(names(all_modules)[3], length(all_modules[[3]]))
    ))
  td <- data.frame(corrs = td, group = group)
  a <- aov(corrs ~ group, data = td)
  t <- TukeyHSD(a)

  res <- list(anova = a, tukeyhsd = t)

  if (plot) {
    # function to group variables that are not different.
    t_labs <- function(t, v){
      levs <- t[[v]][,4]
      labs <- data.frame(multcompLetters(levs)['Letters'])
      labs$treatment <- rownames(labs)
      labs <- labs[order(labs$treatment) , ]
      return(labs)
    }
    labels <- t_labs(t, "group")

    # add labels to td for colouring.
    td$color <- NA
    for (i in seq_len(nrow(labels))) {
      td$color[td$group == labels$treatment[i]] <- as.character(labels$Letters[i])
    }

    pall <- c("#E69F00", "#56B4E9", "#7BB31A")
    p <- ggplot(td, aes(x = group, y = corrs, fill = color)) +
      geom_boxplot() + 
      scale_fill_manual(values = pall) + 
      xlab("") + 
      ylab("Correlation") + 
      guides(fill = guide_legend(title = "Significance")) + 
      theme(
        legend.position = "none"
      )+
    stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed")
    print(p)
  }
  return(res)  
}

###########################################################################################
#' plotNetwork
#' Plots a network diagram of the output of an EMMLi analysis. Nodes are proportionally
#' sized to within-module rho, and lines are larger and darker according to between-module
#' rho.
#' @name plotNetwork
#' @param rhos The rhos that come out of an EMMLi analysis
#' @param module_names The names of the modules - if absent generic numbers are used. If "rhos" then the nodes are named with the within-module rho
#' @param linecolour The colour of the joining lines.
#' @param title Title for the plot.
#' @param layout Currently not implemented.
#' @export


plotNetwork <- function(rhos, module_names = NULL, linecolour = "#56B4E9", 
  title = NULL, layout = NULL) {
  require(ggplot2)
  require(ggnetwork)
  require(network)
  require(qgraph)

  # get nmodule
  withins <- grep("Module*", colnames(rhos))
  nmodule <- length((grep("Module*", colnames(rhos))))

  # get the rholist
  rholist <- t(rhos)

  # make this object called "words"
  words <- strsplit(rownames(rholist), " ")

  # make a correlation matrix for the network to be plotted off.
  plotcorr <- matrix(data = NA, nrow = nmodule, ncol = nmodule)
  modnums <- unlist(lapply(withins, function(x) strsplit(colnames(rhos)[x], " ")[[1]][2]))

  mods <- sapply(words[withins], function(x) paste(x, collapse = " "))
  mods <- gsub("Module", "M", mods)
  mods <- mods[order(sapply(mods, function(x) 
    as.numeric(strsplit(x, "M ")[[1]][[2]])))]
  colnames(plotcorr) <- rownames(plotcorr) <- mods

  for (i in 1:(length(words) - 1)) {
    if (length(words[[i]]) == 2) {
      module <- paste("M", words[[i]][2])
      plotcorr[module, module] <- rholist[i, "MaxL_p"]
    }

    if (length(words[[i]]) == 3) {
      from_module <- paste("M", (words[[i]][1]))
      to_module <- paste("M", words[[i]][3])
      plotcorr[from_module, to_module] <- rholist[i, "MaxL_p"]
      plotcorr[to_module, from_module] <- rholist[i, "MaxL_p"]
    }
  }

  within <- diag(plotcorr)
  between <- plotcorr

  if (is.null(module_names)) {
    mod.names <- mods
  } else if (module_names == "rhos") {
    mod.names <- within
  }

  if (is.null(layout)) {
    qgraph(between,
      shape = "circle",
      posCol = linecolour,
      labels = mod.names,
      vsize = within * 10,
      diag = FALSE,
      title = title)
  } else {
    qgraph(between,
      shape = "circle",
      posCol = linecolour,
      labels = mod.names,
      vsize = within * 10,
      diag = FALSE,
      title = title,
      layout = layout)    
  }
}

###########################################################################################
#' plotMeanNetwork
#' Plots the mean network inferred from the rhos returned from a set of subsampled analyses.
#' @name plotMeanNetwork
#' @param subsamples The output of subsampleEmmli
#' @param ... Arguments for plotNetwork
#' @export


plotMeanNetwork <- function(subsamples, ...) {
  bestRhos <- subsamples$summary$bestRho
  n <- length(bestRhos)
  par(mfrow = n2mfrow(n))
  mod_names <- names(bestRhos)
  for (i in seq_along(bestRhos)) {
    yy <- colMeans(bestRhos[[i]])
    yy <- yy[1:(length(yy) - 3)]
    yy <- t(data.frame(MaxL_p = yy))
    plotNetwork(yy, title = mod_names[i], ...)
  }
}

###########################################################################################
#' plotRandomSubsamples
#' plots random subsamples from a set of subsampled analyses.
#' @name plotRandomSubsamples
#' @param subsasmples The output of subSampleEmmli
#' @param n The number of random subsamples to plot.
#' @param ... Arguments for plotNetwork
#' @export

plotRandomSubsamples <- function(subsamples, n, ...) {
  samples <- sample(1:length(subsamples$results), n)
  par(mfrow = n2mfrow(n))
  for (i in samples) {
    rh <- subsamples$results[[i]]$rhos_best[[1]]
    plotNetwork(rh, ...)
  }
}


