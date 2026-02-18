# Small synthetic test data for WGCNAplus tests
# Creates a small expression matrix with known module structure

create_test_data <- function(ngenes = 100, nsamples = 20, nmodules = 3) {
  set.seed(42)

  # Create module eigengenes
  ME <- matrix(rnorm(nsamples * nmodules), nrow = nsamples, ncol = nmodules)
  colnames(ME) <- paste0("ME", seq_len(nmodules))

  # Create gene expression with module structure
  genes_per_module <- ngenes %/% nmodules
  X <- matrix(0, nrow = ngenes, ncol = nsamples)

  for (m in seq_len(nmodules)) {
    start <- (m - 1) * genes_per_module + 1
    end <- min(m * genes_per_module, ngenes)
    for (g in start:end) {
      X[g, ] <- ME[, m] + rnorm(nsamples, sd = 0.5)
    }
  }

  # Add noise genes for remaining
  remaining <- (nmodules * genes_per_module + 1):ngenes
  if (length(remaining) > 0 && remaining[1] <= ngenes) {
    for (g in remaining) {
      X[g, ] <- rnorm(nsamples)
    }
  }

  rownames(X) <- paste0("gene", seq_len(ngenes))
  colnames(X) <- paste0("sample", seq_len(nsamples))

  # Create samples data frame
  samples <- data.frame(
    group = rep(c("A", "B"), each = nsamples / 2),
    age = rnorm(nsamples, mean = 50, sd = 10),
    row.names = colnames(X)
  )

  # Create contrasts
  contrasts <- matrix(0, nrow = nsamples, ncol = 1)
  contrasts[1:(nsamples / 2), 1] <- 1
  contrasts[(nsamples / 2 + 1):nsamples, 1] <- -1
  rownames(contrasts) <- colnames(X)
  colnames(contrasts) <- "A_vs_B"

  list(
    X = X,
    samples = samples,
    contrasts = contrasts,
    ME = ME
  )
}
