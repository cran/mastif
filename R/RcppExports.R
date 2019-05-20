# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

kernYrRcpp <- function(dmat, fec, years, seedyear, treeyear, seedrow, treecol) {
    .Call(`_mastif_kernYrRcpp`, dmat, fec, years, seedyear, treeyear, seedrow, treecol)
}

byRcpp <- function(nr, frommat, totmat, summat, minmat, maxmat) {
    .Call(`_mastif_byRcpp`, nr, frommat, totmat, summat, minmat, maxmat)
}

tnormRcpp <- function(lo, hi, mu, sig) {
    .Call(`_mastif_tnormRcpp`, lo, hi, mu, sig)
}

trMVNmatrixRcpp <- function(avec, muvec, smat, lo, hi, whichSample, idxALL) {
    .Call(`_mastif_trMVNmatrixRcpp`, avec, muvec, smat, lo, hi, whichSample, idxALL)
}

betaRcpp <- function(n, X, y, sigma, AI) {
    .Call(`_mastif_betaRcpp`, n, X, y, sigma, AI)
}

randEffectRcpp <- function(gindex, groups, X, y, sigma, AI) {
    .Call(`_mastif_randEffectRcpp`, gindex, groups, X, y, sigma, AI)
}

solveRcpp <- function(A) {
    .Call(`_mastif_solveRcpp`, A)
}

rmvnormRcpp <- function(n, mu, sigma) {
    .Call(`_mastif_rmvnormRcpp`, n, mu, sigma)
}
