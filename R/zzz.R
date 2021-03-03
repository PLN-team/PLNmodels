.onAttach <- function(libname, pkgname) {
    version <- read.dcf(file=system.file("DESCRIPTION", package=pkgname), fields="Version")
    packageStartupMessage(paste(pkgname, version))
    packageStartupMessage('PLNPCA, PLNmixture and PLNnetwork$stability_selection integrates features from the future package.')
    packageStartupMessage('Set up your plan with future::plan(multicore) or future::plan(multisession) for more speed.')
}
