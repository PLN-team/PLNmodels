.onAttach <- function(libname, pkgname) {
    version <- read.dcf(file=system.file("DESCRIPTION", package=pkgname), fields="Version")
    packageStartupMessage(paste0("This is packages '", pkgname,"' version ",version))
    packageStartupMessage('Use future::plan(multicore/multisession) to speed up PLNPCA/PLNmixture/stability_selection.')
}
