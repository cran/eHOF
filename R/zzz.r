.onAttach <- function(lib, pkg)  {
    packageStartupMessage("This is eHOF ",
    utils::packageDescription("eHOF", field="Version"), '; Built: ',
    utils::packageDescription("eHOF", field="Built"),
    paste('\nSee citation(package="eHOF")'),
    appendLF = TRUE)
    options(eHOF.bootselectmessage = TRUE)
    options(repos = c(CRAN="http://cran.r-project.org"))
}
