.onAttach <- function(lib, pkg)  {
    packageStartupMessage("This is eHOF ",
    utils::packageDescription("eHOF", field="Version"),
    appendLF = TRUE)
}
