"Para.HOF.list" <-
    function (resp, ...) 
{
    out <- lapply(resp, Para, ...)
#    class(out) <- "Para.HOF.frame"
    out
}
