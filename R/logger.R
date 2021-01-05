#' Log the message to file
#'
#' @param fn String, file name
#' @param logger String, logger name
#' @param newfile Boolean, new logging file
#'
#' @import futile.logger
#' @export
log_to_file <- function(fn, logger = "ROOT", newfile = False) {
    if (newfile) {
        if (file.exists(fn))
            # Delete file if it exists
            file.remove(fn)
    }
    tmp = futile.logger::flog.appender(futile.logger::appender.file(fn), name = logger)
}


#' Log the message to console
#'
#' @param logger String, logger name
#'
#' @import futile.logger
#' @export
log_to_console <- function(logger = "ROOT") {
    tmp = futile.logger::flog.appender(futile.logger::appender.console(), name = logger)
}


#' Log the message to both file and console
#'
#' @param fn String, file name
#' @param logger String, logger name
#' @param newfile Boolean, new logging file
#'
#' @import futile.logger
#' @export
log_to_both <- function(fn, logger = "ROOT", newfile = F) {
    if (newfile) {
        if (file.exists(fn))
            # Delete file if it exists
            file.remove(fn)
    }
    tmp = futile.logger::flog.appender(futile.logger::appender.tee(fn), name = logger)
}


#' Set the log threshold
#'
#' @param level log streshold
#' @param logger logger name
#'
#' @import futile.logger
#' @export
set_log_level <- function(level = "INFO", logger = "ROOT") {
    tmp = futile.logger::flog.threshold(level, name = logger)
}


#' Log message
#'
#' @param level log streshold
#' @param l... log messages
#'
#' @import futile.logger
#' @export
logmsg <- function(..., level = "info") {
    msg = paste0(...)
    if (level == "debug") {
        futile.logger::flog.debug(msg)
    } else if (level == "info") {
        futile.logger::flog.info(msg)
    } else if (level == "warn") {
        futile.logger::flog.warn(msg)
    } else if (level == "error") {
        futile.logger::flog.error(msg)
    } else {
        stop("Log level not found!")
    }
}


