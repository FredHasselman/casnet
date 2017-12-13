.onLoad <- function(libname, pkgname="casnet") {
  op <- options()
  op.casnet <- list(
    casnet.path = find.package("casnet"),
    casnet.path_to_rp = normalizePath(paste0(find.package("casnet"),"/exec")),
    casnet.rp_prefix = "./",
    casnet.install.args = "",
    casnet.name = "Fred Hasselman", #"A toolbox for studying Complex Adaptive Systems and NETworks",
    casnet.desc.author = '"Fred Hasselman <f.hasselman@bsi.ru.nl> [aut, cre]"',
    casnet.desc.license = "GP-L3",
    casnet.desc.suggests = NULL,
    casnet.desc = list()
  )
  toset <- !(names(op.casnet) %in% names(op))
  if(any(toset)) options(op.casnet[toset])
  invisible()
}


#' Which OS is running?
#'
#' @return A string, "osx", "windows", "linux"
#' @export
#'
get_os <- function(){
  sysinf <- Sys.info()
  os <- NA
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')  os <- "osx"
    if (os == 'Windows') os <- "windows"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}

#' Set command line RQA executable
#'
#' @return nothing
#' @export
#'
set_command_line_rp <- function(){
  os <- get_os()
  execPath <- getOption("casnet.path_to_rp")
  sourcePath <- getOption("casnet.path")

  # macOSX
  if(os%in%"osx"){
    if(Sys.info()[["machine"]]%in%"x86_64"){
      sysCommand <- c("cp",paste(normalizePath(paste0(sourcePath,"/inst/commandline_rp/macOS_x86_64/rp_osxI.dms")), normalizePath(paste0(execPath,"/rp"), mustWork = FALSE)), "chmod",paste("a+x",normalizePath(paste0(execPath,"/rp"), mustWork = FALSE)))
      message("Detected: macOS x86_64 \n  Copied: 'rp_osxI.dms' to 'rp' in package subdirectory 'exec' as the commandline CRP executable")
    } else {
      sysCommand <- c("cp",paste(normalizePath(paste0(sourcePath,"/commandline_rp/macOS_PowerPC/rp_osxPPC.dms")), normalizePath(paste0(execPath,"/rp"), mustWork = FALSE)),"chmod",paste("a+x ",normalizePath(paste0(execPath,"/rp"), mustWork = FALSE)))
      message(paste0("Detected: macOS PowerPC \n  Copied: 'rp_osxPPC.dms' to 'rp' in package subdirectory 'exec' as the commandline CRP executable"))
    }
  }

  # Windows
  if(os%in%"windows"){
    sysCommand <- c("cp",paste(normalizePath(paste0(sourcePath,"/commandline_rp/windows_x86/rp_x86.exe")), normalizePath(paste0(execPath,"/rp.exe"), mustWork = FALSE)),"attrib", paste("+s", normalizePath(paste0(execPath,"/rp.exe"), mustWork = FALSE)))
    message(paste0("Detected: Windows \n  Copied: 'rp_x86.exe' to 'rp.exe' in package subdirectory 'exec' as the commandline CRP executable"))
    options(casnet.rp_prefix="/")
  }

  # Linux
  if(os%in%"linux"){
    if(grepl("amd64",Sys.info()[["machine"]])&grepl("linux-gnu", R.version$os)){
      sysCommand <- c("cp",paste(normalizePath(paste0(sourcePath,"/inst/commandline_rp/linux_AMD_x86_64gnu/rp_x86_64.dms")), normalizePath(paste0(execPath,"/rp"), mustWork = FALSE)),"chmod", paste("a+x ",normalizePath(paste0(execPath,"/rp"), mustWork = FALSE)))
      message(paste0("Detected: Linux AMD gnu \n  Copied: 'rp_x86_64.dms' to 'rp' in package subdirectory 'exec' as the commandline CRP executable"))
    }
    if(grepl("amd64",Sys.info()[["machine"]])&!grepl("linux-gnu", R.version$os)){
      sysCommand <- c("cp",paste(normalizePath(paste0(sourcePath,"/inst/commandline_rp/linux_AMD_x86_64i/rp_x86_64i.dms")), normalizePath(paste0(execPath,"/rp"), mustWork = FALSE)),"chmod",paste("a+x",normalizePath(paste0(execPath,"/rp"), mustWork = FALSE)))
      message(paste0("Detected: Linux AMD gnu \n  Copied: 'rp_x86_64i.dms' to 'rp' in package subdirectory 'exec' as the commandline CRP executable"))
    }
    if(grepl("i686",Sys.info()[["machine"]])){
      sysCommand <- c("cp",paste(normalizePath(paste0(sourcePath,"/inst/commandline_rp/linux_i686/rp_i686.dms")), normalizePath(paste0(execPath,"/rp"), mustWork = FALSE)),"chmod",paste("a+x ",normalizePath(paste0(execPath,"/rp"), mustWork = FALSE)))
      message(paste0("Detected: Linux \n  Copied: 'rp_i686.dms' to 'rp' in package subdirectory 'exec' as the commandline CRP executable"))
    }
  }

  if(all(nchar(sysCommand)>0)){
    devtools::RCMD(sysCommand[[1]], options = sysCommand[[2]], path = getOption("casnet.path"), quiet = TRUE)
    devtools::RCMD(sysCommand[[3]], options = sysCommand[[4]], path = getOption("casnet.path"), quiet = TRUE)
    rio::export(data.frame(sysCommand=c(paste(sysCommand[[1]],sysCommand[[2]]),paste(sysCommand[[3]],sysCommand[[4]]))),normalizePath(paste0(getOption("casnet.path_to_rp"),"/rp_instal_log.txt"), mustWork = FALSE))
  }

}


