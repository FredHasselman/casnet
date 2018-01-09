.onLoad <- function(libname, pkgname="casnet") {
  op <- options()
  op.casnet <- list(
    casnet.path = system.file(package="casnet"),
    casnet.path_to_rp = system.file("exec", package="casnet"),
    casnet.rp_prefix = "./",
    casnet.rp_command = "rp",
    casnet.sysdel = "rm",
    casnet.syscopy = "cp",
    casnet.install.args = "",
    casnet.name = "Fred Hasselman", #"A toolbox for studying Complex Adaptive Systems and NETworks",
    casnet.desc.author = '"Fred Hasselman <f.hasselman@bsi.ru.nl> [aut, cre]"',
    casnet.desc.license = "GP-L3",
    casnet.desc.suggests = NULL,
    casnet.desc = list()
  )
  toset <- !(names(op.casnet) %in% names(op))
  if(any(toset)) options(op.casnet[toset])
  for(p in c("grDevices","graphics","stats","utils")){requireNamespace(p, quietly = TRUE)}
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
#' @return Message informing whether the procedure was succesful.
#' @export
#'
set_command_line_rp <- function(){

  dl_instruction <- c("Download failed!\nCopying failed!\n\nTo install do the following:\n1. Either go to https://github.com/FredHasselman/casnet/tree/master/inst and download and unzip 'commandline_rp.zip', or go to http://tocsy.pik-potsdam.de/commandline-rp.php \n2. Find the executable for your OS\n3. Copy it to the '/exec' directory under 'casnet' \n4. Rename to 'rp' or 'rp.exe' on Windows\n5. Put an empty text file in '/exec' with the following name 'rp_instal_log.txt' \n 6. Run this code to test if everything is ok: crqa_cl(rnorm(100))")

  os         <- get_os()
  execPath   <- getOption("casnet.path_to_rp")
  sourcePath <- getOption("casnet.path")

  TRYSYS <- FALSE

  ZIP = FALSE
  if(file.exists(system.file("commandline_crp","cl_crp.zip", package="casnet"))){
    utils::unzip(system.file("commandline_crp","cl_crp.zip", package="casnet"),exdir=system.file("commandline_crp",package="casnet"))
    ZIP = TRUE
  }

  # macOSX
  if(os%in%"osx"){

    if(Sys.info()[["machine"]]%in%"x86_64"){
      URL <- "http://tocsy.pik-potsdam.de/RP/rp_osxI"
      sys <- "macOS_x86_64"
      exe <- "rp_osxI.dms"
    } else {
      URL <- "http://tocsy.pik-potsdam.de/RP/rp_osxPPC"
      sys <- "macOS_PowerPC"
      exe <- "rp_osxPPC.dms"
    }
  }

  # Windows
  if(os%in%"windows"){
    options(casnet.rp_prefix="")
    options(casnet.rp_command="rp.exe")
    options(casnet.sysdel="del")
    options(casnet.syscopy="copy")

    sys <- "windows_x86"
    exe <- "rp_x86.exe"
    URL <- "http://tocsy.pik-potsdam.de/RP/rp_x86.exe"
  }

  # Linux
  if(os%in%"linux"){

    if(grepl("amd64",Sys.info()[["machine"]])&grepl("linux-gnu", R.version$os)){
      sys <- "linux_AMD_x86_64gnu"
      exe <- "rp_x86_64.dms"
      URL <- "http://tocsy.pik-potsdam.de/RP/rp_x86_64"
    }

    if(grepl("amd64",Sys.info()[["machine"]])&!grepl("linux-gnu", R.version$os)){
      sys <- "linux_AMD_x86_64i"
      exe <- "rp_x86_64i.dms"
      URL <- "http://tocsy.pik-potsdam.de/RP/rp_x86_64"
    }

    if(grepl("i686",Sys.info()[["machine"]])){
      sys <- "linux_i686"
      exe <- "rp_i686.dms"
      URL <- "http://tocsy.pik-potsdam.de/RP/rp_i686"
    }
  }

  # Get the file from internet
  LOG <- try(utils::download.file(url = URL, destfile = normalizePath(paste0(execPath,"/",getOption("casnet.rp_command")), mustWork = FALSE)))

  if(LOG==0){
    devtools::RCMD(cmd="chmod",options=paste0("a+x ",getOption("casnet.rp_command")),path=normalizePath(execPath, mustWork = FALSE))
    message(paste0("Detected: ",sys,"\n  Copied: ",URL," to ",getOption("casnet.rp_command")," in ",execPath," as the commandline CRP executable"))
    rio::export(data.frame(url=URL),normalizePath(paste0(execPath,"rp_install_log.txt")))
  } else {
    if(ZIP){
      message(paste0("Detected: ",sys, "\n  FAILED to Copy: ",URL," to ",getOption("casnet.rp_command")," in ",execPath," as the commandline CRP executable \n Trying .zip file..."))
      sysCommand <- c(getOption("syscopy"),paste(normalizePath(paste0(sourcePath,"/commandline_rp/",sys,"/",exe)), normalizePath(paste0(execPath,"/",getOption("casnet.rp_command")), mustWork = FALSE)),"chmod",paste("a+x ",normalizePath(paste0(execPath,"/",getOption("casnet.rp_command")), mustWork = FALSE)))
      TRYSY = TRUE
    } else {
      message(dl_instruction)
    }
  }
  if(TRYSYS){
    if(all(nchar(sysCommand)>0)){
      devtools::RCMD(sysCommand[[1]], options = sysCommand[[2]], path = getOption("casnet.path"), quiet = TRUE)
      devtools::RCMD(sysCommand[[3]], options = sysCommand[[4]], path = getOption("casnet.path"), quiet = TRUE)
      rio::export(data.frame(sysCommand=c(paste(sysCommand[[1]],sysCommand[[2]]),paste(sysCommand[[3]],sysCommand[[4]]))),normalizePath(paste0(getOption("casnet.path_to_rp"),"/rp_instal_log.txt"), mustWork = FALSE))
    }
  }
}

