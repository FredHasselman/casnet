.onLoad <- function(libname, pkgname="casnet") {
  if(grepl("CNAS.RU.NL",system.file(package="casnet"))){
    cat("\nRU account detected...\n")
    RU <- TRUE
    ppath <- normalizePath(paste0("U:\\",paste0(strsplit(system.file(package="casnet"),"[/]")[[1]][-c(1,2)],collapse="/")), mustWork = TRUE)
    #execpath <- normalizePath(paste0(ppath,"/exec"), mustWork = TRUE)
    #execpath <- normalizePath(paste0("C:\\Users/",paste0(strsplit(system.file(package="casnet"),"[/]")[[1]][2],"/Documents",collapse="/")), mustWork = TRUE)
    execpath <- normalizePath(paste0("C:\\Temp"))
  } else {
    RU <- FALSE
    ppath     <- system.file(package="casnet")
    execpath <- system.file("exec", package="casnet")
  }
  casnet_OS_options <- set_os_options()
  op <- options()
  op.casnet <- list(
    casnet.path = ppath,
    casnet.path_to_rp = execpath,
    casnet.isRU = RU,
    casnet.rp_prefix = casnet_OS_options$rp_prefix,
    casnet.rp_command = casnet_OS_options$rp_command,
    casnet.rp_URL = casnet_OS_options$URL,
    casnet.sysdel = casnet_OS_options$sysdel,
    casnet.syscopy = casnet_OS_options$syscopy,
    casnet.install.args = "",
    casnet.name = "Fred Hasselman", #"A toolbox for studying Complex Adaptive Systems and NETworks",
    casnet.desc.author = '"Fred Hasselman <f.hasselman@bsi.ru.nl> [aut, cre]"',
    casnet.desc.license = "GP-L3",
    casnet.desc.suggests = NULL,
    casnet.desc = list()
  )
  toset <- !(names(op.casnet) %in% names(op))
  if(any(toset)) options(op.casnet[toset])
  for(p in c("grDevices","graphics","stats","utils","ggplot2")){requireNamespace(p, quietly = TRUE)}
  invisible()
}

#' Which OS is running?
#'
#'  Some systems not tested, but based on the cran page: [check flavors](https://cran.r-project.org/web/checks/check_flavors.html)
#'
#' @return A string, "osx", "windows", "linux"
#' @export
#'
get_os <- function(){
  sysinf <- Sys.info()
  os <- NA
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == "Darwin")  os <- "osx"
    if (os == "Windows") os <- "windows"
    if (os == "SunOS")   os <- "sun"
    if (os == "Solaris") os <- "solaris"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    (grepl("linux-gnu", R.version$os))
    os <- "linux-gnu"
  }
  tolower(os)
}

set_os_options <- function(os = get_os()){

  # Default Unix / OSx
  sys        <- "linux"
  exe        <- "rp"
  URL        <- "http://tocsy.pik-potsdam.de/RP/rp_x86_64"
  rp_prefix  <- "./"
  rp_command <- "rp"
  sysdel     <- "rm"
  syscopy    <- "cp"


  # Sun / Solaris
  if(os%in%c("sun","solaris")){
    URL <- "http://tocsy.pik-potsdam.de/RP/rp_sun"
    sys <- "sun_solaris"
    exe <- "rp_sun.dms"
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
    rp_prefix  <- ""
    rp_command <- "rp.exe"
    sysdel     <- "del"
    syscopy    <- "copy"

    sys <- "windows_x86"
    exe <- "rp_x86.exe"
    URL <- "http://tocsy.pik-potsdam.de/RP/rp_x86.exe"
  }

  # Linux
  if(os%in%c("linux","linux-gnu")){

    if(grepl("amd64",Sys.info()[["machine"]])&grepl("linux-gnu", R.version$os)){
      sys <- "linux_AMD_x86_64gnu"
      exe <- "rp_x86_64.dms"
      URL <- "http://tocsy.pik-potsdam.de/RP/rp_x86_64"
    }

    if(grepl("amd64",Sys.info()[["machine"]])&!grepl("linux", R.version$os)){
      sys <- "linux_AMD_x86_64i"
      exe <- "rp_x86_64i.dms"
      URL <- "http://tocsy.pik-potsdam.de/RP/rp_x86_64i"
    }

    if(grepl("i686",Sys.info()[["machine"]])){
      sys <- "linux_i686"
      exe <- "rp_i686.dms"
      URL <- "http://tocsy.pik-potsdam.de/RP/rp_i686"
    }
  }

  return(list(
    sys        = sys,
    exe        = exe,
    URL        = URL,
    rp_prefix  = rp_prefix,
    rp_command = rp_command,
    sysdel     = sysdel,
    syscopy    = syscopy
  ))
}



#' Set command line RQA executable
#'
#' @return Message informing whether the procedure was succesful.
#' @export
#'
set_command_line_rp <- function(){

  copyright_text <- c("Note that the platform specific `rp` command line executables were created by Norbert Marwan and obtained under a Creative Commons License from the website of the Potsdam Institute for Climate Impact Research at: http://tocsy.pik-potsdam.de/ \n\n The full copyright statement on the website is as follows: \n\n  > \u00A9 2004-2017 SOME RIGHTS RESERVED  \n  > University of Potsdam, Interdisciplinary Center for Dynamics of Complex Systems, Germany  \n  > Potsdam Institute for Climate Impact Research, Transdisciplinary Concepts and Methods, Germany  \n  > This work is licensed under a [Creative Commons Attribution-NonCommercial-NoDerivs 2.0 Germany License](https://creativecommons.org/licenses/by-nc-nd/2.0/de/).  \n\n  More information about recurrence quantification analysis can be found on the [Recurrence Plot website](http://www.recurrence-plot.tk).")

  dl_instruction <- c("Download failed!\nCopying failed!\n\nTo install do the following:\n1. Either go to https://github.com/FredHasselman/casnet/tree/master/inst and download and unzip 'commandline_rp.zip', or go to http://tocsy.pik-potsdam.de/commandline-rp.php \n2. Find the executable for your OS\n3. Copy it to the '/exec' directory under 'casnet' \n4. Rename to 'rp' or 'rp.exe' on Windows\n5. Put an empty text file in '/exec' with the following name 'rp_instal_log.txt' \n 6. Run this code to test if everything is ok: crqa_cl(rnorm(100))")


  os                <- get_os()
  casnet_OS_options <- set_os_options()

  if(getOption("casnet.isRU")){
    URL <- "https://darwin.pwo.ru.nl/skunkworks/courseware/1718_DCS/crp_cl/windows_x86/rp_86.exe"
  } else {
    URL <- casnet_OS_options$URL
    }
  sys <- casnet_OS_options$sys
  exe <- casnet_OS_options$exe
  rp_prefix  <- casnet_OS_options$rep_prefix
  rp_command <- casnet_OS_options$rp_command
  sysdel     <- casnet_OS_options$sysdel
  syscopy    <- casnet_OS_options$syscopy

  execPath   <- getOption("casnet.path_to_rp")
  sourcePath <- getOption("casnet.path")

  TRYSYS <- FALSE

  # Check if the zipfile can be unpacked
  ZIP = FALSE
  if(file.exists(system.file("commandline_crp","cl_crp.zip", package="casnet"))){
    utils::unzip(system.file("commandline_crp","cl_crp.zip", package="casnet"),exdir=system.file("commandline_crp",package="casnet"))
    ZIP = TRUE
  }

  # Get the file from internet
  LOG <- try(utils::download.file(url = URL, mode = "wb", cacheOK = FALSE, destfile = normalizePath(paste0(execPath,"/",rp_command), mustWork = FALSE)))

  if(LOG==0){
    if(!os%in%"windows"){
      devtools::RCMD(cmd="chmod",options=paste0("a+x ",rp_command), path=normalizePath(execPath, mustWork = FALSE))
    }
    message(paste0("Detected: ",sys,"\n  Copied: ",URL," to ",rp_command," in ",execPath," as the commandline CRP executable"))
    rio::export(data.frame(url=c(URL,copyright_text)),normalizePath(paste0(execPath,"/rp_install_log.txt"), mustWork = FALSE))
  } else {
    if(ZIP){
      message(paste0("Detected: ",sys, "\n  FAILED to Copy: ",URL," to ",rp_command," in ",execPath," as the commandline CRP executable \n Trying .zip file..."))
      sysCommand <- c(syscopy,paste(normalizePath(paste0(sourcePath,"/commandline_rp/",sys,"/",exe)), normalizePath(paste0(execPath,"/",rp_command), mustWork = FALSE)),"chmod",paste("a+x ",normalizePath(paste0(execPath,"/",rp_command), mustWork = FALSE)))

      if(all(nchar(sysCommand)>0)){
        devtools::RCMD(sysCommand[[1]], options = sysCommand[[2]], path = sourcePath, quiet = TRUE)
        devtools::RCMD(sysCommand[[3]], options = sysCommand[[4]], path = sourcePath, quiet = TRUE)
        utils::write.table(data.frame(sysCommand=c(paste(sysCommand[[1]], sysCommand[[2]]), paste(sysCommand[[3]], sysCommand[[4]]),copyright_text)), normalizePath(paste0(execPath,"/rp_install_log.txt"), mustWork = FALSE))
      }
    } else {
      message(dl_instruction)
    }
  }

  message(paste0("\n==================COPYRIGHT=NOTICE==================\n\n",copyright_text,"\n\n==================COPYRIGHT=NOTICE==================\n"))
}

