##------------------------------------------------------------
## Installation script for Joint SIB / SciLifeLab
## Autumn School Single Cell Analysis 2019
##------------------------------------------------------------
##
## based on a script from Mike L. Smith written for CSAMA 2019
## https://www.huber.embl.de/users/msmith/csama2019/install_packages.R


##-------------------------------------------
## System requirements
##-------------------------------------------
.required_R_version       <- c( "3.6.1" )
.required_Bioc_version    <- "3.9"
.Bioc_devel_version       <- "3.10"
.required_rstudio_version <- "1.2"
.r_url                    <- "https://stat.ethz.ch/CRAN/"
.rstudio_url              <- "https://www.rstudio.com/products/rstudio/download/"
.bioc_url                 <- "http://www.bioconductor.org/install/"
.course_repos             <- character(0)
.pkg_file                 <- "required_packages.md"

mem_pattern = "[0-9]+"
min_mem = 4 # gigabytes


##---------------------------
## Install BiocManager
##---------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE))
    BiocManager::install("remotes", quiet = TRUE, update = FALSE, ask = FALSE)
if (!requireNamespace("Biobase", quietly = TRUE))
    BiocManager::install("Biobase", quiet = TRUE, update = FALSE, ask = FALSE)
options(warn = 1)


##-------------------------------------------
## installation function
##-------------------------------------------
installer_with_progress <- function(pkgs) {
    
    if (length(pkgs) == 0) { invisible(return(NULL)) }
    
    if (!requireNamespace("progress", quietly = TRUE)) {
        suppressMessages(
            BiocManager::install('progress', quiet = TRUE, update = FALSE, ask = FALSE)
        )
    }
    
    toInstall <- pkgs
    bp <- progress::progress_bar$new(total = length(toInstall),
                                     format = "Installed :current of :total (:percent ) - current package: :package",
                                     show_after = 0,
                                     clear = FALSE)
    
    length_prev <- length(toInstall)
    fail <- NULL
    while (length(toInstall)) {
        pkg <- toInstall[1]
        bp$tick(length_prev - length(toInstall),  tokens = list(package = pkg))
        length_prev <- length(toInstall)
        if (pkg %in% rownames(installed.packages())) {
            toInstall <- toInstall[-1]
        } else {
            tryCatch(
                suppressMessages( BiocManager::install(pkg, quiet = TRUE, update = FALSE, ask = FALSE,
                                                       Ncpus = parallel::detectCores() ) ),
                error = function(e) { fail <<- c(fail, pkg) },
                warning = function(w) { fail <<- c(fail, pkg) },
                finally = toInstall <- toInstall[-1]
            )
        }
    }
    bp$tick(length_prev - length(toInstall),  tokens = list(package = "DONE!"))
    
    return(fail)
}


## Check memory size
mem <-
    switch(.Platform$OS.type,         
           unix = 
               if (file.exists("/proc/meminfo")) {
                   ## regular linux
                   res <- system('grep "^MemTotal" /proc/meminfo', intern=TRUE)
                   as.numeric(regmatches(res, regexpr(mem_pattern, res)))/10^6
               } else {
                   if (file.exists("/usr/sbin/system_profiler")) {
                       ## try MAC os
                       res <- system('/usr/sbin/system_profiler SPHardwareDataType | grep "Memory"', intern=TRUE)
                       as.numeric(regmatches(res, regexpr(mem_pattern, res)))
                   } else NULL  
               },
           windows = 
               tryCatch({
                   res = system("wmic ComputerSystem get TotalPhysicalMemory", ignore.stderr=TRUE, intern=TRUE)[2L]
                   as.numeric(regmatches(res, regexpr(mem_pattern, res)))/10^9
               }, error = function(e) NULL),
           NULL)

if (is.null(mem)) {
    warning(sprintf("Could not determine the size of your system memory. Please make sure that your machine has at least %dGB of RAM!", min_mem))
} else {
    mem = round(mem)
    if ( mem < min_mem ) stop(sprintf("Found %dGB of RAM. You need a machine with at least %dGB of RAM for the course!", mem, min_mem))
    else message(sprintf("Found %dGB of RAM", mem))
}


## Check the R version
R_version = paste(R.version$major, R.version$minor, sep=".")
if( !(R_version %in% .required_R_version) )
    stop(sprintf("You are using R-%s, which is not the one required for the course.\nPlease install R-%s\nR can be downloaded from here: %s",
                 R_version, .required_R_version[length(.required_R_version)], .r_url))


## Check Rstudio version
.baseurl <- unique(c(BiocManager::repositories(), .course_repos))
hasApistudio = suppressWarnings(require("rstudioapi", quietly=TRUE))
if( !hasApistudio ){
    BiocManager::install("rstudioapi", update=FALSE, siteRepos = .baseurl, quiet = TRUE, ask = FALSE)
    suppressWarnings(require("rstudioapi", quietly=TRUE))
}

.rstudioVersion = try( rstudioapi::versionInfo()$version, silent=TRUE )
if( inherits( .rstudioVersion, "try-error" ) ){
    .rstudioVersion = gsub("\n|Error : ", "", .rstudioVersion)
    rstudioError = sprintf("The following error was produced while checking your Rstudio version: \"%s\"\nPlease make sure that you are running this script from an Rstudio session. If you are doing so and the error persists, please contact the course organisers.\n", .rstudioVersion)
    stop( rstudioError )
}

if( !( .rstudioVersion >= .required_rstudio_version ) ){
    rstudioVersionError = sprintf("You are using a Rstudio v%s, which is not the one required for the course.\nPlease install Rstudio v%s or higher.\nThe latest version of Rstudio can be found here: %s",
                                  .rstudioVersion, .required_rstudio_version, .rstudio_url)
    stop( rstudioVersionError )
}


## Check BioC version
if( BiocManager::version() != .required_Bioc_version )
    stop(sprintf("You are using Bioconductor %s, which is not the one required for the course.\nPlease install Bioconductor %s\nInstallation instructions for Bioconductor are available here: %s\n",
                 BiocManager::version(), .required_Bioc_version, .bioc_url))


## Get list of packages to install
lns <- readLines(.pkg_file)
lns <- sub("^ *- *", "", lns)
deps <- lns[grep("^[a-zA-Z0-9.]+$", lns)]
deps <- data.frame(name = gsub(x = deps, "^[[:alnum:]]+/", ""),
                   source = deps, 
                   stringsAsFactors = FALSE)


## omit packages not supported on WIN and MAC
#type = getOption("pkgType")
#if ( type == "win.binary" || type == "mac.binary" ) {
#  deps = setdiff(deps, c('gmapR'))
#}
toInstall = deps[which( !deps$name %in% rownames(installed.packages())), "source"]


## set up directory where downloaded packages are stored
destdir = NULL


## do not compile from sources
options(install.packages.compile.from.source = "never")
if(.Platform$OS.type == "windows" || Sys.info()["sysname"] == "Darwin") {
    BiocManager::install(toInstall, ask = FALSE, quiet = TRUE, update = FALSE)
} else {
    fail <- installer_with_progress(toInstall)
}


##---------------------------
## Additional installation commands
##---------------------------
# if(Biobase::package.version("msdata") < "0.24.1") {
#     BiocManager::install(msdata, ask = FALSE, quiet = TRUE, update = FALSE)
# }
# 
# TENxPBMCData::TENxPBMCData(dataset = "pbmc3k")
# TENxPBMCData::TENxPBMCData(dataset = "pbmc4k")
# 
# library("ensembldb")
# ah <- AnnotationHub::AnnotationHub()
# ah_91 <- AnnotationHub::query(ah, "EnsDb.Hsapiens.v91")
# edb <- ah[[names(ah_91)]]


##-------------------------
## Feedback on installation
##---------------------------
if(all( deps$name %in% rownames(installed.packages()) )) {
    cat(sprintf("\nCongratulations! All packages were installed successfully :)\nWe are looking forward to seeing you at the course!\n\n"))
} else {
    notinstalled <- deps[which( !deps$name %in% rownames(installed.packages()) ), "name"]
    
    
    if( .Platform$pkgType == "win.binary" & 'Rsubread' %in% notinstalled ){
        cat("The windows binaries for the package 'Rsubread' are not available. However, this package is not absolutely necessary for the exercises. If this is the only package
    that was not installed, there is no reason to worry. \n")
    }
    
    cat(sprintf("\nThe following package%s not installed:\n\n%s\n\n", if (length(notinstalled)<=1) " was" else "s were", paste( notinstalled, collapse="\n" )))
    
    if( .Platform$pkgType != "source" ){
        message("Please try re-running the script to see whether the problem persists.")
    } else {
        install_command <- paste0("BiocManager::install(c('", paste(notinstalled, collapse = "', '"), "'))")
        message("Please try running the following command to attempt installation again:\n\n",
                install_command, "\n\n")
    }
    
    message("If you need help with troubleshooting, please contact the course organisers.")
    
    #if( .Platform$pkgType == "source" ){
    #    message("Some of the packages (e.g. 'Cairo', 'mzR', rgl', 'RCurl', 'tiff', 'XML') that failed to install may require additional system libraries.*  Please check the documentation of these packages for unsatisfied dependencies.\n A list of required libraries for Ubuntu can be found at http://www.huber.embl.de/users/msmith/csama2019/linux_libraries.html \n\n")
    #}
}
