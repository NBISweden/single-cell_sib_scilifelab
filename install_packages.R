##------------------------------------------------------------
## Installation script for Joint SIB / SciLifeLab
## Autumn School Single Cell Analysis 2019
##------------------------------------------------------------
##
## based on a script from Mike L. Smith created for CSAMA 2019
## https://www.huber.embl.de/users/msmith/csama2019/install_packages.R
##
## invoke by:
## source("https://raw.githubusercontent.com/NBISweden/single-cell_sib_scilifelab/master/install_packages.R")


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
.course_repos             <- "https://sib-course.s3.eu-west-2.amazonaws.com/repo/"
.pkg_file                 <- "https://raw.githubusercontent.com/NBISweden/single-cell_sib_scilifelab/master/required_packages.txt"

mem_pattern = "[0-9]+"
min_mem = 4 # gigabytes


##--------------------------------------
## Obtain required packages/code to run
##--------------------------------------
lns <- readLines(.pkg_file)
lns.pos <- match(c("## START-PACKAGES:", "## END-PACKAGES",
                   "## START-EXPRESSIONS:", "## END-EXPRESSIONS"),
                 lns)
if (any(is.na(lns.pos)) || any(diff(lns.pos) < 0)) {
    stop(.pkg_file, " does not contain expected anchors; please contact the course organizers.")
}
lns.pkgs <- unique(lns[seq(lns.pos[1] + 1L, lns.pos[2] - 1L)])
deps <- data.frame(name = sub(x = lns.pkgs, "^([[:alnum:]]+/)?([a-zA-Z0-9.]+)(-([0-9.]+))?$", "\\2"),
                   source = sub("-[0-9.]+$", "", lns.pkgs),
                   min.version = sub(x = lns.pkgs, "^([[:alnum:]]+/)?([a-zA-Z0-9.]+)(-([0-9.]+))?$", "\\4"),
                   stringsAsFactors = FALSE)
cmds <- unique(lns[seq(lns.pos[3] + 1L, lns.pos[4] - 1L)])


##---------------------------
## Install BiocManager
##---------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE))
    BiocManager::install("remotes", quiet = TRUE, update = FALSE, ask = FALSE, site_repository = .course_repos)
if (!requireNamespace("Biobase", quietly = TRUE))
    BiocManager::install("Biobase", quiet = TRUE, update = FALSE, ask = FALSE, site_repository = .course_repos)
options(warn = 1)


##-------------------------------------------
## installation function
##-------------------------------------------
installer_with_progress <- function(pkgs) {

    if (length(pkgs) == 0) { invisible(return(NULL)) }

    if (!requireNamespace("progress", quietly = TRUE)) {
        suppressMessages(
            BiocManager::install('progress', quiet = TRUE, update = FALSE, ask = FALSE, site_repository = .course_repos)
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
                                                       site_repository = .course_repos,
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
                   res <- system('grep "^MemTotal" /proc/meminfo', intern = TRUE)
                   as.numeric(regmatches(res, regexpr(mem_pattern, res))) / 10^6
               } else {
                   if (file.exists("/usr/sbin/system_profiler")) {
                       ## try MAC os
                       res <- system('/usr/sbin/system_profiler SPHardwareDataType | grep "Memory"', intern = TRUE)
                       as.numeric(regmatches(res, regexpr(mem_pattern, res)))
                   } else NULL
               },
           windows =
               tryCatch({
                   res = system("wmic ComputerSystem get TotalPhysicalMemory", ignore.stderr = TRUE, intern = TRUE)[2L]
                   as.numeric(regmatches(res, regexpr(mem_pattern, res))) / 10^9
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
R_version = paste(R.version$major, R.version$minor, sep = ".")
if ( !(R_version %in% .required_R_version) )
    stop(sprintf("You are using R-%s, which is not the one required for the course.\nPlease install R-%s\nR can be downloaded from here: %s",
                 R_version, .required_R_version[length(.required_R_version)], .r_url))


## Check Rstudio version
hasApistudio = suppressWarnings(require("rstudioapi", quietly = TRUE))
if ( !hasApistudio ) {
    BiocManager::install("rstudioapi", update = FALSE, site_repository = .course_repos, quiet = TRUE, ask = FALSE)
    suppressWarnings(require("rstudioapi", quietly = TRUE))
}

.rstudioVersion = try( rstudioapi::versionInfo()$version, silent = TRUE )
if ( inherits( .rstudioVersion, "try-error" ) ) {
    .rstudioVersion = gsub("\n|Error : ", "", .rstudioVersion)
    rstudioError = sprintf("The following error was produced while checking your Rstudio version: \"%s\"\nPlease make sure that you are running this script from an Rstudio session. If you are doing so and the error persists, please contact the course organisers.\n", .rstudioVersion)
    stop( rstudioError )
}

if ( !(.rstudioVersion >= .required_rstudio_version ) ) {
    rstudioVersionError = sprintf("You are using a Rstudio v%s, which is not the one required for the course.\nPlease install Rstudio v%s or higher.\nThe latest version of Rstudio can be found here: %s",
                                  .rstudioVersion, .required_rstudio_version, .rstudio_url)
    stop( rstudioVersionError )
}


## Check BioC version
if ( BiocManager::version() != .required_Bioc_version )
    stop(sprintf("You are using Bioconductor %s, which is not the one required for the course.\nPlease install Bioconductor %s\nInstallation instructions for Bioconductor are available here: %s\n",
                 BiocManager::version(), .required_Bioc_version, .bioc_url))


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
if (.Platform$OS.type == "windows" || Sys.info()["sysname"] == "Darwin") {
    BiocManager::install(toInstall, ask = FALSE, quiet = TRUE, update = FALSE, site_repository = .course_repos)
} else {
    fail <- installer_with_progress(toInstall)
}


##---------------------------------------
## Checking minimal versions of packages
##---------------------------------------
if (any(needVersionCheck <- deps$min.version != "")) {
    for (i in which(needVersionCheck)) {
        if (packageVersion(deps[i, "name"]) < deps[i, "min.version"]) {
            message("installing newer version of ", deps[i, "name"])
            BiocManager::install(deps[i, "source"], ask = FALSE, quiet = TRUE, update = FALSE, site_repository = .course_repos)
        } else {
            message(deps[i, "name"], "-", packageVersion(deps[i, "name"]),
                    " already fulfills minimal version requirement (", deps[i, "min.version"], ")")
        }
    }
}

##----------------------------------
## Additional installation commands
##----------------------------------
fail.cmds <- NULL
if (length(cmds) > 0) {
    message("Preparing additional requirements for course")
    bp <- progress::progress_bar$new(total = length(cmds),
                                     format = "Processing :current of :total expressions (:percent ) - current: :cmd",
                                     show_after = 0,
                                     clear = FALSE)
    for (i in seq_along(cmds)) {
        bp$tick(tokens = list(cmd = cmds[i]))
        tryCatch(
            suppressMessages( res <- eval(parse(text = cmds[i])) ),
            error = function(e) { fail.cmds <- c(fail.cmds, cmds[i]) },
            warning = function(w) { fail.cmds <- c(fail.cmds, cmds[i]) }
        )
    }
    message("DONE!")
}


##-------------------------
## Feedback on installation
##---------------------------
if(all( deps$name %in% rownames(installed.packages()) )) {
    cat(sprintf("\nCongratulations! All packages were installed successfully :)\nWe are looking forward to seeing you at the course!\n\n"))
} else {
    notinstalled <- deps[which( !deps$name %in% rownames(installed.packages()) ), "name"]


    if ( .Platform$pkgType == "win.binary" & 'Rsubread' %in% notinstalled ) {
        cat("The windows binaries for the package 'Rsubread' are not available. However, this package is not absolutely necessary for the exercises. If this is the only package
    that was not installed, there is no reason to worry. \n")
    }

    cat(sprintf("\nThe following package%s not installed:\n\n%s\n\n",
                if (length(notinstalled) <= 1) " was" else "s were",
                paste( notinstalled, collapse = "\n" )))

    if ( .Platform$pkgType != "source" ) {
        message("Please try re-running the script to see whether the problem persists.")
    } else {
        install_command <- paste0("BiocManager::install(c('", paste(notinstalled, collapse = "', '"), "'), site_repository = '", .course_repos, "')")
        message("Please try running the following command to attempt installation again:\n\n",
                install_command, "\n\n")
    }

    #if ( .Platform$pkgType == "source" ) {
    #    message("Some of the packages (e.g. 'Cairo', 'mzR', rgl', 'RCurl', 'tiff', 'XML') that failed to install may require additional system libraries.*  Please check the documentation of these packages for unsatisfied dependencies.\n A list of required libraries for Ubuntu can be found at http://www.huber.embl.de/users/msmith/csama2019/linux_libraries.html \n\n")
    #}

    if (length(fail.cmds) > 0) {
        message("Some of the commands to download datasets failed:\n   ", paste(fail.cmds, collapse = "\n   "), "\n")
    }

    message("If you need help with troubleshooting, please contact the course organisers.")
}
