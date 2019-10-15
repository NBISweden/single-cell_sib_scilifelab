# <img border="0" src="https://www.svgrepo.com/show/19652/maths-class-materials-cross-of-a-pencil-and-a-ruler.svg" width="40" height="40"> Pre-course Materials

***

<br/>


### <img border="0" src="https://www.svgrepo.com/show/4795/installation-symbol.svg" width="40" height="40"> Package Installation

***

In preparation for the practical sessions, you should install required software on the computer that you will bring to the autumn school.

For this, please follow the steps below:
1. you need write permission to your R library folder (typically `.libPaths()[1]`) and enough free space (>2GB)
2. you need to be connected to the internet 
3. start R-Studio on your computer
4. run the following command on the R console:

```
source("https://raw.githubusercontent.com/NBISweden/single-cell_sib_scilifelab/master/install_packages.R")
```

The script will check your environment, install what is missing, and report any problems.
During the installation, you may be prompted to update existing packages (recommended).

Don’t worry if not everything is installed successfully - please contact the course organisers for help, either by email or on Sunday before dinner.

We are looking forward to meet you all in Leysin!

<br/>

***

Here are a couple of common issues and their solutions:  

- package fails to compile/install on a mac with an error of:  
  `clang: error: unsupported option '-fopenmp'`. This error can be fixed by installing
  an up-to-date version of compilers:  
  * - go to https://stat.ethz.ch/CRAN/ -> download R for mac, follow the link to ‘tools’ (in the ‘Important’ paragraph), and download and install clang-7.0.0 and gfortran-6.1 (for R 3.6). 

  * modify your ~/.R/Makevars (or create it if it does not exist yet) to include the following: 
```
FLIBS=""  
F77="/usr/local/gfortran/bin/gfortran"
FC="/usr/local/gfortran/bin/gfortran"

CC=/usr/local/clang7/bin/clang
CXX=/usr/local/clang7/bin/clang++
CXX11=/usr/local/clang7/bin/clang++
CXX14=/usr/local/clang7/bin/clang++
CXX17=/usr/local/clang7/bin/clang++
CXX1X=/usr/local/clang7/bin/clang++
LDFLAGS=-L/usr/local/clang7/lib
```

   * then try the installation again.

- the package `HDF5Array` fails to compile/install with a `file not found` error (files from `Rhdf5lib`). This is caused by recent changes in `Rhdf5lib` (1.6.2). `HDF5Array` (1.12.3) has been fixed, but is not yet available from bioconductor.org, but it can be downloaded/installed directly from github:   
```
git clone --branch RELEASE_3_9 https://git.bioconductor.org/packages/HDF5Array
R INSTALL HDF5Array
```

## #[Back to main](README.md)
