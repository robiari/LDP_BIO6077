# this script initializes the project using renv() I suggest running it one line
# at a time and thinking about the need to "restore". For example, if all of the
# scripts run without error and give the same results then you can probably just
# work with what you have. If there are errors you can "restore" and see if that
# gets ride of the errors.

# first let's check to see what, if anything, is out of sync
renv::status()

# we can either opt to work with what we have OR restore our
# environment to the same state that the code was developed under.

# If we want to work in the same state as development run the following to make
# sure we have the same packages as renv.lock file
renv::restore()

# Note that this is far from perfect!!! Some packages need to be compile and
# that does not always work. You can try excluding the package from the restore
# but... that does not always work because of dependencies between packages. Two
# common problem packages are MASS ad Matrix. You can try this to see if it
# works.
renv::restore(exclude = c("MASS", "Matrix"))

# if that does not work you can try this:
#require(devtools)
#install_version("Matrix", version = "1.5-1", repos = "http://cran.us.r-project.org")

# or go to the cran archive, get the url for the version you want and cross your fingers.
# install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.5-1.tar.gz", repos=NULL, type="source")

# If the above "things" don't work to resolve your package installation issues
# you should try updating to the latest version of R and installing the gfortran
# compiler (Mac) and / or RTools which is available from here: 
# Mac - https://mac.r-project.org/tools/
# Windows - https://cran.r-project.org/bin/windows/Rtools/

# to update renv.lock:
renv::snapshot()
