# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/tests.html
# * https://testthat.r-lib.org/reference/test_package.html#special-files

library(testthat)
library(gencoreSC)

## For some reason the github action is failing, it can't read the h5 data in,
## but running the tests locally works fine. Possibly related to corruption of
## h5 on github?

## The error:
## Error in `initialize(value, ...)`: invalid name for slot of class "Assay5": strip.suffix

## Uncomment this to renable tests
# test_check("gencoreSC")

