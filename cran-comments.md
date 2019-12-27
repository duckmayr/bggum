## Resubmission

This is a resubmission. In this version I have:

* Added a citation to the paper developing the generalized graded unfolding
  model (GGUM -- the statistical model our package implements) to the
  DESCRIPTION file, as requested by Swetlana Herbrandt.
* Removed all instances of \dontrun{}. Swetlana Herbrandt suggested either
  replacing them with \donttest{} or unwrapping any examples that can be
  executed in less than 5 seconds per Rd-file. As originally constituted,
  not all examples could have done that, but I scaled them down a bit and
  they easily passed that mark by far.
* Unrelated to requests by the CRAN team, when editing the documentation as
  described above, I found various small tweaks to improve clarity. None of
  the R or C++ functions have been altered, just tweaks to the documentation.

## Test environments

* local Manjaro Linux 18.1.4 install, R 3.6.1, 3.6.2, and devel
* Ubuntu 16.04.6 (on Travis CI), R 3.5.3, 3.6.1, and devel
* macOS High Sierra 10.13.6 (on Travis CI), R 3.5.3 and 3.6.2
* Windows Server 2012 R2 (on AppVeyor) R 3.6.2 Patched (2019-12-26 r77621)
* win-builder (R-devel and 3.6.2)

## R CMD check results

0 errors | 0 warnings | 1 note

* The only note was:
  
checking CRAN incoming feasibility ... NOTE

Maintainer: 'JBrandon Duck-Mayr <j.duckmayr@gmail.com>'

New submission

Possibly mis-spelled words in DESCRIPTION:
  Donoghue (23:60)
  Laughlin (24:5)

So there was, of course, the new submission note.
The flagged words for potential misspelling are last names and are not
misspelled (they are last names involved in the citation in the DESCRIPTION
file mentioned above).

## Reverse dependencies

This is a new release, so there are no reverse dependencies.

