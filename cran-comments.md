## Test environment
* local Ubuntu 18.04.2 LTS
* devtools::check_rhub()
* devtools::check_win_devel()

## R CMD CHECK results
There were no ERRORS or NOTES

There was 1 WARNING:
  Compilation used the following non-portable flag(s):
    ‘-Wdate-time’ ‘-Werror=format-security’ ‘-Wformat’

I do not find the file I need to chenge to get rid of this error.

## Downstream dependencies

FastWilcoxTest only depends on cran software. Should I really test them, too?
I get the error message

revdepcheck::revdep_check()
── INSTALL ─────────────────────────────────────────────────────── 2 versions ──
Installing CRAN version of FastWilcoxTest
Error: (converted from warning) package ‘FastWilcoxTest’ is not available (for R version 3.5.3)

I assume I do not have downstream dependencies?

