
## Test environments
* local OS X install, R 3.6.0
* ubuntu 16.04 (on travis-ci), R 3.6.0
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

New maintainer:
  Robrecht Cannoodt <rcannood@gmail.com>
Old maintainer(s):
  ORPHANED

## Reverse dependencies

R CMD check was run on 10 Sep 2019 using the following command:

```
revdepcheck::revdep_check(timeout = as.difftime(60, units = "mins"), num_workers = 8)
```

Summary:
```
── CHECK ────────────────────────────-----─────────────── 2 packages ──
✔ dimRed 0.2.3                           ── E: 2     | W: 0     | N: 1                                                                             
I dyndimred 1.0.2                        ── E: 1     | W: 0     | N: 0                                                                             
OK: 2                                                                                                                                            
BROKEN: 0
Total time: 11 min
```
