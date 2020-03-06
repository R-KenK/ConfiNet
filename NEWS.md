# ConfiNet 0.7.1
* Added progress bar through pblapply (pbapply package) to Boot_scan()

# ConfiNet 0.7.0
* Implemented focal scan with definable probs by giving output option to do.scan()
* fixed typo in previous NEWS.md
* added Boot.focal_scan() function (quick implementation, require some proofing)
* added scaling option to both group and focal bootstrap functions
* merged Boot.group_scan() and Boot.focal_scan() into Boot_scan() with "method" option
* Wrapped previous Boot_scan into higher bootstrap function to avoid multiple makeCluster()/stopCluster() call (huge impact on perfs otherwise), renamed Boot_scan() into iterate_scans() for better clarity, and moved some functional tests into Boot_scans() rather than iterate_scans()
* added irrelevant NULLs assignement to dummy variables to remove NOTEs in R CMD check (till better fix emerge in the R community...)

# ConfiNet 0.6.1
* Bug fixes and renaming: scale.to.binary.prob -> Binary.prob, moved simu.scan into .WIP folder
* other minor syntax/bug fixes till build doesn't return any error.

# ConfiNet 0.6.0
* Added focal scan simulation function simu.focal (need cleaning)
* Pre-implementation of user-defined probabilities for simu.scan and simu.focal
* Added option to output either list of scans or weighted adjacency matrix in simu.scan
* implemented scale.to.binary.prob (will haev to rename probably) and do.scan

# ConfiNet 0.5.1
* Added group scan simulation function simu.scan

# ConfiNet 0.5.0
* Change package name KenNet -> ConfiNet + script and function names changes (KenNet -> ConfiNet and KenBoot -> ConfiBoot).

# ConfiNet 0.4.1
* Bug fixes related to Documents and Check

# ConfiNet 0.4.0
* Reading through Whitehead (2008) to implement the different association indices (WIP)...
* Implementing adjacency to association index for the case of focal sampling and egocentric networks

# ConfiNet 0.3.1
* Implemented edge -> adj, probably useful even if only for internal use...

# ConfiNet 0.3.0
* Progress toward bootstrap and network building implementation (still WIP).

# ConfiNet 0.2.0
* Preparations for RTWS algorithm: (rk.)closest.dates, extend.window and random.dates.

# ConfiNet 0.1.3
* Added Adj -> edge list with repetitions representing edge weights.

# ConfiNet 0.1.2
* Added Obs list -> list of subsets to prepare for automatic dynamic SNA.

# ConfiNet 0.1.1
* Formatting everything related to versioning


# ConfiNet 0.1.0
* Added a `NEWS.md` file to track changes to the package.
