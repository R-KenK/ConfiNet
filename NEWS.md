# ConfiNet (development version)

# ConfiNet 0.10.1
* Cleaned expected time modelling of standard and optimized method, use the glm models as internal data.
* Proper implementation of the decide_use.rare.opti() using those models.

# ConfiNet 0.10.0
* HUGE revamp from the base up to implement diverse optimization according to succesive benchamrk results. Notably relying solely on rbinom, implementing an alternative algorithm in the case of rare event (cf. Optimization for rare events.R), and homogenization of outputs according to the different method (avoiding additional list levels).
* Also now benefitting from the use of "..." to clean code, although R CMD check require nullifying a few optional parameters...
* few checks and need figuring out the expected time according to use.rare.opti, then should be pushable to master. 

# ConfiNet 0.9.1
* Revamp and huge code cleaning of simulation script. Now, considering only an unested list per variable, the parameter combinations are automatically and flexibly generated, and cleanly integrated in the different bootstraps.
* added Boot_get.param() function to retrieve the stored parameters as a data.frame of factors, to ease further data gathering
* added adjacency_cor() to calculate the correlation coefficient between "flattened" adjacency matrices
* added boot_progress.param() progress report for the iteration through the parameters.list
* prototype of data collection and graphical exploration in .WIP/Simulation_script.R

# ConfiNet 0.9.0
* Added sections to the matrix.tool.R script
* Added obs.prob_bias() as a wrapper to implement observation probability biases (either trait- or network-based), via the use of user-defined functions
* Considering how flexible the obs.prob_bias() function is, added relevant if statements in observable_edges() to test if input obs.prob matrice or vector contains negative values, and to scale its values to ]0,1[ (zeros and ones are avoided, because of their poor simulation relevance, according the the smallest non-null value of obs.prob)
* Fixed previous code error in sum_up.scans() (line 32 of b63fcc6ece0476362daa6e8ab9a8534357d68617) for specific parameter choice.

# ConfiNet 0.8.1
* Structuring simulation script and proof of concept for multiple parameters testing
* Added proper handling of adjacency mode in the form of adjacency_mode() and binary_adjacency_mode(), included in do.scan() and sum_up.scans()
* Added a wrapper for the frequently used do.call(rbind,lapply()) sequence

# ConfiNet 0.8.0
* Added observable_edge() to subset group scans according to observation probability. Supports probability matrix, vector, or single value.
* Added output="all" option to keep track of both list and adjacency matrix in Boot_scans()
* Updated consequently nested functions
* Added internal function Bootstrap_add.attributes() to store Boot_scans() parameters as attributes in the Bootstrap output
* Added Boot_get.list() function to cleanly retrieve some data of a Boot_scans() output. Relies on the Bootstrap attributes to avoid users inputting them

# ConfiNet 0.7.1
* Added progress bar through pblapply (pbapply package) to Boot_scan()

# ConfiNet 0.7.0
* Implemented focal scan with definable probs by giving output option to do.scan()
* fixed typo in previous NEWS.md
* Added Boot_scan into higher bootstrap function with method option (group or focal scans)

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
