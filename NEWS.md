# `aftgee` 1.2.0
	* Faster ns estimator for `aftsrr()`
	* Faster GEE in `aftgee()` when `margin` is all 1
	* Faster resampling with GEE in `aftgee()` when `margin` is all 1
	* Export `gee()` when `family = gaussian()`
  * Fix the print in ini.conv
  * Add S3 method confint()
# `aftgee` 1.1.6
  * Fixed examples that depends on `kidney` from survival package, which no longer exist
  * Change 'http' to 'https'
# `aftgee` 1.1.5
  * Add QIC
  * Fixed bugs from `non-conformable arguments`  
# `aftgee` 1.1.4
  * Fixed `PW` and `GP` weights
  * Added `gp.pwr` for `rankWeights = GP` so that `PW` and `logrank` can be called with `rankWeights = GP` with `gp.pwr = 1` and `gp.pwr = 0`, respectively.
# `aftgee` 1.1.3
  * rebuild with oxygen2
  * hosted on GitHub
# `aftgee` 1.1.2
  * Fixed memory leak messages
  * Use S4 functions to set and dispatch methods in aftsrr
# `aftgee` 1.1.1
  * Version updated to published version with changes made based on comments from the JSS paper