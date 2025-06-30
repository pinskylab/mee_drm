## `drmr`: A Bayesian approach to Dynamic Range Models in R

This repository provides the R code and data to reproduce the results for the
paper, "`drmr`: A Bayesian approach to Dynamic Range Models in R."


### Installation

Before running any analysis, please follow these steps to set up the required
environment.

1.  **Install CmdStan**: This project depends on `CmdStan` (version \>=
    2.36). Please install it and the R package `cmdstanr` by following the
    [official
    instructions](https://mc-stan.org/cmdstanr/articles/cmdstanr.html).

2.  **Restore R Environment**: This repository uses `renv` to manage R package
    dependencies (including the version of `drmr` in the `pkg/` directory). To
    install all required packages, open this project in R and run:y
    ```r
    renv::restore()
    ```

3.  **System Dependencies (Optional)**: If you wish to regenerate the raw
    environmental data for the red-bellied woodpecker case study, you must also
    have [GDAL](https://gdal.org/en/stable/) installed on your system. This is
    not required if you use the provided, pre-processed data.


### Reproducing the Results

The analyses for the two case studies can be run using the scripts in the
`code/` directory.

#### Case Study 1: Summer Flounder

To reproduce the full analysis for the summer flounder, run the
`code/01-summer-flounder.R` script.y


#### Case Study 2: Red-bellied Woodpecker

For convenience, the final processed dataset for this study is already provided
at `data/birds/processed.parquet`.

  * **To run the analysis using the provided data**, execute the `code/04-rbw.R`
    script.
  
  * **(Optional) To regenerate the processed data from raw sources**, run the
    following scripts in order:
    1.  `code/02-download-bbs.R`: Downloads Breeding Bird Survey (BBS) data via
        the [`bbsBayes2` package](https://github.com/bbsBayes/bbsBayes2).
    2.  `chelsa-download.sh`: A `bash` script to download and crop
        [CHELSA](https://chelsa-climate.org/) air temperature data.
    3.  `code/03-chelsa-env4bbs.R`: Aggregates the temperature data and creates
        the final `processed.parquet` file.


### Notes on Reproducibility

Please be aware that results from `CmdStan` may vary slightly across different
operating systems and hardware platforms due to minor differences in random
number generator implementations. Full reproducibility is only guaranteed when
using an identical computational environment.

#### System information

The analyses were conducted on a macOS machine with the following specifications:

  * **Processor**: Apple M2 Pro
  * **Memory**: 16GB RAM

The R session information is provided below for full reproducibility.
```
> sessionInfo()
#> R version 4.4.1 (2024-06-14)
#> Platform: aarch64-apple-darwin20
#> Running under: macOS 15.5
#>
#> Matrix products: default
#> BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Vers#> ions/A/libBLAS.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
```
