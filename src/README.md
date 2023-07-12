# RelRe
## H1N1 variation file


## Usage example

```sh
julia --threads 10 RelRe.jl -c data.csv -u M -b 6B.1A.5a1 -a 4.5 -t 0.6 -l 7 -s 2019-10-01 -e 2022-03-01 -w 2020-04-01 2020-10-01 2021-04-01 2021-10-01 -p 1e-4
```

> **Note**
> This version of the code is optimized for Sars-Cov-2. Parameters should be adjusted as needed.

## Params

Currently the following parameters are recognized

|	parameter	| variable |description																									|																					
|	-----------	|----------|-----------																									|																					
| `-a`, `--alpha` 		| ALPHA     	| shape parameter of gamma distribution for generation time (type: Float64, default: 2.03)		|
| `-b`, `--baseline` 	| BASELINE		| variant used as the baseline of relative                                                      |
| `-c`, `--estimate_CI` |   			| estimate 95% confidence intervals                                                             |
| `-d`, `--delta` 		| DELTA    		| unit time of calculation (in days) (type:Float64, default: 0.5)                               |
| `-D`, `--Dirichlet`   |  		 		| use Dirichlet multinomial as the observation model                                            |
| `-e`, `--end` 		| END         	| end date of the analysis (default: "")                                                        |
| `-f`, `--future` 		| FUTURE   		| duration in days for predicting variant frequencies (type: Int64, default: 0)                 |
| `-g`, `--estimate_GT` |   			| estimate relative generation times of variants                                                |
| `-h`, `--help`        |   			| show this help message and exit                                                               |
| `-i`, `--in` 			| IN 			| input file containing temporal count data of variants                                         |
| `-j`, `--subjects` 	| [SUBJECTS...]	| list of variants to calculate relative reproduction numbers (type: Symbol)                    |
| `-l`, `--len` 		| LEN         	| trancation point of gamma distribution for generation time (type: Int64, default: 16)         |
| `-n`, `--undetected`  |  				| assume all variants exist undetected from the start date                                      |
| `-o`, `--out` 		| OUT         	| prefix of output files (default: "")                                                          |
| `-p`, `--precision` 	| PRECISION		| stopping criterion used as ftol_abs in NLopt (type: Float64, default: 0.0001)                 |
| `-q`, `--frequency`   |   			| calculate the time course of variant frequencies                                              |
| `-s`, `--start` 		| START     	| start date of the analysis (default: "")                                                      |
| `-t`, `--theta`		| THETA     	| scale parameter of gamma distribution for generation time (type: Float64, default: 2.32)      |
| `-u`, `--unit` 		| UNIT       	| unit time of observations: D (Daily),W (Weekly),or M (Monthly) (type: Symbol, default: :D)    |
