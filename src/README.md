# RelRe
## H1N1 variation file


## TODO and Issues

* Currently all strains are being used as subject. Need to allow this to be subset, as the confidence interval calculation can be lengthy.
* Code needs to be merged to minimize cross file dependencies.
* Other segments of the code can be combined under a single file.
* Correct variable scoping.
 
## Usage example

```sh
julia --threads 5 advantage-GT-Rt-mtly.jl -m yearseason_freq.csv -b 6B.1 -s 2020-01-01 -e 2021-01-01
```

> **Note**
> This version of the code is optimized for human H1N1 sequences. Alpha and theta parameters need to be fitted for other viruses. 

## Params

Currently the following parameters are recognized

|	subcommands	|	description								|
|	-----------	|	-----------								|
|	`-m`		|	Matrix file as a path. *Required*		|
|	`-b`		|	Baseline clade as string. *Required*	|
|	`-p`		|	Precision. Default 1e-4. 				|
|	`-s`		|	Start date of analysis. *Required*		|
|	`-e`		|	End Date of analysis. *Required*		|
