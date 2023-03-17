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

|	subcommands	|	description																						|
|	-----------	|	-----------																						|
|	`-i/--in`			|	*Required* Input filename. Format should be as a tab delimited table					|
|	`-o/--out`			|	Output filename. 																		|
|	`-s/--start`		|	Start date of analysis.	Truncates prior dated data.										|
|	`-e/--end`			|	End Date of analysis. Truncates proceeding data.										|
|	`-f/--future`		|	Future. Number of days to predict the trajectory into the future. 						|
|	`-b/--baseline`		|	*Required* Baseline clade as string. 													|
|	`-j/--subjects`		|	The subjects as a string to calculate the relative reproduction numbers from. 			|
|	`-w/--breaks`		|	A list of window dates to divide the analysis by. 										|
|	`-p/--precision`	|	Precision. Default 1e-4.  																|
|	`-l/--len`			|	Maximum length of time distribution for an infection. 									|
|	`-a/--alpha`		|	Alpha parameter for generation function. Default = 16, use 7 for H1-pdm influenza. 		|
|	`-t/--theta`		|	Theta parameter for generation function.  Default = 2.03, use 4.5 for H1-pdm influenza.	|
|	`-d/--delta`		|	Delta parameter. Default = 2.32, use 0.6 for H1-pdm influenza.							|
|	`-u/--unit`			|	*Required* unit of analysis as day, week, or month [D,W,M]. Default is :D/day.			|
|	`-q/--frequency`	|	True/False Predict a trajectory.														|
|	`-g/--estimate_GT`	|	True/False Estimate generation time. 													|
|	`-c/--estimate_CI`	|	True/False Estimate confidence intervals. 												|
