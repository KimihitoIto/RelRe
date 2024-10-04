ALPHA=4.109
THETA=0.73
THREADS=16
N=3
n=100
all:RelRe.jl counts.csv pop_freq.pdf estimates.csv parameters.csv	
RelRe.jl:
	ln -s ../$@ .

parameters.csv:
	julia params_gen.jl $(N)

counts.csv:parameters.csv
	julia RelRe_simul.jl -a $(ALPHA) -t $(THETA) -i $< -d 1 -e 2024-06-30 -n $(n) 

estimates.csv:counts.csv
	julia --threads=$(THREADS) RelRe.jl -b v000 -a $(ALPHA) -t $(THETA) -i $< -d 1 -u 
#-c 	

pop_freq.csv:counts.csv

pop_freq.pdf:pop_freq.csv
	R -f plot_pop_freq.R

clean:
	rm -f *counts.csv
	rm -f *pdf 
	rm -f *pop_freq.csv
	rm -f *loglikelihood.csv
	rm -f *estimates.csv
	rm -f parameters.csv
