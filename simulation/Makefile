ALPHA=4.109
THETA=0.73

all:RelRe.jl counts.csv pop_freq.pdf	

RelRe.jl:
	ln -s ../$@ .

counts.csv:parameters.csv
	julia RelRe_simul.jl -b v001 -a $(ALPHA) -t $(THETA) -i $< -d 1 -e 2024-06-30 

estimates.csv:count.csv
	julia --threads RelRe.jl -b v001 -a $(ALPHA) -t $(THETA) -i $< -d 1 -c 	

pop_freq.pdf:pop_freq.csv
	R -f plot_pop_freq.R

clean:
	rm -f *counts.csv
	rm -f *pdf 
	rm -f *pop_freq.csv
