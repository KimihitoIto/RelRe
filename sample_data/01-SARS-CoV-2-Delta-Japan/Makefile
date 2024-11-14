SHAPE=3.42
SCALE=1.36
# Nishiura, H.; Linton, N.M.; Akhmetzhanov, A.R. Serial interval of novel coronavirus (COVID-19) infections. Int. J. Infect. Dis. 2020, 93, 284–286.

all:RelRe.jl estimates.csv frequencies.pdf average_k.pdf

RelRe.jl:
	ln -s ../../$@ .

estimates.csv:Japan_Linenage_Counts.csv
	julia --threads 10 RelRe.jl -b other -a $(SHAPE) -t $(SCALE) -i Japan_Linenage_Counts.csv -c -q -f 90

frequencies.csv:estimates.csv

loglikelihood.csv:estimates.csv

frequencies.pdf:frequencies.csv
	R -f plot_frequency.R

average_k.pdf:frequencies.csv
	R -f plot_average_k.R

clean:
	rm -f estimates.csv
	rm -f frequencies.csv
	rm -f loglikelihood.csv
	rm -f Dirichlet.csv
	rm -f frequencies.pdf
	rm -f average_k.pdf
