ALPHA=3.42
THETA=1.36
# Nishiura, H.; Linton, N.M.; Akhmetzhanov, A.R. Serial interval of novel coronavirus (COVID-19) infections. Int. J. Infect. Dis. 2020, 93, 284–286.

all:RelRe.jl estimates.csv frequencies.pdf average_k.pdf

RelRe.jl:
	ln -s ../../$@ .

estimates.csv:England_Alpha_Delta.csv
	julia --threads 10 RelRe.jl -b Alpha -a $(ALPHA) -t $(THETA) -i England_Alpha_Delta.csv -c -q -D

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
