all:RelRe.jl estimates.csv frequencies_bar.pdf average_k.pdf

RelRe.jl:
	ln -s ../../$@ .

estimates.csv:Denmark_all_lineages.csv
	julia --threads 16 RelRe.jl -b Delta -a 2.03 -t 2.32 -i Denmark_all_lineages.csv -c -q -f 30 -g -D

frequencies.csv:estimates.csv

loglikelihood.csv:estimates.csv

Dirichlet.csv:estimates.csv

frequencies_bar.pdf:frequencies.csv
	R -f plot_frequency_bar.R

average_k.pdf:frequencies.csv
	R -f plot_average_k.R

clean:
	rm -f estimates.csv
	rm -f frequencies.csv
	rm -f loglikelihood.csv
	rm -f Dirichlet.csv
	rm -f *.pdf
	rm -f *.svg
