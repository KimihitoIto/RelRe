all:RelRe.jl estimates.csv frequencies.pdf average_k.pdf

RelRe.jl:
	ln -s ../../$@ .

estimates.csv:Denmark_Delta_Omicron.csv
	julia --threads 10 RelRe.jl -b Delta -a 2.03 -t 2.32 -i Denmark_Delta_Omicron.csv -c -q -d 1 -f 30

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
