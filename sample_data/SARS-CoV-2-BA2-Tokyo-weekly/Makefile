all:RelRe.jl estimates.csv frequencies.pdf average_k.pdf

RelRe.jl:
	ln -s ../../$@ .

estimates.csv:Tokyo_BA1_BA2.csv
	julia --threads 10 RelRe.jl -b Omicron_BA1 -a 4.03 -t 0.735 -i Tokyo_BA1_BA2.csv -c -q -d 1 -f 30 -D -u

frequencies.csv:estimates.csv

loglikelihood.csv:estimates.csv

frequencies.pdf:frequencies.csv
	R -f plot_frequency_wkly.R

average_k.pdf:frequencies.csv
	R -f plot_average_k.R

clean:
	rm -f estimates.csv
	rm -f frequencies.csv
	rm -f loglikelihood.csv
	rm -f Dirichlet.csv
	rm -f frequencies.pdf
	rm -f average_k.pdf
