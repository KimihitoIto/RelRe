all:RelRe.jl estimates.csv

RelRe.jl:
	ln -s ../../src/$@ .

estimates.csv:Tokyo_BA1_BA2.csv
	julia --threads 10 RelRe.jl -u W -b Omicron_BA1 -a 2.03 -t 1.392 -p 1e-4 -w 2022-03-01 2022-04-01 -c Tokyo_BA1_BA2.csv -i

clean:
	rm -f estimates.csv
