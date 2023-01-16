all:RelRe.jl estimates_ml.csv

RelRe.jl:
	ln -s ../../src/$@ .

estimates_ml.csv:Tokyo_BA1_BA2.csv
	julia RelRe.jl -u W -b Omicron_BA1 -a 2.03 -t 1.392 -p 1e-12 Tokyo_BA1_BA2.csv

estimates_95CI.csv:estimates_ml.csv
	touch $@

clean:
	rm -f estimates_ml.csv
	rm -f estimates_95CI.csv
