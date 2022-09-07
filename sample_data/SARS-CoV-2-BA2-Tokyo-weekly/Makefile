all:RelRe.jl estimates_ml.csv

RelRe.jl:
	ln -s ../../src/$@ .

2022-04-19estimates_ml.csv:Tokyo_BA1_BA2.csv
	julia RelRe.jl -c Tokyo_BA1_BA2.csv -u W -b Omicron_BA1 -a 2.03 -t 1.392 -e 2022-04-19

2022-04-19estimates_95CI.csv:estimates_ml.csv
	touch $@

clean:
	rm -f 2022-04-19estimates_ml.csv
	rm -f 2022-04-19estimates_95CI.csv
