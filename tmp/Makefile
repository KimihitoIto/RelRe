all:estimates_ml.csv

estimates_ml.csv:Tokyo_BA1_BA2.csv
	julia advantage-SI-Rt-wkly.jl

estimates_95CI.csv:estimates_ml.csv
	touch $@

clean:
	rm -f estimates_95CI.csv
	rm -f estimates_ml.csv
