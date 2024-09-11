all:RelRe_simul.jl counts.csv frequencies.pdf	

RelRe_simul.jl:
	ln -s ../../$@ .

counts.csv:parameters.csv
	julia RelRe_simul.jl -b v001 -a 4.109 -t 0.73 -i $< -d 1 -e 2024-06-30 

frequencies.pdf:frequencies.csv
	R -f plot_frequency.R

clean:
	rm -f counts.csv
	rm -f *pdf 
