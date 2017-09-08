compare_cylinder:
	python benchmark_sym.py
	head -n 12 c.log
	# grep matrices c.log
	head -n 12 symc.log
	# grep matrices symc.log

diff:
	diff results.csv sym_results.csv

clean:
	rm -f c.log symc.log
