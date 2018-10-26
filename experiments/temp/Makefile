all:
	make clean
	make run
	make compare

run:
	python runner.py run

compare:
	../../src/flowstar f_id_dis.txt flow.txt left_after_precond 8

clean: 
	rm -rf csvs/ experiments.html f_id_dis.txt flow.txt images/ infos/ outputs/ plts/

