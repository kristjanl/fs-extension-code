This directory contains the script for running experiments for comparing the compositional approach to 
the non-compositional approach on hybrid systems which had some form of compositionality on them naturally.

Running CFlow can be done with the command  
```
python comps.py run
```
which produces CSV data files in the directory `csvs`.

Using the CSV data files, HTML file depicting the flowpipes can be generated with the command
```
python comps.py compare
```
which produces `experiments_comps.html`.

Both of these commands can be issued together with the command
```
python comps.py both
```