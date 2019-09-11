This directory contains the script for running experiments for different processing methods 
and comparing the compositional approach to the non-compositional approach (on the systems which 
are listed in the thesis).

For processing methods, running CFlow can be done with the command  
```
python method.py run
```
which produces CSV data files in the directory `csvs`.

Using the CSV data files, HTML file depicting the flowpipes can be generated with the command
```
python method.py compare
```
which produces `experiments_method.html`.

Both of these commands can be issued together with the command
```
python method.py both
```

For comparing the compositional approach with the non-compositional one, running CFlow can be done 
with the command  
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