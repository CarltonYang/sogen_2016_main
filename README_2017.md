
********************
    This is a readme file for 2017 segmentation clock project
    It contains information about how to complie, run simulation and use plotting scripts

0. Copy everthing onto cluster

1. To compile:
1.1
--- simulation:
cd into simulation folder
call scons
1.2
--- sres:
cd into sres folder
for parallel sres: scons mpi=1


2. To run:
on CPU:
./simulation -x 50 -w 10 -y 4 -i 2017refine2.params -s 1000 -S 0.01 -V 600 -Y 600 -Z 600 -Q 600 -K 600 -r gradient1.txt -u pert.txt -t -D output -L 1

call this line for a normal simulation
for more details, refer to the general readme file.



on Cluster:
use qsub sres.run (sres.run - sres5.run)
change test name, output file name and ranges file accordingly



3. 

3.1
Use extract_para_set.py to extraction all parameter sets from the output above.
example:
python extract_para_sets.py ./s51.txt ./ test1 0.0075
details can be found in extract_para_sets.py

3.2
****
NOTE: need to rm features.txt under simulation folder everytime before you run this step.
****
For a certain parameter, run it using code from 2. on CPU on cluster (Mac and cluster have different outputs)
To run multiple sets, use -p flag
(remember change .params to the output file in 3.1)
Outputs will be written in features.txt under the same directory as simulation


NOTE: need to rm features.txt under simulation folder everytime before you run this step.


3.3
Use scp to copy features.txt from cluster to local 

Convert features.txt to features.csv using the following command line:

python feature_to_csv.py features.txt

New features.csv will be written in the same directory







4. 

Use scripts to plot



For 4.1 and 4.2,
Make sure wanted concentration levels are output to file (to change species, change code in mutant create in init.cpp)
——————————————————————————————————————————————
4.1
plot mh1 synchronization 

python plot-snap-sync.py ../simulation/output/wildtype/set_0.cons ./check wt_sync_mh1 0.01 600
——————————————————————————————————————————————
4.2
plot mh1 vs. mespa/mespb complementary

python plot-snap-comp.py ../simulation/output/wildtype/set_0.cons ../simulation/output/MESPAOVER/set_0.cons ./check compl-mh1-mespa 0.01 600
python plot-snap-comp.py ../simulation/output/wildtype/set_0.cons ../simulation/output/MESPBOVER/set_0.cons ./check compl-mh1-mespb 0.01 600





——————————————————————————————————————————————

NOTE:



——————————————————————————————————————————————
4.3
plot mh1 synchronization bar graph

python plot-bar-sync_mh1.py ./check sync_score_mh1
——————————————————————————————————————————————
4.4
plot mesp synchronization bar graph

python plot-bar-sync_mesp.py ./check sync_score_mesp
——————————————————————————————————————————————
4.5
plot her1 mRNA amplitude (WT/Delta, WT/DAPT) bar graph

python plot-bar-amp_mh1-5 ./check amp_score_mh1
——————————————————————————————————————————————
4.6
plot mespa amplitude comparison (WT/Delta, WT/DAPT) bar graph

python plot-bar-amp-mmespa-5 ./check amp_score_mespa
——————————————————————————————————————————————
4.7
plot mespb mRNA in Mespa/Mespb amplitude

python plot-bar-amp-mesp-1.py ./check amp_score_mesp-1
——————————————————————————————————————————————
4.8
plot wildtype/her1over mRNA amplitude

python plot-bar-amp-her1over-1.py ./ testing features.csv






