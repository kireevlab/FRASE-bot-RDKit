#!/bin/bash

grep '::' *.out > if.txt
awk -F ' :: ' '{print$1}' if.txt  | awk -F ':' '{print$2}' > label
sed -i 's/$/.pdb/g' label
cp label label1
chmod +x label
sed -i -e 's/^/mv /g' -e 's/$/ target_ML/g' label

cd newFragments && mkdir target_ML
../label
#rm * 
mv target_ML .. 


cd .. && awk -F ' :: ' '{print$2}' if.txt  > IF_target
paste label1 IF_target > if1.txt
cat -n if1.txt > IF_target_sort.txt

sed -i -e 's/^/    /' -e 's/$/,/' -e '$ s/.$//' -e '1i [' -e '$a ]'  IF_target && mv IF_target IF_target.json

mkdir slurmout && mv slurm-*.out slurmout
rm if.txt if1.txt label* -r triplet_db
