#!bin/bash
#./tcpml| grep Ez | sed 's/^.*://' | tee td.txt &
#./testCPML| grep Ez | sed 's/^.*://' | tee tc.txt &

# make 
make clean && make all

# small domain size 
mkdir small && cd small
mkdir small && cd small
../testCPML 2>&1 > outl.txt &
cd ..

# large domain size 
mkdir large && cd large
mkdir large && cd large
../tcpml 2>&1 > outs.txt &
cd ..

