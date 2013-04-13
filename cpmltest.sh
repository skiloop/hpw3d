#!bin/bash
#./tcpml| grep Ez | sed 's/^.*://' | tee td.txt &
#./testCPML| grep Ez | sed 's/^.*://' | tee tc.txt &

# make 
make clean && make all

# small domain size 
mkdir -p old && cd old && sh ../clean.sh
../origProgram 2>&1 > opml.txt &
cd ..

# large domain size 
mkdir -p new && cd new && sh ../clean.sh
../hpw3d 2>&1 > npml.txt &
cd ..

