#!bin/bash
#./tcpml| grep Ez | sed 's/^.*://' | tee td.txt &
#./testCPML| grep Ez | sed 's/^.*://' | tee tc.txt &

# make 
make clean && make all

# small domain size 
mkdir -p old && cd old && sh ../clean.sh
../origProgram | grep Ez | sed -s 's/^.*://'  > opml.txt &
cd ..

# large domain size 
mkdir -p new && cd new && sh ../clean.sh
../hpw3d | grep Ez | sed -s 's/^.*://' > npml.txt &
cd ..

