#!bin/bash
./tcpml | tee tlarge.txt &
./testCPML| tee tsmall.txt &
#./tcpml | grep Ez | sed 's/^.*://' | tee tlarge.txt &
#./testCPML| grep Ez | sed 's/^.*://' | tee tsmall.txt &


