#!bin/bash
./dmain | grep Ez | sed 's/^.*://' | tee td.txt &
./testCPML| grep Ez | sed 's/^.*://' | tee tc.txt &

