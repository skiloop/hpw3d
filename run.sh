#!bin/bash
./dmain | grep Ez | sed 's/^.*://' | tee td.txt &

