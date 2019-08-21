#!/bin/bash
pid=$1 
echo $pid
interval=100  
while true
do
    echo $(date +"%y-%m-%d %H:%M:%S") >> $2
    cat  /proc/$pid/status|grep -e VmRSS >> $2
    sleep $interval
done
