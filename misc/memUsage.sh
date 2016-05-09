#!/bin/bash
while [ 1 ]
do
   total=`cat /proc/meminfo | grep MemTotal | cut -d ":" -f 2 | cut -d "k" -f 1 | tr -d " "`
   free=`cat /proc/meminfo | grep MemFree | cut -d ":" -f 2 | cut -d "k" -f 1 | tr -d " "`
   currentTime=`date +%s`
   echo -n $currentTime
   echo -n " "
   echo "$total - $free" | bc
   sleep 1
done
