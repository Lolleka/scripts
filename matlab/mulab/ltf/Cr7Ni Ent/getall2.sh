#!/bin/sh
 for ((i=$3; i<=$4; i++))
 do
 	echo "/afs/.psi.ch/project/bulkmusr/data/$1/d$2/tdc/deltat_tdc_$1_$i.bin" >> list.dat
 done
 rsync -avR --files-from=./list.dat l_musr_ltf@pc9627:/ .
 rm list.dat
 mv ./afs/.psi.ch/project/bulkmusr/data/$1/d2013/tdc/* ./
 rm -r ./afs
