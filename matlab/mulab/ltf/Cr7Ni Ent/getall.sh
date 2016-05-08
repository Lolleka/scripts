#!/bin/sh
 rsync -avR l_musr_ltf@pc9627:'`find /afs/.psi.ch/project/bulkmusr/data/'$1'/d2013/tdc/ -mtime -'$2' -type f`' .
 mv ./afs/.psi.ch/project/bulkmusr/data/$1/d2013/tdc/* ./
 rm -r ./afs