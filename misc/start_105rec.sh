#!/bin/bash
cd ~
streamripper http://shoutcast.unitedradio.it:1101 -a "./down/Zoo105/Zoo105-$(date +%m-%d-%y).mp3" -l $1 -s