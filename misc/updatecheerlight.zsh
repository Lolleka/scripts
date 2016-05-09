#!/usr/bin/zsh
curl http://api.thingspeak.com/channels/1417/field/1/last.txt > /tmp/cheerlight
CHEERLIGHT=$(cat /tmp/cheerlight | cut -d "%" -f 1)
