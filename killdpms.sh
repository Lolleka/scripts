#!/bin/bash
XAUTHORITY=/home/lolleka/.Xauthority DISPLAY=:0 xset dpms 0 0 0
XAUTHORITY=/home/lolleka/.Xauthority DISPLAY=:0 xset s noblank
sleep 1
XAUTHORITY=/home/lolleka/.Xauthority DISPLAY=:0 xset -dpms
XAUTHORITY=/home/lolleka/.Xauthority DISPLAY=:0 xset q > /home/lolleka/scripts/killdpms.log
date >> /home/lolleka/scripts/killdpms.timestamp
#xset s noblank
#sleep 2
#xset -dpms
#xset s noblank
#sleep 3
#xset -dpms
#xset s noblank
#sleep 4
#xset -dpms
#xset s noblank
#sleep 5
#xset -dpms
#xset s noblank