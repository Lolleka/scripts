#!/bin/bash

res=$(rofi -width 200 -dmenu $(i3-color-rofi) < ~/.dmenu-i3exit)

if [ $res = "logout" ]; then
    i3-msg exit
fi
if [ $res = "restart" ]; then
    reboot
fi
if [ $res = "shutdown" ]; then
    poweroff
fi
exit 0 
