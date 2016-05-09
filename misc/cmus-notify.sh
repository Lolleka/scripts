#!/bin/bash

# SETTINGS
NOW_PLAYING_PIPE="$HOME/tmp/now-playing"


title=""
artist=""
stat="stopped"
url=""

while [ "$1" '!=' "" ]
do
	case "$1" in
		title)
			title="$2"
		;;
		artist)
			artist="$2"
		;;
		status)
			stat="$2"
		;;
		file)
			file="$2"
		;;
		url)
			url="$2"
		;;
		*)
		;;
	esac
	shift
	shift
done

msg="Music:"
if [ "$stat" '=' 'stopped' ]
then
	msg="$msg stopped"
else
	if [ "$stat" '=' 'paused' ]
	then
		msg="$msg [paused]"
	fi

	if [ -n "$title" ]
	then
		msg="$msg \"$title\""
	else
		if [ -n "$file" ]
		then
			msg="$msg `basename "$file"`"
		else
			msg="$msg <noname>"
		fi
	fi

	if [ -n "$artist" ]
	then
		msg="$msg by $artist"
	fi
fi

notify-send -i multimedia-volume-control -t 1800 "$Cur_song"
#echo "$msg" > "$NOW_PLAYING_PIPE"