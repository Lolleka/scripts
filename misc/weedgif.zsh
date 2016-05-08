#!/usr/bin/zsh
BASEDIR="/mnt/nasstorage/data/stills/$1" #$(date +'%d.%m.%y')"
echo $BASEDIR
for img in $BASEDIR/*.jpg ; do
   echo "Converting $img ..."
   imgname=$(basename $img | cut -d"-" -f 1,2)
   convert                  \
     -background '#0008'    \
     -gravity center        \
     -fill white            \
     -format gif            \
     -pointsize 40          \
     -geometry +10+10       \
      caption:"${imgname}"      \
      "${img}"              \
     +swap                  \
     -gravity NorthWest     \
     -composite        \
      "/tmp/wc-${imgname}.giff"

   #convert                  \
   #  -format gif            \
   #   "${img}"              \
   #   "/tmp/wc-${imgname}.giff"

done
echo "Creating animated gif..."
convert -monitor -loop 0 -delay 15 /tmp/*.giff /tmp/weed.gif
rm /tmp/*.giff
mv /tmp/weed.gif $BASEDIR
