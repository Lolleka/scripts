USER="lorenzo.bordonali@gmail.com"
PASS="\$avulleh\$"
 
COUNT=`curl -su $USER:$PASS https://mail.google.com/mail/feed/atom || echo "<fullcount>unknown number of</fullcount>"`
COUNT=`echo "$COUNT" | grep "<fullcount>" | sed -e s/"<fullcount>"// -e s/"<\/fullcount>"//`
echo $COUNT
if [ "$COUNT" != "0" ]; then
    if [ "$COUNT" = "1" ];then
        WORD="mail";
    else
        WORD="mails";
    fi
fi