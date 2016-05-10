#!/bin/zsh

BASE_BCK_DIR="$HOME/git/dotfiles"
mkdir -p $BASE_BCK_DIR/home/.config
mkdir -p $BASE_BCK_DIR/misc
mkdir -p $BASE_BCK_DIR/etc

echo "Base directory: $BASE_BCK_DIR"

# --------- /etc files --------- #
#Backup ssh
if [ -d "/etc/ssh" ]
then
    echo "Saving ssh config files..."
    sudo cp -r /etc/ssh $BASE_BCK_DIR/etc
else
    echo "ssh not found."
fi

#Backup Samba config
if [ -e "/etc/samba/smb.conf" ]
then
    echo "Saving Samba config files..."
    mkdir -p $BASE_BCK_DIR/etc/samba
    cp /etc/samba/smb.conf $BASE_BCK_DIR/etc/samba
else
    echo "Samba not found."
fi

#Backup bash.bashrc
if [ -e "/etc/bash.bashrc" ]
then
    echo "Saving bash.bashrc..."
    cp /etc/bash.bashrc $BASE_BCK_DIR/etc
else
    echo "bash.bashrc not found."
fi

#Backup inadyn.conf
if [ -e "/etc/inadyn.conf" ]
then
    echo "Saving inadyn.conf..."
    cp /etc/inadyn.conf $BASE_BCK_DIR/etc
else
    echo "inadyn.conf not found."
fi

#Backup iptables.rules
if [ -e "/etc/iptables/iptables.rules" ]
then
    echo "Saving iptables rules..."
    mkdir -p $BASE_BCK_DIR/etc/iptables
    cp /etc/iptables/iptables.rules $BASE_BCK_DIR/etc/iptables
else
    echo "iptables rules not found."
fi

#Backup minidlna
if [ -e "/etc/minidlna.conf" ]
then
    echo "Saving minidlna.conf..."
    cp /etc/minidlna.conf $BASE_BCK_DIR/etc/minidlna.conf
else
    echo "minidlna.conf not found."
fi

#Backup sudoers
if [ -e "/etc/sudoers" ]
then
    echo "Saving sudoers..."
    sudo cp /etc/sudoers $BASE_BCK_DIR/etc
else
    echo "sudoers not found."
fi

#Backup fstab
if [ -e "/etc/fstab" ]
then
    echo "Saving fstab..."
    cp /etc/fstab $BASE_BCK_DIR/etc/fstab
else
    echo "fstab not found."
fi

# --------- home files --------- #
#Backup ssh
if [ -d "$HOME/.ssh" ]
then
    echo "Saving ssh config files..."
    cp -r $HOME/.ssh $BASE_BCK_DIR/home
else
    echo "ssh not found in home."
fi
#Backup i3 dotfiles
if [ -d "$HOME/.i3" ]
then
    echo "Saving i3 config files..."
    cp -r $HOME/.i3 $BASE_BCK_DIR/home
else
    echo "i3 not found."
fi
#Backup i3blocks dotfiles
if [ -d "$HOME/.config/i3blocks" ]
then
    echo "Saving i3blocks files..."
    cp -r $HOME/.config/i3blocks $BASE_BCK_DIR/home/.config
    mkdir -p $BASE_BCK_DIR/usr/lib/i3blocks
    cp -r /usr/lib/i3blocks/* $BASE_BCK_DIR/usr/lib/i3blocks
else
    echo "i3blocks not found."
fi
#Backup transmission config
if [ -e "$HOME/.config/transmission/settings.json" ]
then 
    echo "Saving transmission config files..."
    mkdir -p $BASE_BCK_DIR/home/.config/transmission
    cp $HOME/.config/transmission/settings.json $BASE_BCK_DIR/home/.config/transmission
else
    echo "transmission not found."
fi

#Backup .Xdefaults
if [ -e "$HOME/.Xdefaults" ]
then
    echo "Saving .Xdefaults..."
    cp $HOME/.Xdefaults $BASE_BCK_DIR/home
else
    echo ".Xdefaults not found."
fi

#Backup .zshrc
if [ -e "$HOME/.zshrc" ]
then
    echo "Saving .zshrc..."
    cp $HOME/.zshrc $BASE_BCK_DIR/home
else
    echo ".zshrc not found."
fi

#Backup .zprofile
if [ -e "$HOME/.zprofile" ]
then
    echo "Saving .zprofile..."
    cp $HOME/.zprofile $BASE_BCK_DIR/home
else
    echo ".zprofile not found."
fi

#Backup .lscolors
if [ -e "$HOME/.lscolors" ]
then
    echo "Saving .lscolors..."
    cp $HOME/.lscolors $BASE_BCK_DIR/home
else
    echo ".lscolors not found."
fi

#Backup .vimrc
if [ -e "$HOME/.vimrc" ]
then
    echo "Saving .vimrc..."
    cp $HOME/.vimrc $BASE_BCK_DIR/home
else
    echo ".vimrc not found"
fi

#Backup .vim
#sudo cp -r $HOME/.vim $BASE_BCK_DIR/home/

#Backup ranger config
if [ -d "$HOME/.config/ranger" ]
then
    echo "Saving ranger config files..."
    mkdir -p $BASE_BCK_DIR/home/.config/ranger
    cp -r $HOME/.config/ranger/* $BASE_BCK_DIR/home/.config/ranger
else
    echo "ranger config files not found."
fi

#Backup ncmpcpp
if [ -e "$HOME/.ncmpcpp" ]
then
    echo "Saving .ncmpcpp..."
    mkdir -p $BASE_BCK_DIR/home/.ncmpcpp
    cp -r $HOME/.ncmpcpp/* $BASE_BCK_DIR/home/.ncmpcpp
else
    echo ".ncmpcpp not found."
fi

#Backup .xinitrc
if [ -e "$HOME/.xinitrc" ]
then
    echo "Saving .xinitrc..."
    cp $HOME/.xinitrc $BASE_BCK_DIR/home/.xinitrc
else
    echo ".xinitrc not found."
fi

#Backup mldonkey
if [ -d "$HOME/.mldonkey" ]
then
    echo "Saving mldonkey config files..."
    mkdir -p $BASE_BCK_DIR/home/.mldonkey
    cp -r $HOME/.mldonkey/* $BASE_BCK_DIR/home/.mldonkey
    echo "mldonkey config files not found."
fi

#Backup mplayer
if [ -d "$HOME/.mplayer" ]
then
    echo "Saving mplayer config files..."
    mkdir -p $BASE_BCK_DIR/home/.mplayer 
    cp -r $HOME/.mplayer/* $BASE_BCK_DIR/home/.mplayer 
else
    echo "mplayer config files not found."
fi

#Backup .gnupg
if [ -d "$HOME/.gnupg" ]
then
    echo "Saving gnupg config files..."
    mkdir -p $BASE_BCK_DIR/home/.gnupg
    cp -r $HOME/.gnupg/* $BASE_BCK_DIR/home/.gnupg
else
    echo "gnupg config files not found."
fi

#Backup .vnc
if [ -d "$HOME/.vnc" ]
then
    echo "Saving vnc config files..."
    mkdir -p $BASE_BCK_DIR/home/.vnc
    cp -r $HOME/.vnc/* $BASE_BCK_DIR/home/.vnc
else
    echo "vnc config files not found."
fi
# --------- misc files ---------- #

#Backup crontab
crontab -l > $BASE_BCK_DIR/misc/cron.bck
#Backup self
cp $0 $BASE_BCK_DIR/misc

sudo zip -r -q $BASE_BCK_DIR/dotfiles $BASE_BCK_DIR/* 

