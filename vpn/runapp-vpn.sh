#!/usr/bin/env zsh

# Initialize VPN
sudo vpnns up
sudo vpnns start_vpn

# Popcorn time!
sudo ip netns exec frootvpn sudo -u $USER $1

# Cleanup
sudo ip netns pids frootvpn | xargs -rd'\n' sudo kill
sudo vpnns down
