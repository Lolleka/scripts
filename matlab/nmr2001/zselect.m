%ZSELECT: select a number of 2-D zones (experiments) from workspace, delete the other ones. 
%USAGE:
%  zselect(indices)

function zselect(ix);
zkill(ix,1);