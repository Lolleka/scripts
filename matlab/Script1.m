load matrix.dat
[row,col]=size(matrix);
x=matrix(1,2:col);
y=matrix(2:row,1);
z=matrix(2:row,2:col);
figure(1);
imagesc(x,y,z);
d=[];
row=row-1;
col=col-1;
for i=1:row
    u=[0 diff(z(i,1:col))./diff(x)];
    d=[d; u];
end
figure(2);
imagesc(x,y,d);