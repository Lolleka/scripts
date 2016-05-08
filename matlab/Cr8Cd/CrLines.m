
function [mat_lines] = CrLines(newA)
    global A;
    %global h_lines;
    %global f_lines;
    A = newA;
    Gamma = 2.406;
    
    %carica momenti
    moments = load('.\Cr8Cd_moments.dat');

    
    h = (0:0.1:9)';
    
    %Spin totale
    ST = moments(:,1);    
    %momenti
    Sites_S0 = repmat(moments (1,2:5),size(h,1),1);
    Sites_S1 = repmat(moments (2,2:5),size(h,1),1);
    Sites_S2 = repmat(moments (3,2:5),size(h,1),1);
    Sites_S3 = repmat(moments (4,2:5),size(h,1),1);


    %calcola le frequenze di larmor
    Line_S0 = abs((Sites_S0*A+cat(2,h,h,h,h))*Gamma);
    Line_S1 = abs((Sites_S1*A+cat(2,h,h,h,h))*Gamma);
    Line_S2 = abs((Sites_S2*A+cat(2,h,h,h,h))*Gamma);
    Line_S3 = abs((Sites_S3*A+cat(2,h,h,h,h))*Gamma);
    mat_lines = [h Line_S0 h Line_S1 h Line_S2 h Line_S3];
    save('mat_lines.dat','mat_lines','-ascii');
    %figure(1)
    %plot(mat_lines(:,1),mat_lines(:,2),'o');
end



