%Loads sample and ladder curve and extract fragment lenght distribution
function [res] = BALoad(Sample_FN, Ladder_FN)

%Load Ladder first
Ladder = LoadFile(Ladder_FN);

%find marker peaks and tag them with the appropriate length value (in bp)
[pks, locs] = findpeaks (Ladder(:,2), 'MINPEAKHEIGHT', max(Ladder(:,2))/5, 'NPEAKS', 15, 'MINPEAKDISTANCE',20);
%marker length values (note: 0 means untagged)

figure(2);
plot(Ladder(:,1),Ladder(:,2));
hold on;
plot(Ladder(locs,1),pks,'o','Color','red');
hold off;

bp = [35 0 100 150 200 300 400 500 600 0 1000 2000 0 0 10380];

%Find and remove untagged marker peaks
l = find (bp > 0);
bp = bp(l);
locs = locs(l);

%find times corresponding to each marker
times = Ladder(locs,1);

%plot (times, bp,'-o');
%Load Sample
Sample = LoadFile(Sample_FN);
Sample_min = min(Sample(:,2));
if Sample_min < 0
    Sample(:,2) = Sample(:,2) - Sample_min ;
end

bp_tot = interp1(times,bp,Sample(:,1));
amp_tot = Sample(:,2);

%build function return value 
res = [bp_tot amp_tot];

%trim curve to exclude marker tails ([40 9600] looks like a good length
%range (in bp))
locs = find(res(:,1) > 45 & res(:,1) < 9000);
res = res(locs,:);
end

function [data] = LoadFile(filename)
    %Open file
    fid = fopen(filename,'rt');
    %Skip header and parse file
    datacell = textscan(fid, '%f%f', 'Delimiter', ',', 'HeaderLines', 18, 'CollectOutput', 1);
    %Close file
    fclose(fid);
    %Return data
    data = datacell{1};
    %Cleanup
    clear datacell;
end