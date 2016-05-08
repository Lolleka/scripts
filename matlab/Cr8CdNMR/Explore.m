function [ output_args ] = Explore( pars )
fn = pars{1};
data = load(fn);
figure(1);
data = data(find(data(:,4) == 1),:);
data = data(find(data(:,3)>=0),:);
plot(data(:,1),data(:,3),'o');
hold on;
plot(data(:,1),data(:,3),'-');
hold off;
setup;
pars{1} = strrep(pars{1},'.dat','.pna');
figure(2);
parmaecho (pars{:});
end

