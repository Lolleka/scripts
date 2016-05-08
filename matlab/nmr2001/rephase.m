function rw=rephase(mx,cf,nn)
%%phase correction of a FD fid with a linear phase modulation

%% input validation
if (nargin < 2),
  error('Missing arguments')
end
if ~(isnumeric(mx) & isnumeric(cf)),
  error('Inputs must be numeric')
end
if ((isvec(cf) & length(cf)~=2) | (~isvec(cf) & size(cf,1)~=2)),
  error('Input #2 must be either a 2-element vector or a 2-row matrix');
end
row_flag = isrow(mx);
if (row_flag),
  mx=mx(:);
end

if (isvec(cf)),
  cf=cf(:)*ones(1,size(mx,2));
elseif (size(cf,2) ~= size(mx,2)),
  error('Size mismatch between 1st and 2nd argument'); 
end

if (nargin < 3),
  nn=ones(1,size(mx,2))*size(mx,1);
elseif ~isnumeric(nn),
  error('Inputs must be numeric')
elseif (max(size(nn))==1), %isscal(nn),
  nn=ones(1,size(mx,2))*min(nn,size(mx,1));
elseif (prod(size(nn)) ~= size(mx,2)),
  error('Size mismatch between 1st and 3rd argument');
else
  nn=min([ones(1,size(mx,2))*size(mx,1);nn(:)']);
end

%%Do rephase
rw = zeros(size(mx));
for k=1:length(nn);
  rw(1:nn(k),k)=mx(1:nn(k),k).*exp(i * (cf(1,k) +  cf(2,k) * ...
      ((0:nn(k)-1)'-nn(k)/2)));
end
if (row_flag),
  rw=rw.';
end

%%local macros
function r=isvec(a);
r=isnumeric(a) & (min(size(a))==1);

function r=isrow(a);
r=isnumeric(a) & (size(a,1)==1);



