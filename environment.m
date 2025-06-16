function environment(base_dir)
%
% Copyright (C) 2019-2025, Dmytro Sytnyk. All rights reserved. 

if nargin==0
    prefix  = [pwd(),filesep];
else
    prefix  = [base_dir,filesep];
end
pdir = [pwd(),filesep];

addpath(prefix);

% Initialize contributed packages
cd([prefix,'contrib',filesep]);
cdir = [pwd(),filesep];
addpath(strcat(cdir,'ml'));

cd(pdir);

end
