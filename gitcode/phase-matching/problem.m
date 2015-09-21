%% problem

% carico quello simmetrico
[ans1, ans2] = uigetfile({'../data/*.mat'});     % choose file to load
load( strcat(ans2, ans1) );
sym     = data;
fprintf(['File ', ans1, ' loaded.\n']);

% carico quello asimmetrico
[ans1, ans2] = uigetfile({'../data/*.mat'});     % choose file to load
load( strcat(ans2, ans1) );
symA    = data;
fprintf(['File ', ans1, ' loaded.\n']);

% carico quello bisimmetrico
[ans1, ans2] = uigetfile({'../data/*.mat'});     % choose file to load
load( strcat(ans2, ans1) );
symB    = data;
fprintf(['File ', ans1, ' loaded.\n']);

clear data ans1 ans2

% create matrix @ T_amb for TE only, all orders, all
Msym = squeeze(sym.neff(:,1,1,1,1,:));
MsymA = squeeze(symA.neff(:,1,1,1,1,:));
MsymB = squeeze(symB.neff(:,1,1,1,1,:));

% sym vs asym
SvsA = Msym - MsymA;
% sym vs bsym
SvsB = Msym - MsymB;
% asym vs bsym
AvsB = MsymA - MsymB;