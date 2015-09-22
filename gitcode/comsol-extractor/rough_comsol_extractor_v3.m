%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   fork of stanamodi v04 of M.Borghi, M.Mancinelli, M.Bernard from
%
%                                                         D. Bazzanella   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

close all;
clc;

fprintf('Memory Cleared\n');

%% INITIAL PARAMETERS

% setting up initial parameters
fprintf('\nSetting up initial parameters\n');

[ans1, ans2] = uigetfile({'../comsol/*.mph'});  % choose file to load
setting.comsolfile = strcat(ans2,ans1);         % filepath of the comsol model
setting.data_set = 'dset3';                     % name of the dataset in the consol model
setting.save_path = '../data/';                 % filepath of filesaves

% file-saving parameters
setting.FILE_EXT	= '.mat';         	% Determines the extension of the files: .mat or .dat
setting.FILE_TYP	= '-mat';        	% Determines the file type of the files: -mat or -ascii

% processig parameters
setting.SAVE_H		= false;        	% 1: saves also the H fields
setting.COMPUTE_H	= true;            	% 1: Compute H from matlab though impedance, 0 loads from COMsol.
setting.IMAGINARY	= false;            % Extract also imaginary part of the fields
setting.FIT_ORD     = true;             % Order of the fit in the effective index dispersion as a function of lambda
setting.TM_TOO      = true;             % classify also TM modes [1=yes]
setting.SAVE_FIELDS = true;             % 1: saves fields as matrix

% sampling parameters
setting.smpx		= 1300;             % Total samples per x direction for all the box
setting.smpy        = 1300;             % Total samples per y direction for all the box
setting.BOX_CONST   = true;             % If the box is fixed, the grid is generated only once
setting.x_cut	= 0.5;                  % um to add to the width of the waveguide core (evanescent wave)
setting.y_cut	= 0.5;                  % um to add to the height of the waveguide core (evanescent wave)

const.dB    = -30;                      % Check on the mode losses
const.c     = 299792458;            	% Speed of light [m/s]
const.T     = 293.15;                   % Ambient temperature[K]

clear ans1 ans2

%% MODEL LOADING

tic
fprintf('\nLoading the model...');
model = mphload(setting.comsolfile);
fprintf('\b\b\b:\tModel Loaded.\n');
toc

%% PARAMETERS LOADING

tic
fprintf('\nLoading model parameters...');

% get parameters
par.n_wlen	= mphglobal(model,{'n_la'},'dataset', setting.data_set,'outersolnum',1);
par.n_wlen	= par.n_wlen(1);
par.n_wg_wid	= mphglobal(model,{'n_w'},'dataset', setting.data_set,'outersolnum',1);
par.n_wg_wid	= 1;	%par.n_wg_wid(1);
par.n_wg_hgt	= mphglobal(model,{'n_h'},'dataset', setting.data_set,'outersolnum',1);
par.n_wg_hgt	= 1;	%par.n_wg_hgt(1);
par.n_temp	= mphglobal(model,{'n_T'},'dataset', setting.data_set,'outersolnum',1);
par.n_temp	= par.n_temp(1);
par.n_modi	= numel(mphglobal(model,{'emw.neff'},'dataset', setting.data_set,'outersolnum',1));

% allocate vectors
vec.wlen	= zeros(par.n_wlen,1);				% Vector for the lambdas
vec.wid		= zeros(par.n_wg_wid,1);			% Vector for the widths
vec.hgt		= zeros(par.n_wg_hgt,1);			% Vector for the heights
vec.temp    = zeros(par.n_temp,1);				% Vector for the temperatures

% get parameters
vec.wlen    = get_par(model, 'la', setting.data_set)'.*1e-6;  % [m]
vec.wid     = get_par(model, 'w', setting.data_set)';         % [um]
vec.hgt     = get_par(model, 'hg', setting.data_set)';        % [um]
vec.temp    = get_par(model, 'Temp', setting.data_set)';      % [K] - 293.15

% get materials
par.neff_si     = mphglobal(model,{'n_si'},'dataset', setting.data_set,'outersolnum',1);
par.neff_si     = par.neff_si(1);
par.neff_sio2	= mphglobal(model,{'n_sio2'},'dataset', setting.data_set,'outersolnum',1);
par.neff_sio2	= par.neff_sio2(1);

% Total number of solutions
par.n_sol	= par.n_wlen*par.n_temp*par.n_wg_wid*par.n_wg_hgt;
fprintf('\b\b\b:\tParameters Loaded.\n');
toc

%% INITIALIZE DATA STRUCTURE
tic
fprintf('\nInitialize structure for data...\n');

% Generate the structure full of NaNs, but with the right dimensions (par.n_wlen* ...)
data.neff       = zeros(par.n_wlen, par.n_temp, par.n_wg_wid, par.n_wg_hgt, 2, par.n_modi).*NaN;
data.Aeff       = zeros(par.n_wlen, par.n_temp, par.n_wg_wid, par.n_wg_hgt, 2, par.n_modi).*NaN;
%data.lossdB     = zeros(par.n_wlen, par.n_temp, par.n_wg_wid, par.n_wg_hgt, 1, par.n_modi).*NaN;

sol.wlen        = zeros(par.n_sol,par.n_modi);
sol.wid         = zeros(par.n_sol,par.n_modi);
sol.hgt         = zeros(par.n_sol,par.n_modi);
sol.temp        = zeros(par.n_sol,par.n_modi);
sol.neff        = zeros(par.n_sol,par.n_modi);
sol.Aeff        = zeros(par.n_sol,par.n_modi);
sol.loss_dB     = zeros(par.n_sol,par.n_modi);
toc

% fill in some data
fprintf('\nFill initial data from the %d solutions...\n', par.n_sol);
tic
for jj=1:par.n_sol
    fprintf('%d ',jj);
    sol.wlen(jj,:)  = mphglobal(model,{'la'},'dataset', setting.data_set,'outersolnum',jj).*1e-6; % [m]
    sol.wid(jj,:)   = mphglobal(model,{'w'}, 'dataset', setting.data_set,'outersolnum',jj);       % [um]
    sol.hgt(jj,:)   = mphglobal(model,{'hg'},'dataset', setting.data_set,'outersolnum',jj);       % [um]
    sol.temp(jj,:)  = mphglobal(model,{'Temp'},'dataset', setting.data_set,'outersolnum',jj);     % [K]-293.15
    sol.neff(jj,:)  = mphglobal(model,{'emw.neff'},'dataset', setting.data_set,'outersolnum',jj);
    % save effective area
    sol.Aeff(jj,:)  = mphglobal(model,{'Aeff'},'dataset', setting.data_set,'outersolnum',jj);
end
clear tmp jj
fprintf('\n');
toc

%% SAMPLING SOLUTIONS
tic
fprintf('\nFilling the data structures...\n');
for jj=1:par.n_sol

    index.wlen    = find( vec.wlen == sol.wlen(jj) );
    index.wid     = find( vec.wid  == sol.wid(jj)  );
    index.hgt     = find( vec.hgt  == sol.hgt(jj)  );
    index.temp    = find( vec.temp == sol.temp(jj) );
    
    TEc = fliplr( real( sol.neff(jj,:) ) );
    TMc = TEc(8:end);
    TEc = TEc(1:7);
	TMc(find(TMc) < 1.9) = [];
    
	data.neff(index.wlen, index.temp, index.wid, index.hgt, 1, :)     = TEc;
    data.neff(index.wlen, index.temp, index.wid, index.hgt, 2, :)     = TMc;
    
    TEc = fliplr( sol.Aeff(jj,:) );
    TMc = TEc(8:end);
    TEc = TEc(1:7);
	TMc(find(TMc) < 1.9) = [];
    
	data.neff(index.wlen, index.temp, index.wid, index.hgt, 1, :)     = TEc;
    data.neff(index.wlen, index.temp, index.wid, index.hgt, 2, :)     = TMc;

end
toc

%% SAVING

tic
fp = datestr(datetime,'yyyy-mm-dd_HH:MM:ss');
fp = strcat(setting.save_path, fp, setting.FILE_EXT);

fprintf('Saving data in %s\n',fp);
save(fp, 'data', 'vec', 'par');

clear fp ans
toc

%% END
fprintf('Data saved.\nProgram finished.');
