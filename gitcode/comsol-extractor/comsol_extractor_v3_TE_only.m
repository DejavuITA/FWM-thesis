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
data.neff       = zeros(par.n_wlen, par.n_temp, par.n_wg_wid, par.n_wg_hgt, 1, par.n_modi).*NaN;
data.Aeff       = zeros(par.n_wlen, par.n_temp, par.n_wg_wid, par.n_wg_hgt, 1, par.n_modi).*NaN;
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
	
	data.neff(index.wlen, index.temp, index.wid, index.hgt, :)     = fliplr( real( sol.neff(jj,:) ) );
    data.Aeff(index.wlen, index.temp, index.wid, index.hgt, :)     = fliplr( sol.Aeff(jj,:) );
	%data.lossdB(index.wlen, index.temp, index.wid, index.hgt, :)   = fliplr( sol.loss_dB(jj,:) );

end
toc
%{

tic
fprintf('\nFilling the data structures...\n');

% it generates coordinates grid (x & y) of the box
% and coordinates grid (x2 & y2) of the waveguide core
par.ord_max = 0;
par.neff_max   = 1.4;

% take the measures of the box
if (setting.BOX_CONST)
    par.box_w   = mphglobal(model,{'wt'},'dataset', setting.data_set,'outersolnum',1);    % [um]
    par.box_w   = par.box_w(1);
    par.box_ho  = mphglobal(model,{'ho'},'dataset', setting.data_set,'outersolnum',1);    % [um]
    par.box_ho  = par.box_ho(1);
    par.box_hts = mphglobal(model,{'hts'},'dataset', setting.data_set,'outersolnum',1);   % [um]
    par.box_hts = par.box_hts(1);
    par.box_h   = par.box_ho+par.box_hts;
end

% mode index
index.progress_solution = 0;
index.max_n_modi        = 9;%par.n_modi+1;
for jj=1:par.n_sol
   
    par.step_x   = sol.wid(jj)/setting.smpx;            % Horizontal step is box_w / sampling (um/step)   
    par.step_y   = par.box_h/setting.smpy;              % Vertical   step is box_h / sampling (um/step)
    sample.x = linspace(-sol.wid(jj)/2,sol.wid(jj)/2,setting.smpx);
    sample.ddd = [sample.x;zeros(1,length(sample.x))];    
    
    % calculate dB losses
    sol.loss_dB(jj,:) = -20*log10( exp(1) ).*abs( imag( sol.neff(jj,:) ) )*2*pi./( sol.wlen(jj,:)*1e2 );
    
    % get indexes in the ordered arrays of the current data
    % find(A==b) gets index 'j' in A where A(j)==b
    index.wlen    = find( vec.wlen == sol.wlen(jj) );
    index.wid     = find( vec.wid  == sol.wid(jj)  );
    index.hgt     = find( vec.hgt  == sol.hgt(jj)  );
    index.temp    = find( vec.temp == sol.temp(jj) );
    
    fprintf('Sampling and loading solution %3d of %3d. %3.1d%% ...\n', jj, par.n_sol,round(jj/par.n_sol*100));
    
    % total amount of points where to evaluate the field(s)
    sample.x_tot = linspace(-sol.wid(jj)/2-setting.x_cut,sol.wid(jj)/2+setting.x_cut,setting.smpx);
    sample.y_tot = linspace(-sol.hgt(jj)/2-setting.y_cut,sol.hgt(jj)/2+setting.y_cut,setting.smpx);    
    sample.points = zeros(2,length(sample.y_tot)*length(sample.x_tot));
    for kk=1:length(sample.y_tot)
        for ll=1:length(sample.x_tot)
            sample.points(:, (kk-1)*length(sample.y_tot)+ll) = [sample.x_tot(ll) sample.y_tot(kk)];
        end
    end % si potrebbe fare con mesh grid

    % load field matrix in all model
    [sample.Exfield,sample.d2TE] = mphinterp(model,{'real(emw.Ex)'},'coord',sample.ddd,'dataset', setting.data_set,'outersolnum',jj);
    [sample.Eyfield,sample.d2TM] = mphinterp(model,{'real(emw.Ey)'},'coord',sample.ddd,'dataset', setting.data_set,'outersolnum',jj);
%   [Ezfield,d2TM] = mphinterp(model,{'real(emw.Ez)'},'coord',ddd,'dataset', data_set,'outersolnum',jj);

    % for cycle on the orders of the modes
    sample.ratio = zeros(par.n_modi,1);
    for kk=1:par.n_modi
        fprintf('Sampling and loading mode %3d of %3d ...\n', kk, par.n_modi);
        
        % check if the imaginary part of n_eff is too big:
        % in that case it skips
        if( sol.loss_dB(jj,kk) > const.dB )

            % rapp_controllo_livello
            sample.ratio(kk) = max(abs(sample.Exfield(kk,:)))/max(abs(sample.Eyfield(kk,:))); %Rapporto tra i massimi globali dello stesso modo
            index.progress_solution = index.progress_solution + 1;
            if ( sample.ratio(kk)>1 )
                
                % classify TE modes

                sample.npksTE = size( fnzeros( spline( sample.ddd(1,:), sample.Exfield(kk,:) ) ), 2 );
                fprintf('\tsampling a TE mode; npksTE = %d\n',sample.npksTE);

                if( (sample.npksTE + 1)<index.max_n_modi && (sample.npksTE+1)>0 )
                    
                    % SE abbiamo dei picchi, e sono meno dei modi... ?
                    data.neff(index.wlen, index.temp, index.wid, index.hgt, 1, sample.npksTE+1)     = real( sol.neff(jj,kk) );
                    data.Aeff(index.wlen, index.temp, index.wid, index.hgt, 1, sample.npksTE+1)     = sol.Aeff(jj,kk);
                    data.lossdB(index.wlen, index.temp, index.wid, index.hgt, 1, sample.npksTE+1)   = sol.loss_dB(jj,kk);

%?                    index(index.progress_solution,:)=[sol.hgt(jj),sol.wid(jj),kk,sample.npksTE,1,sol.wlen(jj)];
%?                    index2(index.progress_solution,:)=[jj,kk,1,sample.npksTE];

                    % maxi effective index
                    if( par.neff_max < sol.neff(jj,kk) ) 
                        par.neff_max = sol.neff(jj,kk); 
                    end
                    % max order found
                    if( par.ord_max < sample.npksTE+1 )
                        par.ord_max = sample.npksTE+1;
                    end
                end
            else
                % classify TM modes
                
                sample.npksTM = size( fnzeros( spline( sample.ddd(1,:), sample.Eyfield(kk,:) ) ), 2 );
                fprintf('\tsampling a TM mode; npksTM = %d\n',sample.npksTM);
                
                if( sample.npksTM+1<index.max_n_modi && sample.npksTM+1>0 && setting.TM_TOO)
                    data.neff(index.wlen, index.temp, index.wid, index.hgt, 2, sample.npksTM+1)     = real( sol.neff(jj,kk) );
                    data.Aeff(index.wlen, index.temp, index.wid, index.hgt, 2, sample.npksTM+1)     = sol.Aeff(jj,kk);
                    data.lossdB(index.wlen, index.temp, index.wid, index.hgt, 2, sample.npksTM+1)   = sol.loss_dB(jj,kk);
                    
%?                    FTM(:,pos_lam,npksTM+1)=Eyfield(kk,:);
%?                    index(index.progress_solution,:)=[sol.hgt(jj),sol.wid(jj),kk,sample.npksTM,-1,sol.wlen(jj)];
%?                    index2(index.progress_solution,:)=[jj,kk,-1,sol.npksTM];

                    % maxi effective index
                    if( par.neff_max < sol.neff(jj,kk) )
                        par.neff_max = sol.neff(jj,kk);
                    end
                    % max order found
                    if( par.ord_max < sample.npksTM+1 )
                        par.ord_max = sample.npksTM+1;
                    end
                    
                end
            end
        end
    end % end of for cycle on modes (par.n_modi)
end % end of for cycle on solutions (par.n_sol)
index.progress_solution = 0;

% save unit of measure
units.lossdBunit='dB/cm';
units.neffunit='Effective index';

clear jj kk ll index ans sample
toc

%}

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
