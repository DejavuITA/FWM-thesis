%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   refacing of stanamodi v04 from
%
%                                         M.Borghi, M.Mancinelli, M.Bernard
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

close all;
clear all;
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
vec.wlen    = get_par(model, 'la', setting.data_set)'.*1e-6;                     % 'la' would be in [um]
vec.wid     = get_par(model, 'w', setting.data_set)';
vec.hgt     = get_par(model, 'hg', setting.data_set)';                      % {'hn'} ?
vec.temp    = get_par(model, 'Temp', setting.data_set)';

% get materials
par.neff_si     = mphglobal(model,{'n_si'},'dataset', setting.data_set,'outersolnum',1);
par.neff_si     = par.neff_si(1);
par.neff_sio2	= mphglobal(model,{'n_sio2'},'dataset', setting.data_set,'outersolnum',1);
par.neff_sio2	= par.neff_sio2(1);

% Indices for the parameters sweep
%index.wlen	= 1;
%index.wid   = 1;
%index.hgt	= 1;
%index.temp	= 1;

% Total number of solutions
par.n_sol	= par.n_wlen*par.n_wg_wid*par.n_wg_hgt*par.n_temp;
fprintf('\b\b\b:\tParameters Loaded.\n');
toc

%% INITIALIZE DATA STRUCTURE
tic
fprintf('\nInitialize structure for data...\n');

% Generate the structure full of zeros, but with the right dimensions (par.n_wlen)
data.neff       = zeros(par.n_wlen, par.n_temp, par.n_wg_wid, par.n_wg_hgt, 2, par.n_modi).*NaN;
data.fit_neff   = zeros(par.n_wg_wid, par.n_wg_hgt, 2, par.n_modi).*NaN;
data.Aeff       = zeros(par.n_wlen, par.n_temp, par.n_wg_wid, par.n_wg_hgt, 2, par.n_modi).*NaN;
data.fit_Aeff   = zeros(par.n_wg_wid, par.n_wg_hgt, 2, par.n_modi).*NaN;
data.lossdB     = zeros(par.n_wlen, par.n_temp, par.n_wg_wid, par.n_wg_hgt, 2, par.n_modi).*NaN;


sol.wlen        = zeros(par.n_sol,par.n_modi);
sol.wid         = zeros(par.n_sol,par.n_modi);
sol.hgt         = zeros(par.n_sol,par.n_modi);
sol.temp        = zeros(par.n_sol,par.n_modi);
sol.neff        = zeros(par.n_sol,par.n_modi);
sol.Aeff        = zeros(par.n_sol,par.n_modi);
sol.loss_dB     = zeros(par.n_sol,par.n_modi);

% fill in some data
for jj=1:par.n_sol
    sol.wlen(jj,:)  = mphglobal(model,{'la'},'dataset', setting.data_set,'outersolnum',jj).*1e-6;   % 'la' would be in [um]
    sol.wid(jj,:)   = mphglobal(model,{'w'}, 'dataset', setting.data_set,'outersolnum',jj);
    sol.hgt(jj,:)   = mphglobal(model,{'hg'},'dataset', setting.data_set,'outersolnum',jj);  % {'hn'} ?
    sol.temp(jj,:)  = mphglobal(model,{'Temp'},'dataset', setting.data_set,'outersolnum',jj);
    sol.neff(jj,:)  = mphglobal(model,{'emw.neff'},'dataset', setting.data_set,'outersolnum',jj);
    % save effective area
    sol.Aeff(jj,:)  = mphglobal(model,{'Aeff'},'dataset', setting.data_set,'outersolnum',jj);
end
clear tmp jj
toc

%% SAMPLING SOLUTIONS
tic
fprintf('\nFilling the data structures...\n');

% it generates coordinates grid (x & y) of the box
% and coordinates grid (x2 & y2) of the waveguide core
par.ord_max = 0;
par.neff_max   = 1.4;

% take the measures of the box
if (setting.BOX_CONST)
    par.box_w   = mphglobal(model,{'wt'},'dataset', setting.data_set,'outersolnum',1);
    par.box_w   = par.box_w(1);
    par.box_ho  = mphglobal(model,{'ho'},'dataset', setting.data_set,'outersolnum',1);
    par.box_ho  = par.box_ho(1);
    par.box_hts = mphglobal(model,{'hts'},'dataset', setting.data_set,'outersolnum',1);
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
                    data.neff(index.wlen, index.temp, index.wid, index.hgt, 1, sample.npksTM+1)     = real( sol.neff(jj,kk) );
                    data.Aeff(index.wlen, index.temp, index.wid, index.hgt, 1, sample.npksTM+1)     = sol.Aeff(jj,kk);
                    data.lossdB(index.wlen, index.temp, index.wid, index.hgt, 1, sample.npksTM+1)   = sol.loss_dB(jj,kk);
                    
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

clear jj kk ll ans
toc

%% CALCULATING DISPERSION
% it generates dispersion as a function of wavelenght
% dispersion near zero
fprintf('\nGenerating dispersion( wavelength )...\n');

tic
% silicon n2 chi3 real part in um^2/W
disp.n2 = (0.45*10^-17); %*10^12;

disp.lp = 2.1;

% Order of the fit on effective index data
setting.min_fit = 3;
setting.fit_str = 'loess';

[xdata, ydata] = meshgrid(vec.wlen, vec.temp);
xdata = xdata(:);
ydata = ydata(:);

for ww=1:par.n_wg_wid
    for hh=1:par.n_wg_hgt
        for pp=1:2
            for mm=1:par.n_modi
                
                % TE dispersion
                % it controls that the number of points is greater than the
                % order of the fit
                
                % it searches for the points wher neff is greater than 0
                index.good  = find( data.neff(:,:,ww,hh,pp,mm) > 0 );
                index.zero  = find( data.neff(:,:,ww,hh,pp,mm) == 0  );
                index.nan   = find( isnan( data.neff(:,:,ww,hh,pp,mm) ) );

                % eliminates zeros and NaNs
                disp.pte = [];      % points to evaluate
                if ~(isempty(index.good))
                    disp.pte = setdiff( index.good, index.zero);
                    disp.pte = setdiff( disp.pte, index.nan);
                end
                
                % check if there are at least 3 points (to fit a flat surface)
                % and that those points are not on a line
                disp.ok_fit = false;
                if length(index.good)>2
                    if (length( unique( xdata(disp.pte) ) ) > 1 && length( unique( ydata(disp.pte) ) ) > 1 )
                        disp.ok_fit = true;
                    end
                end
                
                % it fit the points to a polynom, if it's possible
                if (disp.ok_fit)                    
                    zdata = data.neff(:,:,ww,hh,pp,mm);
                    zdata = zdata(:);
                    
                    [fit1, gof1] = fit( [xdata(disp.pte),ydata(disp.pte)],zdata(disp.pte),setting.fit_str);
                    data.fit_neff(ww,hh,pp,mm).fit = fit1;
                    data.fit_neff(ww,hh,pp,mm).gof = gof1;
                    
                    zdata = data.Aeff(:,:,ww,hh,pp,mm).*1e12;   % why is there .*1e12 ?
                    zdata = zdata(:);
                    
                    [fit1, gof1] = fit( [xdata(disp.pte),ydata(disp.pte)],zdata(disp.pte),setting.fit_str);
                    data.fit_Aeff(ww,hh,pp,mm).fit = fit1;
                    data.fit_Aeff(ww,hh,pp,mm).gof = gof1;

                    %gamma ??
                    for ll=1:par.n_wlen
                        for tt=1:par.n_temp
                            data.Gamma(ww,hh,pp,mm,ll,tt) = ( 2*pi*disp.n2 )./( vec.wlen(ll).*data.Aeff(ll,tt,ww,hh,pp,mm) );
                        end
                    end
                end
                
            end
        end
    end
end

vec.xdata = xdata;
vec.ydata = ydata;

toc;
clear ww hh pp mm ll ans fit1 gof1 xdata ydata zdata

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
