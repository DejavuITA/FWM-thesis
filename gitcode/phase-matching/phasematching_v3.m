%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   This file will search for phase matching condition                    %
%                                                                         %
%                                                         D. Bazzanella   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Import data structures

ans1 = 0;
while ( ans1 == 0 )
    % Construct a questdlg with three options
    choice = questdlg('Import a file or to exit the program?', 'Import File', 	'Import File', 'Exit program', 'Go on', ...
        'Import File');
    % Handle response
    switch choice
        case 'Import File'
            disp('File chosen.')
            ans1 = 0;
        case 'Exit program'
            %disp('\n\nExiting right now.')
            clear ans1 choice
            error('Exiting right now.');
        %{
        case 'Go on'
            disp('\nGoing on...\n');
            ans1 = 1;
        %}
    end
    
    if (ans1 == 0 )
        fprintf('\nChoose the file to import.\nShould be something like yyyy-MM-dd_HH:mm:ss.mat\n');
        [ans1, ans2] = uigetfile({'../data/*.mat'});     % choose file to load
        
        load( strcat(ans2, ans1) );
    end
end

% physical constants
const.c     = 299792458;        % [m/s]
const.T     = 293.15;           % [K]

% rescale unit of measure
vec.wid  = vec.wid .* 1e-6;     % [m]
vec.hgt  = vec.hgt .* 1e-6;     % [m]
vec.temp = vec.temp + const.T;  % [K]

clear ans1 ans2 choice

%% Set initial parameters

tic
fprintf('\nSetting initial parameters\n');

% phase matching parameters
% wlens
pm.pump_w   = 1.55e-6;          % [m]
pm.signal_w = 1.45e-6;          % [m]   it's an example wl
pm.idler_w  = 1.65e-6;          % [m]   it's an example wl

% low & high wlen
pm.low_w = 1.45e-6;             % [m]
pm.high_w = 1.65e-6;            % [m]

% number of points in wlen interval [low_w, high_w]
pm.n_sample_wl = 101;              % should be an odd number, to center the pump wlen.
pm.step = (pm.high_w - pm.low_w)/(pm.n_sample_wl - 1);

vec.sample_wlen = zeros(pm.n_sample_wl,1);
for ww=1:pm.n_sample_wl
    vec.sample_wlen(ww) = pm.low_w + (ww-1)*pm.step;
end

% sampling temperature
pm.n_sample_t = 11;
pm.step_t = (vec.temp(end) - vec.temp(1) )/(pm.n_sample_t - 1);

vec.sample_temp = zeros(pm.n_sample_t,1);
for tt=1:pm.n_sample_t
    vec.sample_temp(tt) = vec.temp(1) + (tt-1)*pm.step_t;
end

clear ww tt
toc

%% Temperature dependence analysis
% it generates dispersion as a function of wavelenght
% dispersion near zero
fprintf('\nCalculating temperature and wavelenght depencence...\n');

tic
% silicon n2 chi3 real part in um^2/W
disp.n2 = (0.45*10^-17); %*10^12;
disp.lp = 2.1;

% initialise data structure
%data.fit_neff.exist       = zeros(par.n_wg_wid, par.n_wg_hgt, 2, par.ord_max).*NaN;
%data.fit_Aeff.exist       = zeros(par.n_wg_wid, par.n_wg_hgt, 2, par.ord_max).*NaN;

% Order of the fit on effective index data
%setting.min_fit = 3;
%setting.fit_str = 'loess';
setting.fit_str = 'poly52';
%setting.fit_str = 'poly25';
%setting.fit_str = 'cubicinterp';

[xdata, ydata, ~] = prepareSurfaceData( vec.wlen, vec.temp, ones(par.n_wlen,par.n_temp) );

for ww=1:par.n_wg_wid
    for hh=1:par.n_wg_hgt
        for pp=1:2
            for mm=1:par.n_modi
                
                data.fit_neff(ww,hh,pp,mm).exist = false;
                data.fit_Aeff(ww,hh,pp,mm).exist = false;
                
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
                if length(index.good)>14    % you need 15 points to fit a poly52
                    if (length( unique( xdata(disp.pte) ) ) > 1 && length( unique( ydata(disp.pte) ) ) > 1 )
                        disp.ok_fit = true;
                    end
                end
                
                % it fit the points to a polynom_xy, if it's possible
                if (disp.ok_fit)                    
                    zdata = data.neff(:,:,ww,hh,pp,mm);
                    zdata = zdata(:);
                    
                    [fit1, gof1] = fit( [xdata(disp.pte),ydata(disp.pte)],zdata(disp.pte),setting.fit_str, 'Normalize', 'on');
                    data.fit_neff(ww,hh,pp,mm).exist = true;                % 'Normalize', 'on' means it center & scales data
                    data.fit_neff(ww,hh,pp,mm).fit = fit1;
                    data.fit_neff(ww,hh,pp,mm).gof = gof1;
                    
                    zdata = data.Aeff(:,:,ww,hh,pp,mm).*1e12;   % why is there .*1e12 ?
                    zdata = zdata(:);
                    
                    [fit1, gof1] = fit( [xdata(disp.pte),ydata(disp.pte)],zdata(disp.pte),setting.fit_str, 'Normalize', 'on');
                    data.fit_Aeff(ww,hh,pp,mm).exist = true;
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
clear ww hh pp mm ll tt ans fit1 gof1 xdata ydata zdata index disp

%% Order choosing cycle

% order choosing: cycles between pump, signal and idler's orders to choose
% the physical combinations. ex: P,S,I = 4,2,6
fprintf('\nListing the configuration of permitted orders\n');

tic
pm.comb = [];
for TE=1:2
    for pp=1:par.ord_max       % DA CORREGGERE! corretto?
        for ss=1:par.ord_max
            for ii=1:par.ord_max
                if ~mod(pp,2)
                    if ~mod(ss+ii,2)
                        if (data.fit_neff(1,1,TE,pp).exist && data.fit_neff(1,1,TE,ss).exist && data.fit_neff(1,1,TE,ii).exist)
                            pm.comb = vertcat(pm.comb, [TE pp TE ss TE ii] );
                            % 1 is for TE modes, 2 is for TM modes
                        end
                    end
                else
                    if mod(ss+ii,2)
                        if (data.fit_neff(1,1,TE,pp).exist && data.fit_neff(1,1,TE,ss).exist && data.fit_neff(1,1,TE,ii).exist)
                            pm.comb = vertcat(pm.comb, [TE pp TE ss TE ii] );
                        end
                    end
                end
            end
        end
    end
end

pm.n_results = length(pm.comb);

clear pp ss ii ans TE
toc

%% Creating structure of data for results

fprintf('\nCreating structure of data for results...');
% define results variable

results.orders  = zeros(pm.n_results,6).*NaN;

results.data    = zeros(pm.n_sample_wl, pm.n_sample_t, pm.n_results, par.n_wg_wid, par.n_wg_hgt).*NaN;

fprintf('\b\b\b:\tstructure created.\n');
clear jj kk oo

%% Main cycle
fprintf('\nRunning main cycle...\n');
tic

% the pump wlen is static 1.55 [um], so i can pre-compute its values out of
% the for loops

for ww=1:par.n_wg_wid                                   % number of widths ?
    for hh=1:par.n_wg_hgt                               % number of heights ?
        for oo=1:pm.n_results
            
            tic
            results.orders(oo,:) = pm.comb(oo,:);

            fprintf('mode orders:\tpump %d,\tsignal %d,\tidler %d (sol %d/%d)\n', pm.comb(oo,2), pm.comb(oo,4), pm.comb(oo,6), oo, pm.n_results );
            for tt=1:pm.n_sample_t
                %kp = 2*pi*feval( data.fit_neff( ww, hh, pm.comb(oo,1), pm.comb(oo,2)).fit, pm.pump_w, vec.sample_temp(tt))/pm.pump_w;
                %kp = 2*pi*data.fit_neff( ww, hh, pm.comb(oo,1), pm.comb(oo,2)).fit(pm.pump_w, vec.sample_temp(tt))/pm.pump_w;
                kp = pm.pump_k(tt,pm.comb(oo,2),pm.comb(oo,1),ww,hh);

                for ll=1:pm.n_sample_wl
                    pm.signal_w = vec.sample_wlen(ll);
                    % conservation of energy
                    pm.idler_w  = 2*pm.pump_w - pm.signal_w;
                
                    % conservation of momentum
                    %ks = 2*pi*feval( data.fit_neff( ww, hh, pm.comb(oo,3), pm.comb(oo,4)).fit, pm.signal_w, vec.sample_temp(tt))/pm.signal_w;
                    ks = 2*pi*data.fit_neff( ww, hh, pm.comb(oo,3), pm.comb(oo,4)).fit(pm.signal_w, vec.sample_temp(tt))/pm.signal_w;
                    %ki = 2*pi*feval( data.fit_neff( ww, hh, pm.comb(oo,5), pm.comb(oo,6)).fit, pm.idler_w, vec.sample_temp(tt))/pm.idler_w;
                    ki = 2*pi*data.fit_neff( ww, hh, pm.comb(oo,5), pm.comb(oo,6)).fit(pm.idler_w, vec.sample_temp(tt))/pm.idler_w;
                    results.data(ll,tt,oo,ww,hh) = 2*kp -ks -ki ;
                end
            end
            toc
        end
    end
end
    
fprintf('\nMain cycle runned and ended.\n');

toc
clear ww hh oo ll ff tt kp ks ki ans

%% surface plot w/ temperature
% it plots only surfaces that cross the plane z=0
xdata = vec.sample_wlen;
ydata = vec.sample_temp;
[xData, yData, ~] = prepareSurfaceData( xdata, ydata, results.data(:,:,1) );

f(pm.n_results) = figure(pm.n_results);
close(f(pm.n_results));

%setting.fit_str = 'linearinterp'; % Linear interpolation          
%setting.fit_str = 'nearestinterp'; % Nearest neighbor interpolation         
setting.fit_str = 'cubicinterp'; % Cubic spline interpolation
%setting.fit_str = 'biharmonicinterp'; % Biharmonic (MATLAB griddata) interpolation % ALSO NOT GOOD
%setting.fit_str = 'thinplateinterp'; % Thin-plate spline interpolation % very computation needy DO NOT USE

%setting.fit_str = 'poly55';    % doesn't work DO NOT USE

for rr=1:pm.n_results;
    fprintf('\n%d\n',rr);
    if ~isnan(results.data(:,:,rr))    
        zdata = results.data(:,:,rr);
        
        % check if the maximum value of l_coh is greater than 1
        if min( abs( zdata(:) ) ) < 200*pi      % min(zdata(:))*max(zdata(:))<=0
            
            zData = zdata(:)./100;              % zData in [1/cm]
            [fitresult, ~] = fit( [xData, yData], zData, setting.fit_str, 'Normalize', 'on' );
            
            f(rr) = figure(rr);
            set(f(rr), 'Position', [100, 100, 960, 480]);
            
            % fit & plot Delta_k
            subplot(1,2,1);
            title('delta_k');
            plot( fitresult, [xData, yData], zData );
            grid on;            
            hold on;
            
            zData = 2*pi./abs(zData);            % zData = [1/cm] --> 2pi./abs(zData) = [cm]
            % eliminates +inf points
            if max(zData) == inf
                indexes = find(max(zData) == zData);
                zData(indexes) = -inf;
                zData(indexes) = max(zData)*10;
            end
            [fitresult, ~] = fit( [xData, yData], zData, setting.fit_str, 'Normalize', 'on' );
            
            % fit & plot l_coh
            subplot(1,2,2);
            title('L_coh [cm]');
            plot( fitresult, [xData, yData], zData );            
            grid on;
        end
    end
end

clear xdata ydata zdata xData yData zData ans rr fitresult indexes str_title

%% SAVING

tic
% file-saving parameters
setting.save_path = '../data/';         % filepath of filesaves
setting.FILE_EXT	= '.mat';         	% Determines the extension of the files: .mat or .dat

fp = datestr(datetime,'yyyy-mm-dd_HH:MM:ss');
fp = strcat(setting.save_path,'phasematch_', fp, setting.FILE_EXT);

fprintf('Saving data in %s\n',fp);
save(fp);

clear fp ans
toc

 %% various code drafts
 %{

[xData, yData, zData] = prepareSurfaceData( xdata, ydata, zdata );
[fitresult, gof] = fit( [xData, yData], zData, 'loess', 'Normalize', 'on' );
plot( fitresult, [xData, yData], zData );

fittypes: 'cubicinterp', 'loess', ...

this things don't work { 
    str_title = strcat('Solution %d: combination (p,s,i) %d %d %d',rr, pm.comb(rr,2), pm.comb(rr,4), pm.comb(rr,6));
    title(str_title);
}

 %}