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

clear ans1 ans2 choice

%% Set initial parameters

tic
fprintf('\nSetting initial parameters\n');
% physical constants
const.c     = 299792458;        % [m/s]

% phase matching parameters
% wlens
pm.pump_w   = 1.55e-6;          % [m]
pm.signal_w = 1.45e-6;          % [m]
pm.idler_w  = 1.65e-6;          % [m]

% low & high wlen
pm.low_w = 1.45e-6;             % [m]
pm.high_w = 1.65e-6;            % [m]

% number of points in wlen interval [low_w, high_w]
pm.n_sample = 101;              % should be an odd number, to center the pump wlen.
pm.step = (pm.high_w - pm.low_w)/(pm.n_sample - 1);

vec.sample_wlen = zeros(pm.n_sample,1);
for ww=1:pm.n_sample
    vec.sample_wlen(ww) = pm.low_w + ww*pm.step;
end

% sampling temperature
pm.n_sample_t = 11;
pm.step_t = (vec.temp(end) - vec.temp(1) )/(pm.n_sample_t - 1);

vec.sample_temp = zeros(pm.n_sample_t,1);
for tt=1:pm.n_sample_t
    vec.sample_temp(tt) = vec.temp(1) + tt*pm.step_t;
end

clear ww tt
toc

%% Temperature dependence analysis
%{

tic
fprintf('\nAnalysing the temperature dependence of the data\n');
% temperature dependence of HTE modes
for jj=1:par.n_wg_wid               % number of widths
    for kk=1:par.n_wg_hgt           % number of heights
        for mm=1:par.n_modi         % number of modes
            for ll=1:par.n_wlen     % number of wavelengths
%fprintf('%d %d', mm, ll);
                %% alpha
%fprintf(' a ');
                x = vec.temp;
                y = zeros(par.n_temp, 1).*NaN;
                % insert the parte where it eliminate the zero terms here
                for tt=1:length(y)
                    if (HTE.DIM(jj,kk).T(tt).O(mm).neff(ll) ~= 0)
                        y(tt) = HTE.DIM(jj,kk).T(tt).O(mm).neff(ll);
                    else
                        y(tt) = NaN;
                        x(tt) = NaN;
                    end
                end
                
                y = y(isfinite(y));
                x = x(isfinite(x));
                
                TOC.TE.DIM(jj,kk).O(mm).WL(ll).data      = length(y);
                
                if (length(y) > 1)
                    TOC.TE.DIM(jj,kk).O(mm).WL(ll).fit   = fit(x,y,'poly1');
                    TOC.TE.DIM(jj,kk).O(mm).WL(ll).alpha = TOC.TE.DIM(jj,kk).O(mm).WL(ll).fit.p1;
                else
                    TOC.TE.DIM(jj,kk).O(mm).WL(ll).fit   = NaN;
                    TOC.TE.DIM(jj,kk).O(mm).WL(ll).alpha = NaN;
                end
            end
        end
    end
end
toc
tic
% temperature dependence of HTM modes
for jj=1:par.n_wg_wid               % number of widths
    for kk=1:par.n_wg_hgt           % number of heights
        for mm=1:par.n_modi         % number of modes
            for ll=1:par.n_wlen     % number of wavelengths
%fprintf('%d %d', mm, ll);
                
                %% alpha
%fprintf(' a ');
                x = vec.temp;
                y = zeros(length(HTM.DIM(jj,kk).T), 1).*NaN;
                % insert the parte where it eliminate the zero terms here
                for tt=1:length(y)
                    if (HTM.DIM(jj,kk).T(tt).O(mm).neff(ll) ~= 0)
                        y(tt) = HTM.DIM(jj,kk).T(tt).O(mm).neff(ll);
                    else
                        y(tt) = NaN;
                        x(tt) = NaN;
                    end
                end
                
                y = y(isfinite(y));
                x = x(isfinite(x));
                
                TOC.TM.DIM(jj,kk).O(mm).WL(ll).data  = length(y);
                
                if (length(y) > 1)
                    TOC.TM.DIM(jj,kk).O(mm).WL(ll).fit   = fit(x,y,'poly1');
                    TOC.TM.DIM(jj,kk).O(mm).WL(ll).alpha = TOC.TM.DIM(jj,kk).O(mm).WL(ll).fit.p1;
                else
                    TOC.TM.DIM(jj,kk).O(mm).WL(ll).fit   = NaN;
                    TOC.TM.DIM(jj,kk).O(mm).WL(ll).alpha = NaN;
                end
                
                %% beta1
                %{
%fprintf(' b1 ');
                emp =  par.n_temp;
                for tt=1:par.n_temp
                    if isempty( HTM.DIM(jj,kk).T(tt).O(mm).beta1 )
                        emp = emp -1;
                    end
                end

                if ( emp > 1 )
                    x = vec.temp;
                    y = zeros(length(HTM.DIM(jj,kk).T), 1).*NaN;
                    % insert the parte where it eliminate the zero terms here
                    for tt=1:length(y)
                        if (HTM.DIM(jj,kk).T(tt).O(mm).beta1(ll) ~= 0 || isempty(HTM.DIM(jj,kk).T(tt).O(mm).beta1(ll)) )
                            y(tt) = HTM.DIM(jj,kk).T(tt).O(mm).beta1(ll);
                        else
                            y(tt) = NaN;
                            x(tt) = NaN;
                        end
                    end

                    y = y(isfinite(y));
                    x = x(isfinite(x));

                    TOC.TM.DIM(jj,kk).O(mm).WL(ll).data1     = length(y);

                    if (length(y) > 1)
                        TOC.TM.DIM(jj,kk).O(mm).WL(ll).fitb1 = fit(x,y,'poly1');
                        TOC.TM.DIM(jj,kk).O(mm).WL(ll).beta1 = TOC.TM.DIM(jj,kk).O(mm).WL(ll).fitb1.p1;
                    else
                        TOC.TM.DIM(jj,kk).O(mm).WL(ll).fitb1 = NaN;
                        TOC.TM.DIM(jj,kk).O(mm).WL(ll).beta1 = NaN;
                    end
                end
                
                %% beta2
%fprintf(' b2 ');
                if ~( isempty(HTM.DIM(jj,kk).T(1).O(mm).beta2) )
                    x = vec.temp;
                    y = zeros(length(HTM.DIM(jj,kk).T), 1).*NaN;
                    % insert the parte where it eliminate the zero terms here
                    for tt=1:length(y)
                        if (HTM.DIM(jj,kk).T(tt).O(mm).beta2(ll) ~= 0)
                            y(tt) = HTM.DIM(jj,kk).T(tt).O(mm).beta2(ll);
                        else
                            y(tt) = NaN;
                            x(tt) = NaN;
                        end
                    end

                    y = y(isfinite(y));
                    x = x(isfinite(x));

                    TOC.TM.DIM(jj,kk).O(mm).WL(ll).data2     = length(y);

                    if (length(y) > 1)
                        TOC.TM.DIM(jj,kk).O(mm).WL(ll).fitb2 = fit(x,y,'poly1');
                        TOC.TM.DIM(jj,kk).O(mm).WL(ll).beta2 = TOC.TM.DIM(jj,kk).O(mm).WL(ll).fitb2.p1;
                    else
                        TOC.TM.DIM(jj,kk).O(mm).WL(ll).fitb2 = NaN;
                        TOC.TM.DIM(jj,kk).O(mm).WL(ll).beta2 = NaN;
                    end
                end
                
                %% beta3
%fprintf(' b3 ');
                if ~( isempty(HTM.DIM(jj,kk).T(1).O(mm).beta3) )
                    x = vec.temp;
                    y = zeros(length(HTM.DIM(jj,kk).T), 1).*NaN;
                    % insert the parte where it eliminate the zero terms here
                    for tt=1:length(y)
                        if (HTM.DIM(jj,kk).T(tt).O(mm).beta3(ll) ~= 0)
                            y(tt) = HTM.DIM(jj,kk).T(tt).O(mm).beta3(ll);
                        else
                            y(tt) = NaN;
                            x(tt) = NaN;
                        end
                    end

                    y = y(isfinite(y));
                    x = x(isfinite(x));

                    TOC.TM.DIM(jj,kk).O(mm).WL(ll).data3     = length(y);

                    if (length(y) > 1)
                        TOC.TM.DIM(jj,kk).O(mm).WL(ll).fitb3 = fit(x,y,'poly1');
                        TOC.TM.DIM(jj,kk).O(mm).WL(ll).beta3 = TOC.TM.DIM(jj,kk).O(mm).WL(ll).fitb3.p1;
                    else
                        TOC.TM.DIM(jj,kk).O(mm).WL(ll).fitb3 = NaN;
                        TOC.TM.DIM(jj,kk).O(mm).WL(ll).beta3 = NaN;
                    end
                end
                
                %% beta4
%fprintf(' b4 ');
                if ~( isempty(HTM.DIM(jj,kk).T(1).O(mm).beta4) )
                    x = vec.temp;
                    y = zeros(length(HTM.DIM(jj,kk).T), 1).*NaN;
                    % insert the parte where it eliminate the zero terms here
                    for tt=1:length(y)
                        if (HTM.DIM(jj,kk).T(tt).O(mm).beta4(ll) ~= 0)
                            y(tt) = HTM.DIM(jj,kk).T(tt).O(mm).beta4(ll);
                        else
                            y(tt) = NaN;
                            x(tt) = NaN;
                        end
                    end

                    y = y(isfinite(y));
                    x = x(isfinite(x));

                    TOC.TM.DIM(jj,kk).O(mm).WL(ll).data4     = length(y);

                    if (length(y) > 1)
                        TOC.TM.DIM(jj,kk).O(mm).WL(ll).fitb4 = fit(x,y,'poly1');
                        TOC.TM.DIM(jj,kk).O(mm).WL(ll).beta4 = TOC.TM.DIM(jj,kk).O(mm).WL(ll).fitb4.p1;
                    else
                        TOC.TM.DIM(jj,kk).O(mm).WL(ll).fitb4 = NaN;
                        TOC.TM.DIM(jj,kk).O(mm).WL(ll).beta4 = NaN;
                    end
                end
               
%fprintf('\n');
                %}
            end
        end
    end
end

toc
clear ans x y jj kk mm ll tt
%}

%% Order choosing cycle

% order choosing: cycles between pump, signal and idler's orders to choose
% the physical combinations. ex: P,S,I = 4,2,6
fprintf('\nListing the configuration of permitted orders\n');

tic
pm.comb = [];
for pp=1:(par.n_modi-1)       % DA CORREGGERE!
    for ss=1:par.n_modi
        for ii=1:par.n_modi
            if ~mod(pp,2)
                if ~mod(ss+ii,2)
                    pm.comb = vertcat(pm.comb, [1 pp 1 ss 1 ii] );
                    % 1 is for TE modes, 2 is for TM modes
                end
            else
                if mod(ss+ii,2)
                    pm.comb = vertcat(pm.comb, [1 pp 1 ss 1 ii] );
                end
            end
        end
    end
end

pm.n_results = length(pm.comb);

clear pp ss ii ans
toc

%% Creating structure of data for results

fprintf('\nCreating structure of data for results...');
% define results variable

results.orders  = zeros(pm.n_results,6).*NaN;

results.data    = zeros(pm.n_sample,pm.n_sample_t, pm.n_results, par.n_wg_wid, par.n_wg_hgt).*NaN;

fprintf('\b\b\b:\tstructure created.\n');
clear jj kk oo

%% Main cycle
fprintf('\nRunning main cycle...\n');
tic

for ww=1:par.n_wg_wid                                   % number of widths ?
    for hh=1:par.n_wg_hgt                               % number of heights ?
        for oo=1:pm.n_results
            
            % eliminates the combination with fit_neff not available
            tic
            pm.ok_fit = 0;
            if ~( max( pm.comb(oo,:) )>length( data.fit_neff( ww, hh, 1, :) ) ) % cambiare l'1
                for ff=1:3
                    pm.ok_fit = pm.ok_fit + isempty(data.fit_neff( ww, hh, pm.comb(oo,ff*2-1), pm.comb(oo,ff*2) ).fit );
                end
            else
                pm.ok_fit = 1;
            end
            if pm.ok_fit < 1
                results.orders(oo,:) = pm.comb(oo,:);
                
                fprintf('mode orders:\tpump %d,\tsignal %d,\tidler %d (sol %d/%d)\n', pm.comb(oo,2), pm.comb(oo,4), pm.comb(oo,6), oo, pm.n_results );
                for ll=1:pm.n_sample
                    pm.signal_w = vec.sample_wlen(ll);
                    % conservation of energy
                    pm.idler_w  = 2*pm.pump_w - pm.signal_w;
                    
                    % conservation of momentum
                    for tt=1:pm.n_sample_t
%                        fprintf('\tp ');
                        %kp = 2*pi*feval( data.fit_neff( ww, hh, pm.comb(oo,1), pm.comb(oo,2)).fit, pm.pump_w, vec.sample_temp(tt))/pm.pump_w;
                        kp = 2*pi*data.fit_neff( ww, hh, pm.comb(oo,1), pm.comb(oo,2)).fit(pm.pump_w, vec.sample_temp(tt))/pm.pump_w;
                        %ks = 2*pi*feval( data.fit_neff( ww, hh, pm.comb(oo,3), pm.comb(oo,4)).fit, pm.signal_w, vec.sample_temp(tt))/pm.signal_w;
                        ks = 2*pi*data.fit_neff( ww, hh, pm.comb(oo,3), pm.comb(oo,4)).fit(pm.signal_w, vec.sample_temp(tt))/pm.signal_w;
                        %ki = 2*pi*feval( data.fit_neff( ww, hh, pm.comb(oo,5), pm.comb(oo,6)).fit, pm.idler_w, vec.sample_temp(tt))/pm.idler_w;
                        ki = 2*pi*data.fit_neff( ww, hh, pm.comb(oo,5), pm.comb(oo,6)).fit(pm.idler_w, vec.sample_temp(tt))/pm.idler_w;
                        results.data(ll,tt,oo,ww,hh) = 2*kp -ks -ki ;
%                        fprintf('\n');
                    end
                end
%                fprintf('\n');
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

f(pm.n_results) = figure(pm.n_results);
close(f(pm.n_results));

for rr=1:pm.n_results-100;
    fprintf('\n%d\n',rr);
    if ~isnan(results.data(:,:,rr))    
        zdata = results.data(:,:,rr);
        if min(zdata(:))*max(zdata(:))<=0
            f(rr) = figure(rr);
            str_title = strcat('Solution %d: combination (p,s,i) %d %d %d',rr, pm.comb(rr,2), pm.comb(rr,4), pm.comb(rr,6));
            title(str_title);
            set(f(rr), 'Position', [100, 100, 960, 480]);
            [xData, yData, zData] = prepareSurfaceData( xdata, ydata, zdata );
            % zData = [m] --> 2pi*100./zData = [cm]
            % fit & plot Delta_k
            subplot(1,2,1)
            title('delta_k');
            %[fitresult, ~] = fit( [xData, yData], zData, 'loess', 'Normalize', 'on' );
            [fitresult, ~] = fit( [xData, yData], zData./100, 'cubicinterp', 'Normalize', 'on' ); % zData in [1/cm]
            plot( fitresult, [xData, yData], zData./100 );
            grid on;
            
            hold on;
            
            % fit & plot l_coh
            subplot(1,2,2);
            title('L_coh [cm]');
            zData = 2*pi*100./abs(zData);            % 2pi./zData = [m] --> 2pi*100./zData = [cm]
            % eliminates +inf points
            if max(zData) == inf
                indexes = find(max(zData) == zData);
                zData(indexes) = -inf;
                zData(indexes) = max(zData)*10;
            end
            %[fitresult, ~] = fit( [xData, yData], zData, 'loess', 'Normalize', 'on' );
            [fitresult, ~] = fit( [xData, yData], zData, 'cubicinterp', 'Normalize', 'on' );
            plot( fitresult, [xData, yData], zData );
            %zdata = 2*pi*100./abs(zdata);            % 2pi./zData = [m] --> 2pi*100./zData = [cm]
            
            %surf( zdata );
            
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
save(fp, 'const', 'data', 'f', 'par', 'pm', 'results', 'vec');

clear fp ans
toc

 %% various code drafts
 %{
 
[xData, yData, zData] = prepareSurfaceData( xdata, ydata, zdata );
[fitresult, gof] = fit( [xData, yData], zData, 'loess', 'Normalize', 'on' );
plot( fitresult, [xData, yData], zData );
 
fittypes: 'cubicinterp', 'loess', ...
 
 %}