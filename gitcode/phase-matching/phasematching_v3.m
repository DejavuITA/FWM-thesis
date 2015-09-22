%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   This script will search for phase matching condition                    %
%                                                                         %
%                                                         D. Bazzanella   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Import data structures

ans1 = 0;
while ( ans1 == 0 )
    % Construct a questdlg with three options
    choice = questdlg('Import a file or to exit the program?', 'Import File', 'Import File', 'Exit program', 'Go on', ...
        'Import File');
    % Handle response
    switch choice
        case 'Import File'
            disp('Choose the file to import.\nShould be something like yyyy-MM-dd_HH:mm:ss.mat\n');
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
        [ans1, ans2] = uigetfile({'../data/*.mat'});     % choose file to load
        
        fprintf('File chosen: ', strcat(ans2, ans1));
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
pm.n_sample_wl = 601;              % should be an odd number, to center the pump wlen.
pm.step = (pm.high_w - pm.low_w)/(pm.n_sample_wl - 1);

vec.sample_wlen = zeros(pm.n_sample_wl,1);
for ww=1:pm.n_sample_wl
    vec.sample_wlen(ww) = pm.low_w + (ww-1)*pm.step;
end

% sampling temperature
pm.n_sample_t = 201;
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
% silicon n2 chi3 real part in m^2/W
disp.n2 = (0.45*10^-17)*1e12;   % in um^2/W
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
                    
                    zdata = data.Aeff(:,:,ww,hh,pp,mm).*1e12;   % it comes in m^2 and i want it in um^2
                    zdata = zdata(:);
                    
                    [fit1, gof1] = fit( [xdata(disp.pte),ydata(disp.pte)],zdata(disp.pte),setting.fit_str, 'Normalize', 'on');
                    data.fit_Aeff(ww,hh,pp,mm).exist = true;
                    data.fit_Aeff(ww,hh,pp,mm).fit = fit1;
                    data.fit_Aeff(ww,hh,pp,mm).gof = gof1;

                    %gamma in rad/(m*W)
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
par.ord_max = 7;
pm.comb = [];
for TE=1:2
    for pp=1:par.ord_max       % DA CORREGGERE! corretto?
        for ss=1:par.ord_max
            for ii=1:par.ord_max
                if ~mod(2*pp+ss+ii,2)
                    if (data.fit_neff(1,1,TE,pp).exist && data.fit_neff(1,1,TE,ss).exist && data.fit_neff(1,1,TE,ii).exist)
                            %pm.comb = vertcat(pm.comb, [TE pp TE ss TE ii] );
                            pm.comb = cat(3, pm.comb, [TE, TE, TE; pp ss ii] );
                            % 1 is for TE modes, 2 is for TM modes
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

results.orders  = zeros(2, 3, pm.n_results).*NaN;

results.data    = zeros(pm.n_sample_wl, pm.n_sample_t, pm.n_results, par.n_wg_wid, par.n_wg_hgt).*NaN;
results.maxima  = zeros(pm.n_sample_wl, pm.n_sample_t, pm.n_results, par.n_wg_wid, par.n_wg_hgt).*NaN;
results.bwid    = zeros(pm.n_sample_wl, pm.n_sample_t, pm.n_results, par.n_wg_wid, par.n_wg_hgt).*NaN;

results.max.q 	= zeros(pm.n_results,10).*NaN;
results.max.m 	= zeros(pm.n_results,10).*NaN;
results.bw.q 	= zeros(pm.n_results,10).*NaN;
results.bw.m 	= zeros(pm.n_results,10).*NaN;

results.indexes = zeros(pm.n_results,1).*NaN;

fprintf('\b\b\b:\tstructure created.\n');
clear jj kk oo

%% Main cycle
fprintf('\nRunning main cycle...\n');
tic

% the pump wlen is static 1.55 [um], so i can pre-compute its values out of
% the for loops

% creat grid of evaluation
[~, pgrid]      = meshgrid(vec.sample_temp, pm.pump_w.*ones(pm.n_sample_wl,1));
[ygrid, sgrid]  = meshgrid(vec.sample_temp, vec.sample_wlen);
[~, igrid]      = meshgrid(vec.sample_temp, 1./(2/pm.pump_w - 1./vec.sample_wlen));

for ww=1:par.n_wg_wid                                   % number of widths
    for hh=1:par.n_wg_hgt                               % number of heights
        for oo=1:pm.n_results                           % number of possible combinations
            
            results.orders(:,:,oo) = pm.comb(:,:,oo);

            fprintf('mode orders:\tpump %d, signal %d, idler %d\t\t(sol %d/%d)\n', pm.comb(2,1,oo), pm.comb(2,2,oo), pm.comb(2,3,oo), oo, pm.n_results );
            
            kp = 2*pi*data.fit_neff( ww, hh, pm.comb(1,1,oo), pm.comb(2,1,oo)).fit(pgrid, ygrid)./pgrid;
            ks = 2*pi*data.fit_neff( ww, hh, pm.comb(1,2,oo), pm.comb(2,2,oo)).fit(sgrid, ygrid)./sgrid;
            ki = 2*pi*data.fit_neff( ww, hh, pm.comb(1,3,oo), pm.comb(2,3,oo)).fit(igrid, ygrid)./igrid;
            
            results.data(:,:,oo,ww,hh) = 2.*kp - ks - ki ;
        end
    end
end
    
fprintf('\nMain cycle runned and ended.\n');

toc
clear ww hh oo kp ks ki pgrid sgrid igrid ygrid ans

%% 3D surface plot w/ temperature
% it plots only surfaces that cross the plane z=0
fprintf('\nGenerating plots...\n');
tic

[xGrid, yGrid]  = meshgrid(vec.sample_temp, vec.sample_wlen);

f(pm.n_results) = figure(pm.n_results);
close(f(pm.n_results));

for rr=1:pm.n_results;
    fprintf('\n%d\n',rr);
    if ~isnan(results.data(:,:,rr))    
        zdata = results.data(:,:,rr);
        
        % check if the maximum value of l_coh is greater than 1
        if min( abs( zdata(:) ) ) < 200*pi      % min(zdata(:))*max(zdata(:))<=0
            results.indexes(rr) = 1;
            
            %% 3D surface plots of Dk (wlen, temps)
            
            %{
            name = strcat( 'sol ',num2str(rr),' p ', num2str(results.orders(2,1,rr)), ' s ', num2str(results.orders(2,2,rr)), ' i ', num2str(results.orders(2,3,rr)), ' 3D' );
            f = figure('name', name, 'OuterPosition',[100, 100, 960, 480]);
            
            % fit & plot Delta_k
            
            zData = zdata./100;              % zData in [1/cm]
            
            subplot1 = subplot(1,2,1, 'Parent', f );
            title('\Delta_k [1/cm]');
            hold(subplot1,'on');
            mesh( xGrid, yGrid, zData );
            grid on;
            colorbar('peer',subplot1);
            
            % fit & plot l_coh
            
            zData = 200*pi./abs(zData);            % zdata = [1/m] --> 200pi./abs(zdata) = [cm]
            
            subplot2 = subplot(1,2,2, 'Parent', f );
            title('L_{coh} [cm]');
            hold(subplot2,'on');
            zlim(subplot2,[-0.1 10]);
            mesh( xGrid, yGrid, zData );
            grid on;
            colorbar('peer',subplot2);
            %}
            
            %% 2D plots of peak and bandwidth
            
            % find local maxima & bandwidth, for each temperature
            
            name = strcat( 'sol ',num2str(rr),' p ', num2str(results.orders(2,1,rr)), ' s ', num2str(results.orders(2,2,rr)), ' i ', num2str(results.orders(2,3,rr)), ' 2D' );
            f = figure('name', name, 'OuterPosition',[100, 100, 960, 480]);
            
            subplot1 = subplot(1,2,1, 'Parent', f );
            title('position of maxima');
            
            zData = 200*pi./abs(zdata);            % zData = [1/cm] --> 200pi./abs(zData) = [cm]
            
            maxL = 0;
            y = zeros(pm.n_sample_t, idivide(int32(pm.n_sample_wl),int32(2)) ).*NaN;
            y2 = zeros(pm.n_sample_t, idivide(int32(pm.n_sample_wl),int32(2)) ).*NaN;
            
            for tt=1:pm.n_sample_t
                [~, locs] = findpeaks(zData(:,tt));
                locs = locs( find(zData(locs,tt) >= 1 ) );
                
                maxL = max([maxL, length(locs)]);
                y(tt, 1:length(locs)) = vec.sample_wlen(locs);
                
                zData_0 = zData(:,tt);
                zData_0(1) = 0;
                zData_0(end) = 0;
                
                a = zData_0(:) >= 1;
                b = abs( diff(a) );
                c = find(b == 1);
                %d = vec.sample_wlen( c );
                % we could implement a function which find the zero with
                % bisection, it's a pain in the ass but it's reliable
                d = (vec.sample_wlen( c+1 ) - vec.sample_wlen( c ))./(zData(c+1,tt) - zData(c,tt)).*(1 - zData(c,tt)) + vec.sample_wlen( c );
                e = d(2:2:end) - d(1:2:end);
                
                y2(tt, 1:length(e)) = e;
                clear a b c d e zData_0
            end
            
            y(:, maxL+1:end) = [];
            y2(:, maxL+1:end) = [];
            
            for tt=1:maxL
                a = y(:,tt);
                b = find(isnan(a) == 0);
                if length(b)>1
                    [results.max.ft, ~] = fit( vec.sample_temp(b), a(b), 'poly1');
                    results.max.m(rr, tt)   = results.max.ft.p1;
                    results.max.q(rr, tt)   = y(1,tt);
                else
                    results.max.m(rr, tt)   = NaN;
                    results.max.q(rr, tt)   = NaN;
                end
                
                a = y2(:,tt);
                b = find(isnan(a) == 0);
                if length(b)>1
                    [results.bw.ft, ~] = fit( vec.sample_temp(b), a(b), 'poly1');
                    results.bw.m(rr, tt)	= results.bw.ft.p1;
                    results.bw.q(rr, tt)	= y2(1,tt);
                else
                    results.bw.m(rr, tt)   = NaN;
                    results.bw.q(rr, tt)   = NaN;
                end
                
                clear a b
            end
            
            clear results.max.ft results.bw.ft
            
            hold(subplot1,'on');
            plot(vec.sample_temp, y );
            axis([250 700 1.4e-6 1.7e-6]);
            annotation( f,'textbox',...
                        [0.2 0.82 0.2 0.1],...
                        'String',{strcat('q = ', num2str(results.max.q(rr,1)) ), strcat( 'm = ', num2str(1e11*results.max.m(rr,1)) )});

            
            top = round( max([1,max(y2(:)*2e9)] ),1,'significant');
            
            subplot2 = subplot( 1,2,2, 'Parent', f, ...
                                'YTick',[0 : top/10 : top]);
            title('bandwidth');
            hold(subplot2,'on');
            plot(vec.sample_temp, y2.*1e9);
            ylim(subplot2,[0 top ]);
            annotation( f,'textbox',...
                        [0.65 0.82 0.2 0.1],...
                        'String',{strcat('q = ', num2str(results.bw.q(rr,1)) ), strcat( 'm = ', num2str(1e11*results.bw.m(rr,1)) )});
            
            results.maxima( 1:size(y,1), 1:size(y,2), rr, 1, 1) = y;
            results.bwid( 1:size(y2,1), 1:size(y2,2), rr, 1, 1) = y2;
            
            saveas(f, strcat('sol',num2str(rr),'_',num2str(pm.comb(2,:,rr))), 'png');
            close 1
            
        end 
    end
end
results.indexes = find(results.indexes == 1);
results.n_final = length(results.indexes);
toc

clear ans xGrid yGrid zdata zData rr jj tt y y2 maxL locs bw fitresult f name subplot1 subplot2 top

%% 2D plots data vs orders

tmp = [''];
for jj=1:length(results.indexes)
    tmp = [tmp; num2str(pm.comb( 2,:,results.indexes(jj) ) ) ];
end
tmp = [tmp; '' ];
tmp = cellstr(tmp)';

f = figure('name', 'Peak position at ambient temperature');
axes1 = axes('Parent',f,'XGrid','on','YGrid','on',...
    'XTickLabelRotation', 60,...
    'XTick', [1:1:length(results.indexes)],...
    'XTickLabel', tmp,...
    'LineStyleOrderIndex',2);
hold(axes1,'on');
title('Peak position at ambient temperature [µm]');
plot( 1e6*results.max.q(results.indexes,1:end), 'o' );
xlim([0 length(results.indexes)+1])
ylabel({'Position [µm]'});
xlabel({'Combination'});

%saveas(1, 'ppaat', 'png');

f = figure('name', 'Change in peaks position [nm/100°C]');
axes1 = axes('Parent',f,'XGrid','on','YGrid','on',...
    'XTickLabelRotation', 60,...
    'XTick', [1:1:length(results.indexes)],...
    'XTickLabel', tmp,...
    'LineStyleOrderIndex',2);
hold(axes1,'on');
title('Change in peak position [nm/100°C]');
plot( 1e11*results.max.m(results.indexes,1:end), 'o' );
xlim([0 length(results.indexes)+1])
ylabel({'Δ-position [nm/100°C]'});
xlabel({'Combination'});

%saveas(2, 'ppc', 'png');

f = figure('name', 'Bandwidth [nm]');
axes1 = axes('Parent',f,'XGrid','on','YGrid','on',...
    'XTickLabelRotation', 60,...
    'XTick', [1:1:length(results.indexes)],...
    'XTickLabel', tmp,...
    'LineStyleOrderIndex',2);
hold(axes1,'on');
title('Bandwidth at ambient temperature [nm]');
plot( 1e9*results.bw.q(results.indexes,1:end), 'o' );
xlim([0 length(results.indexes)+1])
ylabel({'Bandwidth [nm]'});
xlabel({'Combination'});

%saveas(3, 'baat', 'png');

f = figure('name', 'Change in bandwidth [nm/100°C]');
axes1 = axes('Parent',f,'XGrid','on','YGrid','on',...
    'XTickLabelRotation', 60,...
    'XTick', [1:1:length(results.indexes)],...
    'XTickLabel', tmp,...
    'LineStyleOrderIndex',2);
hold(axes1,'on');
title('Change in bandwidth [nm/100°C]');
plot( 1e11*results.bw.m(results.indexes,1:end), 'o' );
xlim([0 length(results.indexes)+1])
ylabel({'Δ-bandwidth [nm/100°C]'});
xlabel({'Combination'});

%saveas(4, 'bc', 'png');

% 2D plots data vs orders, w/out degenerate states

results.indexes2 = [];
for jj=1:results.n_final
    if ~(results.max.q(results.indexes(jj),1) >= 1.549e-6 && results.max.q(results.indexes(jj),1) <= 1.551e-6)
        results.indexes2 = [results.indexes2; results.indexes(jj)];
    end
end
clear jj

tmp = [''];
for jj=1:length(results.indexes2)
    tmp = [tmp; num2str(pm.comb( 2,:,results.indexes2(jj) ) ) ];
end
tmp = [tmp; '' ];
tmp = cellstr(tmp)';

f = figure('name', 'peak position at ambient temperature');
axes1 = axes('Parent',f,'XGrid','on','YGrid','on',...
    'XTickLabelRotation', 60,...
    'XTick', [1:1:length(results.indexes2)],...
    'XTickLabel', tmp,...
    'LineStyleOrderIndex',2);
hold(axes1,'on');
title('Peak position at ambient temperature [µm]');
plot( 1e6*results.max.q(results.indexes2,1:end), 'o' );
xlim([0 length(results.indexes2)+1])
ylabel({'Position [µm]'});
xlabel({'Combination'});

%saveas(5, 'ppaat2', 'png');

f = figure('name', 'Change in peak position [nm/100°C]');
axes1 = axes('Parent',f,'XGrid','on','YGrid','on',...
    'XTickLabelRotation', 60,...
    'XTick', [1:1:length(results.indexes2)],...
    'XTickLabel', tmp,...
    'LineStyleOrderIndex',2);
hold(axes1,'on');
title('Change in peak position [nm/100°C]');
plot( 1e11*results.max.m(results.indexes2,1:end), 'o' );
xlim([0 length(results.indexes2)+1])
ylabel({'Δ-position [nm/100°C]'});
xlabel({'Combination'});

%saveas(6, 'ppc2', 'png');

f = figure('name', 'Bandwidth [nm]');
axes1 = axes('Parent',f,'XGrid','on','YGrid','on',...
    'XTickLabelRotation', 60,...
    'XTick', [1:1:length(results.indexes2)],...
    'XTickLabel', tmp,...
    'LineStyleOrderIndex',2);
hold(axes1,'on');
title('Bandwidth at ambient temperature [nm]');
plot( 1e9*results.bw.q(results.indexes2,1:end), 'o' );
xlim([0 length(results.indexes2)+1])
ylabel({'Bandwidth [nm]'});
xlabel({'Combination'});

%saveas(7, 'baat2', 'png');

f = figure('name', 'Change in bandwidth [nm/100°C]');
axes1 = axes('Parent',f,'XGrid','on','YGrid','on',...
    'XTickLabelRotation', 60,...
    'XTick', [1:1:length(results.indexes2)],...
    'XTickLabel', tmp,...
    'LineStyleOrderIndex',2);
hold(axes1,'on');
title('Change in bandwidth [nm/100°C]');
plot( 1e11*results.bw.m(results.indexes2,1:end), 'o' );
xlim([0 length(results.indexes2)+1])
ylabel({'Δ-bandwidth [nm/100°C]'});
xlabel({'Combination'});

%saveas(8, 'bc2', 'png');

% saveas(1, 'ppaat', 'png');
% saveas(2, 'ppc', 'png');
% saveas(3, 'baat', 'png');
% saveas(4, 'bc', 'png');
% saveas(5, 'ppaat2', 'png');
% saveas(6, 'ppc2', 'png');
% saveas(7, 'baat2', 'png');
% saveas(8, 'bc2', 'png');

clear f

%% max & min
%{
tmp = [''];
for jj=1:length(results.indexes2)
    tmp = [tmp; num2str(pm.comb( 2,:,results.indexes2(jj) ) ) ];
end
tmp = [tmp; '' ];
tmp = cellstr(tmp)';
%}

%% SAVING

tic
% file-saving parameters
setting.save_path = '../data/';         % filepath of filesaves
setting.FILE_EXT	= '.mat';         	% Determines the extension of the files: .mat or .dat

choice = questdlg(  'Save to file or end script?', ...
                    ' ', ...
                    'Save to file', 'End script', ...
                    'Save to file' );
% Handle response
switch choice
	case 'Save to File'
        fp = datestr(datetime,'yyyy-mm-dd_HH:MM:ss');
        fp = strcat(setting.save_path,'phasematch_', fp, setting.FILE_EXT);
        
        fprintf('\nSaving data in %s\n',fp);
        save(fp);
        
	case 'End script'
        fprintf('\nScript End.\n');
end    

clear fp ans choice
toc
