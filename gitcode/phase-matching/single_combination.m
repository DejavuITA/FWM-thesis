%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   This script will search for phase matching conditions                 %
%                                                                         %
%                                                         D. Bazzanella   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Import data structures

ans1 = false;
while ~ans1
    % Construct a questdlg with three options
    choice = questdlg('Import a file or to exit the program?', 'Import File', 'Import File', 'Exit program', 'Go on', ...
        'Import File');
    % Handle response
    switch choice
        case 'Import File'
            disp('Choose the file to import.\nShould be something like yyyy-MM-dd_HH:mm:ss.mat\n');
            ans1 = false;
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
    
    if ~ans1
        [ans1, ans2] = uigetfile({'../data/*.mat'});     % choose file to load
        
        fprintf(strcat('File chosen: ', ans1) );
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

clear ans1 ans2 choice ans

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
pm.n_sample_wl = 2001;              % should be an odd number, to center the pump wlen.
pm.step = (pm.high_w - pm.low_w)/(pm.n_sample_wl - 1);

vec.sample_wlen = zeros(pm.n_sample_wl,1);
for ww=1:pm.n_sample_wl
    vec.sample_wlen(ww) = pm.low_w + (ww-1)*pm.step;
end

% sampling temperature
pm.n_sample_t = 401;
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

% Order of the fit on effective index data
setting.fit_str = 'poly52';

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
                    
                    [fit1, gof1] = fit( [xdata(disp.pte),ydata(disp.pte)],zdata(disp.pte),setting.fit_str);%, 'Normalize', 'on');
                    data.fit_neff(ww,hh,pp,mm).exist = true;                % 'Normalize', 'on' means it center & scales data
                    data.fit_neff(ww,hh,pp,mm).fit = fit1;
                    data.fit_neff(ww,hh,pp,mm).gof = gof1;
                    
                    zdata = data.Aeff(:,:,ww,hh,pp,mm).*1e12;   % it comes in m^2 and i want it in um^2
                    zdata = zdata(:);
                    
                    [fit1, gof1] = fit( [xdata(disp.pte),ydata(disp.pte)],zdata(disp.pte),setting.fit_str);%, 'Normalize', 'on');
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

clear ww hh pp mm ll tt ans fit1 gof1 xdata ydata zdata index disp
toc

%% Order combination choice
% Order combination choice:
% it lets you choose the combination you want to analyze
% or you can define it by writing in the code
% ex: TE P,S,I = [1, 1, 1; 4,2,6] or TM P,S,I = [2, 2, 2; 5, 3, 1]

tic
%overwrite = true;
overwrite = false;

if ~overwrite
    check = false;
    while ~check
        pm.comb = zeros(2,3);

        prompt = 'What type of mode do you want? TE/TM [TE]: ';
        str = input(prompt,'s');

        switch str
            case []
                pm.comb(1,:) = [1,1,1];
                TE = 1;
            case 'TE'
                pm.comb(1,:) = [1,1,1];
                TE = 1;
            case 'TM'
                pm.comb(1,:) = [2,2,2];
                TE = 2;
        end
        
        prompt = 'What is the order of the PUMP (1,2,3...)? ';
        pp = input(prompt);
        pm.comb(2,1) = pp;
        
        prompt = 'What is the order of the SIGNAL (1,2,3...)? ';
        ss = input(prompt);
        pm.comb(2,2) = ss;

        prompt = 'What is the order of the IDLER (1,2,3...)? ';
        ii = input(prompt);
        pm.comb(2,3) = ii;

        if ~mod(2*pp+ss+ii,2)
            fprintf('\nThe configuration is permitted');
            if (data.fit_neff(1,1,TE,pp).exist && data.fit_neff(1,1,TE,ss).exist && data.fit_neff(1,1,TE,ii).exist)
                fprintf('\nThe configuration is permitted');
                check = true;
            end
        end
    end
    
    clear check str order TE
else
    % overwrite here
    TE      = true;
    %TE      = false;
    pump    = 5;
    signal  = 5;
    idler   = 5;
    if TE
        pm.comb = [1, 1, 1; pump, signal, idler];
    else
        pm.comb = [2, 2, 2; pump, signal, idler];
    end
    
    clear TE pump signal idler
end
toc

clear overwrite

%% Main cycle
% phase mismatch analytic evalutaion

tic
fprintf('\nRunning main cycle...\n');

% sets width and height indexes
ww = 1;
hh = 1;

% choses and rename fits
p.fit   = data.fit_neff( ww, hh, pm.comb(1,1), pm.comb(2,1)).fit;
s.fit   = data.fit_neff( ww, hh, pm.comb(1,2), pm.comb(2,2)).fit;
i.fit   = data.fit_neff( ww, hh, pm.comb(1,3), pm.comb(2,3)).fit;

% get coefficients
cn = cell2mat( coeffnames(p.fit) );
cv = coeffvalues(p.fit);
for jj=1:length(cv)
    eval([ 'p.', cn(jj,:), '= cv(jj);']);
end

cn = cell2mat( coeffnames(s.fit) );
cv = coeffvalues(s.fit);
for jj=1:length(cv)
    eval([ 's.', cn(jj,:), '= cv(jj);']);
end

cn = cell2mat( coeffnames(i.fit) );
cv = coeffvalues(i.fit);
for jj=1:length(cv)
    eval([ 'i.', cn(jj,:), '= cv(jj);']);
end

% get formulas
p.fstr      = formula(p.fit);
p.fstr2     = strrep(p.fstr, 'x', 'xp');
p.fstr2     = strrep(p.fstr2, 'y', 't');
p.fstr2     = strrep(p.fstr2, '^', '.^');
p.fstr2     = strrep(p.fstr2, '*', '.*');
%p.fstr2     = strrep(p.fstr2, '+', '.+');
eval([ 'p.fp = @(xp,t,p00,p10,p01,p20,p11,p02,p30,p21,p12,p40,p31,p22,p50,p41,p32) ', p.fstr2],';')
p.fp2       = @(xp,t) p.fp(xp,t,p.p00,p.p10,p.p01,p.p20,p.p11,p.p02,p.p30,p.p21,p.p12,p.p40,p.p31,p.p22,p.p50,p.p41,p.p32);
%p.fp2       = @(xp,t) p.fp((xp-+1.548e-06)*6.01e-08,(t-492.2)/143.1,p.p00,p.p10,p.p01,p.p20,p.p11,p.p02,p.p30,p.p21,p.p12,p.p40,p.p31,p.p22,p.p50,p.p41,p.p32);

s.fstr      = formula(s.fit);
s.fstr2     = strrep(s.fstr, 'x', 'xs');
s.fstr2     = strrep(s.fstr2, 'y', 't');
s.fstr2     = strrep(s.fstr2, '^', '.^');
s.fstr2     = strrep(s.fstr2, '*', '.*');
%s.fstr2     = strrep(s.fstr2, '+', '.+');
eval([ 's.fs = @(xs,t,p00,p10,p01,p20,p11,p02,p30,p21,p12,p40,p31,p22,p50,p41,p32) ', s.fstr2],';')
s.fs2       = @(xs,t) s.fs(xs,t,s.p00,s.p10,s.p01,s.p20,s.p11,s.p02,s.p30,s.p21,s.p12,s.p40,s.p31,s.p22,s.p50,s.p41,s.p32);

i.fstr      = formula(i.fit);
i.fstr2     = strrep(i.fstr, 'x', 'xi');
i.fstr2     = strrep(i.fstr2, 'y', 't');
i.fstr2     = strrep(i.fstr2, '^', '.^');
i.fstr2     = strrep(i.fstr2, '*', '.*');
eval([ 'i.fi = @(xi,t,p00,p10,p01,p20,p11,p02,p30,p21,p12,p40,p31,p22,p50,p41,p32) ', i.fstr2],';')
i.fi2       = @(xi,t) i.fi(xi,t,i.p00,i.p10,i.p01,i.p20,i.p11,i.p02,i.p30,i.p21,i.p12,i.p40,i.p31,i.p22,i.p50,i.p41,i.p32);

xp = 1.55e-6;

f_delta = @(xs,t) 2.*pi.*( 2.*p.fp2(xp,t)./xp - s.fs2(xs,t)./xs - i.fi2(2.*xp-xs,t)./(2.*xp-xs) );
f_Lcoh  = @(xs,t) 2.*pi./abs(f_delta(xs,t));

toc

%% 3D surface plot w/ temperature
% creat grid of evaluation
[wlen_grid, temp_grid]  = meshgrid(vec.sample_wlen, vec.sample_temp);

f = figure(1);
axes1 = axes('Parent',f);
mesh(wlen_grid, temp_grid, f_Lcoh(wlen_grid, temp_grid));
hold(axes1,'on');
zlim(axes1,[-0.1 2]);

clear wlen_grid temp_grid

%% 2D plots of center & bandwidth

