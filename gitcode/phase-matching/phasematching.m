%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   This file will search for phase matching condition                    %
%                                                                         %
%                                                         D. Bazzanella   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Import HTE & HTM structures

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
        [ans1 ans2] = uigetfile({'../data/*.mat'});     % choose file to load
        
        load( strcat(ans2, ans1) );
    end
end
    

% adjust parameters
vec.wlen = vec.wlen'.*1e-6;                     % [m]
vec.temp = vec.temp';
%vec.temp = vec.temp'-vec.temp(1);

clear ans1 ans2 choice

%% Set initial parameters

tic
fprintf('\nSetting initial parameters\n');
% physical constants
const.c     = 299792458;            % [m/s]

% waves parameters
    % pump structure
    ws.pump.wlen   = 1.55e-6;              % [m]
    ws.pump.freq   = const.c/ws.pump.wlen;    % [1/s]

    % signal structure
    ws.signal.wlen = 1.45e-6;               % [m]
    ws.signal.freq = const.c/ws.signal.wlen;  % [1/s]

    % idler structure
    ws.idler.wlen  = 1.65e-6;               % [m]
    ws.idler.freq  = const.c/ws.idler.wlen;   % [1/s]

par.n_sample = 101; % should be an odd number, to center the pump wlen.
par.stepF = (1.65e-6 -1.45e-6)/(par.n_sample - 1);
vec.sample_wlen = zeros(par.n_sample,1);
for ( ww=1:par.n_sample )
    vec.sample_wlen(ww) = 1.45e-6 + ww*par.stepF;
end

clear ww
toc

%% Temperature dependence analysis

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
                
                %% beta1
                %{
%fprintf(' b1 ');
                if ~( isempty(HTE.DIM(jj,kk).T(1).O(mm).beta1) )
                    x = vec.temp;
                    y = zeros(par.n_temp, 1).*NaN;
                    % insert the parte where it eliminate the zero terms here
                    for tt=1:length(y)
                        if (HTE.DIM(jj,kk).T(tt).O(mm).beta1(ll) ~= 0)
                            y(tt) = HTE.DIM(jj,kk).T(tt).O(mm).beta1(ll);
                        else
                            y(tt) = NaN;
                            x(tt) = NaN;
                        end
                    end

                    y = y(isfinite(y));
                    x = x(isfinite(x));

                    TOC.TE.DIM(jj,kk).O(mm).WL(ll).data1     = length(y);

                    if (length(y) > 1)
                        TOC.TE.DIM(jj,kk).O(mm).WL(ll).fitb1 = fit(x,y,'poly1');
                        TOC.TE.DIM(jj,kk).O(mm).WL(ll).beta1 = TOC.TE.DIM(jj,kk).O(mm).WL(ll).fitb1.p1;
                    else
                        TOC.TE.DIM(jj,kk).O(mm).WL(ll).fitb1 = NaN;
                        TOC.TE.DIM(jj,kk).O(mm).WL(ll).beta1 = NaN;
                    end
                end
                
                %% beta2
%fprintf(' b2 ');
                if ~( isempty(HTE.DIM(jj,kk).T(1).O(mm).beta2) )
                    x = vec.temp;
                    y = zeros(par.n_temp, 1).*NaN;
                    % insert the parte where it eliminate the zero terms here
                    for tt=1:length(y)
                        if (HTE.DIM(jj,kk).T(tt).O(mm).beta2(ll) ~= 0)
                            y(tt) = HTE.DIM(jj,kk).T(tt).O(mm).beta2(ll);
                        else
                            y(tt) = NaN;
                            x(tt) = NaN;
                        end
                    end

                    y = y(isfinite(y));
                    x = x(isfinite(x));

                    TOC.TE.DIM(jj,kk).O(mm).WL(ll).data2     = length(y);

                    if (length(y) > 1)
                        TOC.TE.DIM(jj,kk).O(mm).WL(ll).fitb2 = fit(x,y,'poly1');
                        TOC.TE.DIM(jj,kk).O(mm).WL(ll).beta2 = TOC.TE.DIM(jj,kk).O(mm).WL(ll).fitb2.p1;
                    else
                        TOC.TE.DIM(jj,kk).O(mm).WL(ll).fitb2 = NaN;
                        TOC.TE.DIM(jj,kk).O(mm).WL(ll).beta2 = NaN;
                    end
                end
                                
                %% beta3
%fprintf(' b3 ');
                if ~( isempty(HTE.DIM(jj,kk).T(1).O(mm).beta3) )
                    x = vec.temp;
                    y = zeros(par.n_temp, 1).*NaN;
                    % insert the parte where it eliminate the zero terms here
                    for tt=1:length(y)
                        if (HTE.DIM(jj,kk).T(tt).O(mm).beta3(ll) ~= 0)
                            y(tt) = HTE.DIM(jj,kk).T(tt).O(mm).beta3(ll);
                        else
                            y(tt) = NaN;
                            x(tt) = NaN;
                        end
                    end

                    y = y(isfinite(y));
                    x = x(isfinite(x));

                    TOC.TE.DIM(jj,kk).O(mm).WL(ll).data3     = length(y);

                    if (length(y) > 1)
                        TOC.TE.DIM(jj,kk).O(mm).WL(ll).fitb3 = fit(x,y,'poly1');
                        TOC.TE.DIM(jj,kk).O(mm).WL(ll).beta3 = TOC.TE.DIM(jj,kk).O(mm).WL(ll).fitb3.p1;
                    else
                        TOC.TE.DIM(jj,kk).O(mm).WL(ll).fitb3 = NaN;
                        TOC.TE.DIM(jj,kk).O(mm).WL(ll).beta3 = NaN;
                    end
                end
                
                %% beta4
%fprintf(' b4 ');
                if ~( isempty(HTE.DIM(jj,kk).T(1).O(mm).beta4) )
                    x = vec.temp;
                    y = zeros(par.n_temp, 1).*NaN;
                    % insert the parte where it eliminate the zero terms here
                    for tt=1:length(y)
                        if (HTE.DIM(jj,kk).T(tt).O(mm).beta4(ll) ~= 0)
                            y(tt) = HTE.DIM(jj,kk).T(tt).O(mm).beta4(ll);
                        else
                            y(tt) = NaN;
                            x(tt) = NaN;
                        end
                    end

                    y = y(isfinite(y));
                    x = x(isfinite(x));

                    TOC.TE.DIM(jj,kk).O(mm).WL(ll).data4     = length(y);

                    if (length(y) > 1)
                        TOC.TE.DIM(jj,kk).O(mm).WL(ll).fitb4 = fit(x,y,'poly1');
                        TOC.TE.DIM(jj,kk).O(mm).WL(ll).beta4 = TOC.TE.DIM(jj,kk).O(mm).WL(ll).fitb4.p1;
                    else
                        TOC.TE.DIM(jj,kk).O(mm).WL(ll).fitb4 = NaN;
                        TOC.TE.DIM(jj,kk).O(mm).WL(ll).beta4 = NaN;
                    end
                end
               
%fprintf('\n');
                %}

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

%% Order choosing cycle

% order choosing: cycles between pump, signal and idler's orders to choose
% the physical combinations. ex: P,S,I = 4,2,6
fprintf('\nListing the orders'' configuration permitted\n');

tic
vec.comb = [];
for pp=1:par.n_modi
    for ss=1:par.n_modi
        for ii=1:par.n_modi
            if ~mod(pp,2)
                if ~mod(ss+ii,2)
                    vec.comb = vertcat(vec.comb, [pp ss ii] );
                end
            else
                if mod(ss+ii,2)
                    vec.comb = vertcat(vec.comb, [pp ss ii] );
                end
            end
        end
    end
end

clear pp ss ii ans
toc

%% Creating structure of data for results
fprintf('\nCreating structure of data for results\n');
% define results variable
for jj=1:par.n_wg_wid                                   % number of widths ?
    for kk=1:par.n_wg_hgt                               % number of heights ?
        for oo=1:length( vec.comb(:,1) )
            results.DIM(jj,kk).PO( vec.comb(oo,1) ).orders = [];
            results.DIM(jj,kk).PO( vec.comb(oo,1) ).data = [];
        end
    end
end

clear jj kk oo

%% Main cycle
% insert temperature cycle
temp = 0;   % and eliminate this variable


for jj=1:par.n_wg_wid                                   % number of widths ?
    for kk=1:par.n_wg_hgt                               % number of heights ?
        for oo=1:length( vec.comb(:,1) )
            if ( sum( HTE.DIM(jj,kk).T(1).O(vec.comb(oo,1)).neff() ) >1 && sum( HTE.DIM(jj,kk).T(1).O(vec.comb(oo,2)).neff() ) >1 && sum( HTE.DIM(jj,kk).T(1).O(vec.comb(oo,3)).neff() ) >1)
                tmp = [vec.comb(oo,2) vec.comb(oo,3)];
                results.DIM(jj,kk).PO( vec.comb(oo,1) ).orders = vertcat(results.DIM(jj,kk).PO( vec.comb(oo,1) ).orders, tmp);
                tmp = zeros(par.n_sample,par.n_temp);
                fprintf('m_ord:\tpump %d,\tsignal %d,\tidler %d\n', vec.comb(oo,1), vec.comb(oo,2), vec.comb(oo,3) );
%                fprintf('wlen:\n');
                for ww=1:par.n_sample
                    ws.signal.wlen = vec.sample_wlen(ww);
                    % conservation of energy
                    ws.idler.wlen  = 2*ws.pump.wlen - ws.signal.wlen;
                    
                    % conservation of momentum
                    for tt=1:par.n_temp
%                        fprintf('\tp ');
                        kp = wavenumber(1,jj,kk,vec.comb(oo,1), ws.pump.wlen, vec.temp(tt) , vec, HTE, HTM, TOC );
%                        fprintf('s ');
                        ks = wavenumber(1,jj,kk,vec.comb(oo,2), ws.signal.wlen, vec.temp(tt), vec, HTE, HTM, TOC );
%                        fprintf('i ');
                        ki = wavenumber(1,jj,kk,vec.comb(oo,3), ws.idler.wlen, vec.temp(tt), vec, HTE, HTM, TOC );
                        tmp(ww,tt) = 2*kp -ks -ki ;
%                        fprintf('\n');
                    end
                end
                results.DIM(jj,kk).PO( vec.comb(oo,1) ).data = cat(3, results.DIM(jj,kk).PO( vec.comb(oo,1) ).data, tmp);
%                fprintf('\n');
            end

        end
    end
end
    


clear jj kk oo ww ans tmp temp ki kp ks

%%

xdata = vec.sample_wlen;

for ii=1:par.n_modi
    ydata = results.DIM(1,1).PO(ii).data(:,:)';
    if ~isempty(ydata)
        figure(ii)
        plot(xdata,ydata);
        grid on
    end
end

clear xdata ydata ans ii

%% plot w/ temperature
xdata = vec.sample_wlen;
ydata = vec.temp;

for ii=1:par.n_modi
    zdata = results.DIM(1,1).PO(ii).data(:,:,1);
    if ~isempty(ydata)
        figure(ii)
        surf(xdata,ydata,permute(zdata,[2,1]) );
        grid on
    end
end

clear xdata ydata ans ii
