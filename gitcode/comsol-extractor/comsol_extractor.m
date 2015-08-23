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

[ans1 ans2] = uigetfile({'../comsol/*.mph'});  % choose file to load
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
fprintf('\tModel Loaded.\n');
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

% allocate vectors of zero
vec.wlen	= zeros(par.n_wlen,1);				% Vector for the lambdas
vec.wid		= zeros(par.n_wg_wid,1);			% Vector for the widths
vec.hgt		= zeros(par.n_wg_hgt,1);			% Vector for the heights
vec.temp    = zeros(par.n_temp,1);				% Vector for the temperatures

% get parameters
vec.wlen    = get_par(model, 'la', setting.data_set);
vec.wid     = get_par(model, 'w', setting.data_set);
vec.hgt     = get_par(model, 'hg', setting.data_set);                      % {'hn'} ?
vec.temp    = get_par(model, 'Temp', setting.data_set);

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
fprintf('\tParameters Loaded.\n');
toc

%% INITIALIZE DATA STRUCTURE
tic
fprintf('\nInitialize structure for data...\n');
void	= zeros(par.n_wlen,1);  % par.n_wlen is the number of scanned wavelengths

% Generate the structure full of zeros, but with the right dimensions (par.n_wlen)
for ii=1:par.n_wg_wid
    for jj=1:par.n_wg_hgt
        for kk=1:par.n_temp
            for ll=1:par.n_modi
                HTE.DIM(ii,jj).T(kk).O(ll).neff    = void;
                HTM.DIM(ii,jj).T(kk).O(ll).neff    = void;
                HTM.DIM(ii,jj).T(kk).O(ll).Aeff    = void;
                HTE.DIM(ii,jj).T(kk).O(ll).Aeff    = void;
                HTM.DIM(ii,jj).T(kk).O(ll).lossdB  = void;
                HTE.DIM(ii,jj).T(kk).O(ll).lossdB  = void;
                HTE.DIM(ii,jj).T(kk).O(ll).profile_saved = false;
                HTM.DIM(ii,jj).T(kk).O(ll).profile_saved = false;
            end
        end
    end
end
clear ii jj kk ll void

sol.wlen        = zeros(par.n_sol,par.n_modi);
sol.wid         = zeros(par.n_sol,par.n_modi);
sol.hgt         = zeros(par.n_sol,par.n_modi);
sol.temp        = zeros(par.n_sol,par.n_modi);
sol.neff        = zeros(par.n_sol,par.n_modi);
sol.Aeff        = zeros(par.n_sol,par.n_modi);
sol.loss_dB     = zeros(par.n_sol,par.n_modi);

% fill in some data
for jj=1:par.n_sol
    sol.wlen(jj,:)  = mphglobal(model,{'la'},'dataset', setting.data_set,'outersolnum',jj);
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
    sol.loss_dB(jj,:) = -20*log10( exp(1) ).*abs( imag( sol.neff(jj,:) ) )*2*pi./( sol.wlen(jj,:)*10^-4 );
    
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
        if( sol.loss_dB(jj,kk)>const.dB )

            % rapp_controllo_livello
            sample.ratio(kk) = max(abs(sample.Exfield(kk,:)))/max(abs(sample.Eyfield(kk,:))); %Rapporto tra i massimi globali dello stesso modo
            index.progress_solution = index.progress_solution + 1;
            if ( sample.ratio(kk)>1 )                
                % classify TE modes

                sample.npksTE = size( fnzeros( spline( sample.ddd(1,:), sample.Exfield(kk,:) ) ), 2 );
fprintf('\tsampling a TE mode; npksTE = %d\n',sample.npksTE);

                if( (sample.npksTE + 1)<index.max_n_modi && (sample.npksTE+1)>0 )

                    % SE abbiamo dei picchi, e sono meno dei modi... ?
                    HTE.DIM(index.wid,index.hgt).T(index.temp).O(sample.npksTE+1).neff(index.wlen)      = real( sol.neff(jj,kk) );
                    HTE.DIM(index.wid,index.hgt).T(index.temp).O(sample.npksTE+1).Aeff(index.wlen)      = sol.Aeff(jj,kk);
                    HTE.DIM(index.wid,index.hgt).T(index.temp).O(sample.npksTE+1).lossdB(index.wlen)    = sol.loss_dB(jj,kk);

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
                    
                    % Save field profile once
                    %{
                    if (HTE.DIM(index.wid,index.hgt).T(index.temp).O(sample.npksTE+1).profile_saved == false && setting.SAVE_FIELDS)
fprintf('\tsaving field');            
                        %[sample.E_real,sample.d2TE] = mphinterp(model,{'real(emw.Ex)'},'coord',sample.points,'dataset', setting.data_set,'outersolnum',jj,'solnum',kk);
                        %[sample.E_imag,sample.d2TE] = mphinterp(model,{'imag(emw.Ex)'},'coord',sample.points,'dataset', setting.data_set,'outersolnum',jj,'solnum',kk);
                        [sample.E_real,sample.d2TE] = mphinterp(model,{'real(emw.Ex)'},'coord',sample.ddd,'dataset', setting.data_set,'outersolnum',jj,'solnum',kk);
                        [sample.E_imag,sample.d2TE] = mphinterp(model,{'imag(emw.Ex)'},'coord',sample.ddd,'dataset', setting.data_set,'outersolnum',jj,'solnum',kk);
                        sample.E_real = vectom(sample.E_real,setting.smpy,setting.smpx); %Turn to matrix
                        sample.E_imag = vectom(sample.E_imag,setting.smpy,setting.smpx);
                        HTE.DIM(index.wid,index.hgt).T(index.temp).O(sample.npksTE+1).profile = (sample.E_real'+1i.*sample.E_imag').*10^-6; %V/um
                        HTE.DIM(index.wid,index.hgt).T(index.temp).O(sample.npksTE+1).profile_saved = true;
                        HTE.DIM(index.wid,index.hgt).T(index.temp).O(sample.npksTE+1).integ_step = [abs(sample.x_tot(2)-sample.x_tot(1)),abs(sample.y_tot(2)-sample.y_tot(1))];  
fprintf('\tfield saved\n');
                    end
                    %}
                end
            else
                % classify TM modes
                
                sample.npksTM = size( fnzeros( spline( sample.ddd(1,:), sample.Eyfield(kk,:) ) ), 2 );
fprintf('\tsampling a TM mode; npksTM = %d\n',sample.npksTM);
                
                if( sample.npksTM+1<index.max_n_modi && sample.npksTM+1>0 && setting.TM_TOO)
                    HTM.DIM(index.wid,index.hgt).T(index.temp).O(sample.npksTM+1).neff(index.wlen)      = real( sol.neff(jj,kk) );
                    HTM.DIM(index.wid,index.hgt).T(index.temp).O(sample.npksTM+1).Aeff(index.wlen)      = sol.Aeff(jj,kk);
                    HTM.DIM(index.wid,index.hgt).T(index.temp).O(sample.npksTM+1).lossdB(index.wlen)    = sol.loss_dB(jj,kk);
                    
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

                    % Save field profile once
                    %{
                    if HTM.DIM(index.wid,index.hgt).T(index.temp).O(sample.npksTE+1).profile_saved == false
fprintf('\tsaving field');
                        [sample.E_real,sample.d2TM] = mphinterp(model,{'real(emw.Ey)'},'coord',sample.points,'dataset', setting.data_set,'outersolnum',jj,'solnum',kk);
                        [sample.E_imag,sample.d2TM] = mphinterp(model,{'imag(emw.Ey)'},'coord',sample.points,'dataset', setting.data_set,'outersolnum',jj,'solnum',kk);
                        sample.E_real = vectom(sample.E_real,setting.smpy,setting.smpx); %Turn to matrix
                        sample.E_imag = vectom(sample.E_imag,setting.smpy,setting.smpx);
                        HTM.DIM(index.wid,index.hgt).T(index.temp).O(sample.npksTE+1).profile = (sample.E_real'+1i.*sample.E_imag').*10^-6;
                        HTM.DIM(index.wid,index.hgt).T(index.temp).O(sample.npksTE+1).profile_saved = true;
                        HTM.DIM(index.wid,index.hgt).T(index.temp).O(sample.npksTE+1).integ_step = [abs(sample.x_tot(2)-sample.x_tot(1)),abs(sample.y_tot(2)-sample.y_tot(1))];
fprintf('\tfield saved\n');
                    end
                    %}
                end
            end
        end
    end % end of for cycle on modes (par.n_modi)
end % end of for cycle on solutions (par.n_sol)
index.progress_solution = 0;

% roba per plotpara per sapere quanti punti deve fare
%HTE(1).W(1).O(1).w=w;
%HTE(1).W(1).O(1).h=h;
%HTM(1).W(1).O(1).h=h;
%HTM(1).W(1).O(1).w=w;

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
% speed of light in um/s
const.uc = const.c*1e6;
% silicon n2 chi3 real part in um^2/W
disp.n2 = (0.45*10^-17); %*10^12;

disp.lp = 2.1;

% Order of the fit on effective index data
setting.fit_order = 5;
setting.fit_str=['poly',num2str(setting.fit_order)];

for ii=1:par.n_wg_hgt
    for jj=1:par.n_wg_wid
        for kk=1:par.n_temp
            for nn=1:par.n_modi
                
                % TE dispersion
                % it controls that the number of points is greater than the
                % order of the fit
                
                % it searches for the points wher neff is greater than 0
                index.nonzero   = find( HTE.DIM(ii,jj).T(kk).O(nn).neff > 0  );
                index.zero      = find( HTE.DIM(ii,jj).T(kk).O(nn).neff == 0 );

                % Se canna piu punti (anche consecutivi) ?           
                disp.pte = [];      % points to evaluate
                if not(isempty(index.nonzero))
                    disp.pte = setdiff([index.nonzero(1):length(vec.wlen)],index.zero);
                end
                
                % aggiungere interpolazione dove manca il punto
                %{
                togli = [];
                index.found = 0;
                if not(isempty(index.nonzero))
                    if index.nonzero(1)<length(disp.pte)
                    for ll=1:length(disp.pte)-1
                        if abs(HTE.DIM(ii,jj).T(kk).O(nn).neff(disp.pte(ll+1))-HTE.DIM(ii,jj).T(kk).O(nn).neff(disp.pte(ll)))>0.03
                            if index.found==0
                                index.found = 1;
                                togli = [togli,disp.pte(ll+1)];
                            else
                                index.found = 0;
                            end
                        else
                            if index.found==1
                                togli = [togli,disp.pte(ll+1)];
                            end
                        end
                    end
                    end
                end
                if not(isempty(togli))
                    disp.pte = setdiff(disp.pte,togli);
                else
                    if not(isempty(index.nonzero))
                    disp.pte = disp.pte;
                    else
                        disp.pte = [];
                    end
                end
                %}
                
                % if the number of points to evaluate greater than the fit
                % order then it proceeds, otherwise it the fit order will
                % be reduced to a minumim of 3 or the fit will be skipped
                if ( length(disp.pte)<setting.fit_order && length(disp.pte)>2 )
                    setting.fit_str=['poly',num2str( length( disp.pte ) -1 )];
                end
                
                % it fit the points to a polynom, if it's possible
                if(length(disp.pte)>2)
                    
                    fit1 = fit(vec.wlen(disp.pte)',HTE.DIM(ii,jj).T(kk).O(nn).neff(disp.pte),setting.fit_str);
                    HTE.DIM(ii,jj).T(kk).O(nn).fit_neff = fit1;

                    % sixth order polynom
                    if(strcmp(setting.fit_str,'poly6'))
                        p1 = fit1.p1;
                        p2 = fit1.p2;
                        p3 = fit1.p3;
                        p4 = fit1.p4;
                        p5 = fit1.p5;
                        p6 = fit1.p6;
                        p7 = fit1.p7;
                        beta1 = ((-vec.wlen.^.2).*(5.*p1.*vec.wlen.^4 + 4.*(vec.wlen.^3).*p2 + 3.*(vec.wlen.^2).*p3 + 2.*vec.wlen.*p4 + p5) + p7)./const.uc;
                        beta2 = ((vec.wlen.^3).*(vec.wlen.*(vec.wlen.*(5.*vec.wlen.*(3.*vec.wlen.*p1 + 2*p2) + 6*p3) + 3*p4) + p5))./(pi*const.uc^2);
                        beta3 = -(((3.*vec.wlen.^4).*(vec.wlen.*(5.*vec.wlen.*(7.*(vec.wlen.^2).*p1 + 4.*vec.wlen.*p2 + 2*p3) + 4*p4) + p5))./(2*(pi^2)*const.uc^3));
                        beta4 = (3.*(vec.wlen.^5).*(5.*vec.wlen.*(vec.wlen.*(7.*vec.wlen.*(2.*vec.wlen.*p1 + p2) + 3*p3) + p4) + p5))./((pi^3)*const.uc^4);
                    end
                    % fifth order polynom
                    if(strcmp(setting.fit_str,'poly5'))
                        p1 = fit1.p1;
                        p2 = fit1.p2;
                        p3 = fit1.p3;
                        p4 = fit1.p4;
                        p5 = fit1.p5;
                        p6 = fit1.p6;
                        beta1 = ((-vec.wlen.^2).*(4.*(vec.wlen.^3).*p1 + 3.*(vec.wlen.^2).*p2 + 2.*vec.wlen.*p3 + p4) + p6)./const.uc;
                        beta2 = ((vec.wlen.^3).*(10.*(vec.wlen.^3).*p1 + 6.*(vec.wlen.^2).*p2 + 3.*vec.wlen.*p3 + p4))./(pi*const.uc^2);
                        beta3 = -((3.*(vec.wlen.^4).*(2.*vec.wlen.*(5.*vec.wlen.*(2.*vec.wlen.*p1 + p2) + 2*p3) + p4))./(2*(pi^2)*const.uc^3));
                        beta4 = (3.*(vec.wlen.^5).*(5.*vec.wlen.*(7.*(vec.wlen.^2).*p1 + 3.*vec.wlen.*p2 + p3) + p4))./((pi^3)*const.uc^4);
                    end
                    % fourth order polynom
                    if(strcmp(setting.fit_str,'poly4'))
                        p1 = fit1.p1;
                        p2 = fit1.p2;
                        p3 = fit1.p3;
                        p4 = fit1.p4;
                        p5 = fit1.p5;
                        beta1 = ((-vec.wlen.^2).*(3.*(vec.wlen.^2).*p1 + 2.*vec.wlen.*p2 + p3) + p5)./const.uc;
                        beta2 = ((vec.wlen.^3).*(3.*vec.wlen.*(2.*vec.wlen.*p1 + p2) + p3))./(pi*const.uc^2);
                        beta3 = -((3.*(vec.wlen.^4).*(10.*(vec.wlen.^2).*p1 + 4.*vec.wlen.*p2 + p3))./(2*(pi^2)*const.uc^3));
                        beta4 = (3.*(vec.wlen.^5).*(5.*vec.wlen.*(3.*vec.wlen.*p1 + p2) + p3))./((pi^3)*const.uc^4);
                    end
                    % third order polynom
                    if(strcmp(setting.fit_str,'poly3'))
                        p1 = fit1.p1;
                        p2 = fit1.p2;
                        p3 = fit1.p3;
                        p4 = fit1.p4;
                        beta1 = (-(vec.wlen.^2).*(2.*vec.wlen.*p1 + p2) + p4)./const.uc;
                        beta2 = ((vec.wlen.^3).*(3.*vec.wlen.*p1 + p2))./(pi*const.uc^2);
                        beta3 = -((3.*(vec.wlen.^4).*(4.*vec.wlen.*p1 + p2))./(2*(pi^2)*const.uc^3));
                        beta4 = (3.*(vec.wlen.^5).*(5.*vec.wlen.*p1 + p2))./((pi^3)*const.uc^4);
                    end
                    % second order polynom
                    if(strcmp(setting.fit_str,'poly2'))
                        p1 = fit1.p1;
                        p2 = fit1.p2;
                        p3 = fit1.p3;
                        beta1 = (-(vec.wlen.^2).*p1 + p3)./const.uc;
                        beta2 = ((vec.wlen.^3).*p1)./(pi*const.uc^2);
                        beta3 = -((3.*(vec.wlen.^4).*p1)./(2*(pi^2)*const.uc^3));
                        beta4 = (3.*(vec.wlen.^5).*p1)./((pi^3)*const.uc^4);
                    end

                    HTE.DIM(ii,jj).T(kk).O(nn).beta1 = beta1.*10^16;    % ps/cm
                    HTE.DIM(ii,jj).T(kk).O(nn).beta2 = beta2.*10^28;    % ps^2/cm
                    HTE.DIM(ii,jj).T(kk).O(nn).beta3 = beta3.*10^40;    % ps^3/cm
                    HTE.DIM(ii,jj).T(kk).O(nn).beta4 = beta4.*10^52;    % ps^4/cm

                    HTE.DIM(ii,jj).T(kk).O(nn).fit_Aeff     = fit(vec.wlen(disp.pte)',HTE.DIM(ii,jj).T(kk).O(nn).Aeff(disp.pte).*10^12,setting.fit_str);

                    HTE.DIM(ii,jj).T(kk).O(nn).fit_beta1    = fit(vec.wlen(disp.pte)',HTE.DIM(ii,jj).T(kk).O(nn).beta1(disp.pte)',setting.fit_str);
                    HTE.DIM(ii,jj).T(kk).O(nn).fit_beta2    = fit(vec.wlen(disp.pte)',HTE.DIM(ii,jj).T(kk).O(nn).beta2(disp.pte)',setting.fit_str);
                    HTE.DIM(ii,jj).T(kk).O(nn).fit_beta3    = fit(vec.wlen(disp.pte)',HTE.DIM(ii,jj).T(kk).O(nn).beta3(disp.pte)',setting.fit_str);
                    HTE.DIM(ii,jj).T(kk).O(nn).fit_beta4    = fit(vec.wlen(disp.pte)',HTE.DIM(ii,jj).T(kk).O(nn).beta4(disp.pte)',setting.fit_str);

                                                                        % is correct the 2 down here or was correct the .2 before?
                    HTE.DIM(ii,jj).T(kk).O(nn).Disp     = -beta2.*10^28.*2.*pi.*3e8./(vec.wlen.^2).*1e-4;
                    HTE.DIM(ii,jj).T(kk).O(nn).Gamma    = (2*pi*disp.n2)./(vec.wlen'.*1e-6.*HTE.DIM(ii,jj).T(kk).O(nn).Aeff);

                    HTE.DIM(ii,jj).T(kk).O(nn).Popti    = -((2*pi*3e8)*(beta2'.*10^28).*(1./(vec.wlen')-1/disp.lp).^2+(1/12)*(2*pi*3e8).*(beta4'.*10^52).*(1./(vec.wlen')-1/disp.lp).^4)./HTE.DIM(ii,jj).T(kk).O(nn).Gamma;
                    % it finds zero dispersion
                    disp.zr = fnzeros(spline(vec.wlen(disp.pte),HTE.DIM(ii,jj).T(kk).O(nn).beta2(disp.pte)));

                    if(isempty(disp.zr))
                        HTE.DIM(ii,jj).T(kk).O(nn).zero = 0;
                    else
                        HTE.DIM(ii,jj).T(kk).O(nn).zero = disp.zr(1,:); % all the points of zero dispersion
                    end
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%
                % TM dispersion
                % it controls that the number of points is greater than the
                % order of the fit
                
                % it searches for the points wher neff is greater than 0
                index.nonzero   = find( HTM.DIM(ii,jj).T(kk).O(nn).neff > 0  );
                index.zero      = find( HTM.DIM(ii,jj).T(kk).O(nn).neff == 0 );

                % Se canna piu punti (anche consecutivi) ?           
                disp.pte = [];      % points to evaluate
                if not(isempty(index.nonzero))
                    disp.pte = setdiff([index.nonzero(1):length(vec.wlen)],index.zero);
                end
                
                %{
                togli = [];
                index.found = 0;
                if not(isempty(index.nonzero))
                    if index.nonzero(1)<length(disp.pte)
                    for ll=1:length(disp.pte)-1
                        if abs(HTM.DIM(ii,jj).T(kk).O(nn).neff(disp.pte(ll+1))-HTM.DIM(ii,jj).T(kk).O(nn).neff(disp.pte(ll)))>0.03
                            if index.found==0
                                index.found = 1;
                                togli = [togli,disp.pte(ll+1)];
                            else
                                index.found = 0;
                            end
                        else
                            if index.found==1
                                togli = [togli,disp.pte(ll+1)];
                            end
                        end
                    end
                    end
                end
                if not(isempty(togli))
                    disp.pte = setdiff(disp.pte,togli);
                else
                    if not(isempty(index.nonzero))
                    disp.pte = disp.pte;
                    else
                        disp.pte = [];
                    end
                end
                %}
                
                % if the number of points to evaluate greater than the fit
                % order then it proceeds, otherwise it the fit order will
                % be reduced to a minumim of 3 or the fit will be skipped
                if ( length(disp.pte)<setting.fit_order && length(disp.pte)>2 )
                    setting.fit_str=['poly',num2str( length( disp.pte ) -1 )];
                end
                
                % it fit the points to a polynom, if it's possible
                if(length(disp.pte)>2 && setting.TM_TOO)
                    
                    fit1 = fit(vec.wlen(disp.pte)',HTM.DIM(ii,jj).T(kk).O(nn).neff(disp.pte),setting.fit_str);
                    HTM.DIM(ii,jj).T(kk).O(nn).fit_neff = fit1;

                    % sixth order polynom
                    if(strcmp(setting.fit_str,'poly6'))
                        p1 = fit1.p1;
                        p2 = fit1.p2;
                        p3 = fit1.p3;
                        p4 = fit1.p4;
                        p5 = fit1.p5;
                        p6 = fit1.p6;
                        p7 = fit1.p7;
                        beta1 = ((-vec.wlen.^.2).*(5.*p1.*vec.wlen.^4 + 4.*(vec.wlen.^3).*p2 + 3.*(vec.wlen.^2).*p3 + 2.*vec.wlen.*p4 + p5) + p7)./const.uc;
                        beta2 = ((vec.wlen.^3).*(vec.wlen.*(vec.wlen.*(5.*vec.wlen.*(3.*vec.wlen.*p1 + 2*p2) + 6*p3) + 3*p4) + p5))./(pi*const.uc^2);
                        beta3 = -(((3.*vec.wlen.^4).*(vec.wlen.*(5.*vec.wlen.*(7.*(vec.wlen.^2).*p1 + 4.*vec.wlen.*p2 + 2*p3) + 4*p4) + p5))./(2*(pi^2)*const.uc^3));
                        beta4 = (3.*(vec.wlen.^5).*(5.*vec.wlen.*(vec.wlen.*(7.*vec.wlen.*(2.*vec.wlen.*p1 + p2) + 3*p3) + p4) + p5))./((pi^3)*const.uc^4);
                    end
                    % fifth order polynom
                    if(strcmp(setting.fit_str,'poly5'))
                        p1 = fit1.p1;
                        p2 = fit1.p2;
                        p3 = fit1.p3;
                        p4 = fit1.p4;
                        p5 = fit1.p5;
                        p6 = fit1.p6;
                        beta1 = ((-vec.wlen.^2).*(4.*(vec.wlen.^3).*p1 + 3.*(vec.wlen.^2).*p2 + 2.*vec.wlen.*p3 + p4) + p6)./const.uc;
                        beta2 = ((vec.wlen.^3).*(10.*(vec.wlen.^3).*p1 + 6.*(vec.wlen.^2).*p2 + 3.*vec.wlen.*p3 + p4))./(pi*const.uc^2);
                        beta3 = -((3.*(vec.wlen.^4).*(2.*vec.wlen.*(5.*vec.wlen.*(2.*vec.wlen.*p1 + p2) + 2*p3) + p4))./(2*(pi^2)*const.uc^3));
                        beta4 = (3.*(vec.wlen.^5).*(5.*vec.wlen.*(7.*(vec.wlen.^2).*p1 + 3.*vec.wlen.*p2 + p3) + p4))./((pi^3)*const.uc^4);
                    end
                    % fourth order polynom
                    if(strcmp(setting.fit_str,'poly4'))
                        p1 = fit1.p1;
                        p2 = fit1.p2;
                        p3 = fit1.p3;
                        p4 = fit1.p4;
                        p5 = fit1.p5;
                        beta1 = ((-vec.wlen.^2).*(3.*(vec.wlen.^2).*p1 + 2.*vec.wlen.*p2 + p3) + p5)./const.uc;
                        beta2 = ((vec.wlen.^3).*(3.*vec.wlen.*(2.*vec.wlen.*p1 + p2) + p3))./(pi*const.uc^2);
                        beta3 = -((3.*(vec.wlen.^4).*(10.*(vec.wlen.^2).*p1 + 4.*vec.wlen.*p2 + p3))./(2*(pi^2)*const.uc^3));
                        beta4 = (3.*(vec.wlen.^5).*(5.*vec.wlen.*(3.*vec.wlen.*p1 + p2) + p3))./((pi^3)*const.uc^4);
                    end
                    % third order polynom
                    if(strcmp(setting.fit_str,'poly3'))
                        p1 = fit1.p1;
                        p2 = fit1.p2;
                        p3 = fit1.p3;
                        p4 = fit1.p4;
                        beta1 = (-(vec.wlen.^2).*(2.*vec.wlen.*p1 + p2) + p4)./const.uc;
                        beta2 = ((vec.wlen.^3).*(3.*vec.wlen.*p1 + p2))./(pi*const.uc^2);
                        beta3 = -((3.*(vec.wlen.^4).*(4.*vec.wlen.*p1 + p2))./(2*(pi^2)*const.uc^3));
                        beta4 = (3.*(vec.wlen.^5).*(5.*vec.wlen.*p1 + p2))./((pi^3)*const.uc^4);
                    end
                    % second order polynom
                    if(strcmp(setting.fit_str,'poly2'))
                        p1 = fit1.p1;
                        p2 = fit1.p2;
                        p3 = fit1.p3;
                        beta1 = (-(vec.wlen.^2).*p1 + p3)./const.uc;
                        beta2 = ((vec.wlen.^3).*p1)./(pi*const.uc^2);
                        beta3 = -((3.*(vec.wlen.^4).*p1)./(2*(pi^2)*const.uc^3));
                        beta4 = (3.*(vec.wlen.^5).*p1)./((pi^3)*const.uc^4);
                    end

                    HTM.DIM(ii,jj).T(kk).O(nn).beta1 = beta1.*10^16;    % ps/cm
                    HTM.DIM(ii,jj).T(kk).O(nn).beta2 = beta2.*10^28;    % ps^2/cm
                    HTM.DIM(ii,jj).T(kk).O(nn).beta3 = beta3.*10^40;    % ps^3/cm
                    HTM.DIM(ii,jj).T(kk).O(nn).beta4 = beta4.*10^52;    % ps^4/cm

                    HTM.DIM(ii,jj).T(kk).O(nn).fit_Aeff     = fit(vec.wlen(disp.pte)',HTM.DIM(ii,jj).T(kk).O(nn).Aeff(disp.pte).*10^12,setting.fit_str);

                    HTM.DIM(ii,jj).T(kk).O(nn).fit_beta1    = fit(vec.wlen(disp.pte)',HTM.DIM(ii,jj).T(kk).O(nn).beta1(disp.pte)',setting.fit_str);
                    HTM.DIM(ii,jj).T(kk).O(nn).fit_beta2    = fit(vec.wlen(disp.pte)',HTM.DIM(ii,jj).T(kk).O(nn).beta2(disp.pte)',setting.fit_str);
                    HTM.DIM(ii,jj).T(kk).O(nn).fit_beta3    = fit(vec.wlen(disp.pte)',HTM.DIM(ii,jj).T(kk).O(nn).beta3(disp.pte)',setting.fit_str);
                    HTM.DIM(ii,jj).T(kk).O(nn).fit_beta4    = fit(vec.wlen(disp.pte)',HTM.DIM(ii,jj).T(kk).O(nn).beta4(disp.pte)',setting.fit_str);

                                                                        % is correct the 2 down here or was correct the .2 before?
                    HTM.DIM(ii,jj).T(kk).O(nn).Disp     = -beta2.*10^28.*2.*pi.*3e8./(vec.wlen.^2).*1e-4;
                    HTM.DIM(ii,jj).T(kk).O(nn).Gamma    = (2*pi*disp.n2)./(vec.wlen'.*1e-6.*HTM.DIM(ii,jj).T(kk).O(nn).Aeff);

                    HTM.DIM(ii,jj).T(kk).O(nn).Popti    = -((2*pi*3e8)*(beta2'.*10^28).*(1./(vec.wlen')-1/disp.lp).^2+(1/12)*(2*pi*3e8).*(beta4'.*10^52).*(1./(vec.wlen')-1/disp.lp).^4)./HTM.DIM(ii,jj).T(kk).O(nn).Gamma;
                    % it finds zero dispersion
                    disp.zr = fnzeros(spline(vec.wlen(disp.pte),HTM.DIM(ii,jj).T(kk).O(nn).beta2(disp.pte)));

                    if(isempty(disp.zr))
                        HTM.DIM(ii,jj).T(kk).O(nn).zero = 0;
                    else
                        HTM.DIM(ii,jj).T(kk).O(nn).zero = disp.zr(1,:); % all the points of zero dispersion
                    end
                end
            end
        end
    end
end

toc;
clear ii jj kk nn ans p1 p2 p3 p4 p5 p6 p7 beta1 beta2 beta3 beta4 fit1

%% SAVING
%{
%save data: lam, neff, beta2, beta4, Aeff, Gamma
%path=['D:\Documents\Cariplo2\images_wg\waveguide_angle',num2str(w),'.dat'];
path=['/images_wg/waveguide_angle',num2str(w),'.dat'];
export_data=[lam,HTE(1).W(1).O(1).neff,HTE(1).W(1).O(1).beta2,HTE(1).W(1).O(1).beta4,HTE(1).W(1).O(1).Aeff.*1e12,HTE(1).W(1).O(1).Gamma,HTM(1).W(1).O(1).neff,HTM(1).W(1).O(1).beta2,HTM(1).W(1).O(1).beta4,HTM(1).W(1).O(1).Aeff.*1e12,HTM(1).W(1).O(1).Gamma];

save(path,'export_data','-ascii');
%}

tic
fp = datestr(datetime,'yyyy-mm-dd_HH:MM:ss');
fp = strcat(setting.save_path, fp, setting.FILE_EXT);

fprintf('Saving data in %s\n',fp);
save(fp, 'HTE', 'HTM', 'vec', 'par');

clear fp ans
toc

%% END
fprintf('Data saved.\nProgram finished.');
