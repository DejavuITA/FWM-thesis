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
fprintf('Reading preamble for options...\n');

%% INITIAL PARAMETERS

% setting up initial parameters
fprintf('Setting up initial parameters\n');

setting.comsolfile  = '../comsol/FWM_template.mph';    % filepath of the comsol model
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

%% MODEL LOADING

tic
fprintf('\nLoading the model...');
model = mphload(setting.comsolfile);
fprintf('Model Loaded.\n');
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
fprintf('Parameters Loaded.\n');
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
fp = datestr(datetime,'yyyy-MM-dd_HH:mm:ss');
fp = strcat(setting.save_path, fp, setting.FILE_EXT);
fprintf('Saving data in %s',fp);
fprintf('\n');
save(fp, 'HTE', 'HTM', 'vec', 'par');
toc
%%
fp = strcat(setting.save_path, 'HTE', setting.FILE_EXT);
save( fp, 'HTE' );
fp = strcat(setting.save_path, 'HTM', setting.FILE_EXT);
save( fp, 'HTM' );
fp = strcat(setting.save_path, 'vectors', setting.FILE_EXT);
save( fp, 'vec' );
fp = strcat(setting.save_path, 'parameters', setting.FILE_EXT);
save( fp, 'par' );

clear fp ans

fprintf('Data saved.\nProgram finished.');

%% nothing to do from this point on

%% PHASE MATCHING (CHI^2 SPONTANEOUS DOWN CONVERSION)
% it controls phase matching for every mode (?)
% (controllo phase matching modale)
%{
fprintf('Checking modal phase matching...');
tic;

pmatch.pump_wlen        = 1.55; % in [um]
pmatch.pump_wlen_half   = 1.55; % in [um]

pmatch.pump_freq        = (2*pi*const.uc)/(pmatch.pump_wlen);
pmatch.pump_freq_double = (2*pi*const.uc)/(pmatch.pump_wlen_half);

% optical parametric down conversion: f_idl = f_pump/2 + sigma; f_sig = f_pump/2-sigma
% Interval around f_pump/2
sigma_max = 30*1e-3; %in um
sigma_min = 0*1e-3;
sigma_step = 1*1e-3;
sigma = [sigma_min:sigma_step:sigma_max];
delta_sigma_min = 0.001*1e-3;
f_idler = 2*pi*const.uc./(lambda_pump_half+sigma);
f_signal = 2*pi*const.uc./(lambda_pump_half-sigma);

sol_down_conv.data = zeros(11,1);
sol_down_conv.walkoff_idler = 0;
sol_down_conv.walfoff_signal = 0;
sol_down_conv.walkoff = 0;

% per il walk off serve la FWHM impulso
dt=40; %ps
count = 0;
combs = 0;
n_pumps = zeros(1000,1);
    
    for ii=1:n_para_3
        for jj=1:n_para_2
            for pola_p=1:2
                for pola_i=1:2
                    for pola_s=1:2
                        for ord_p=1:n_modi
                            for ord_i=1:n_modi
                                for ord_s=1:n_modi
                                    %Check the mode existence
                                    count = count+1;
                                    combo_found = 0;
                                    
                                    if pola_p==1
                                        H_pump = HTE;
                                    else
                                        H_pump = HTM;
                                    end
                                    if pola_i==1
                                        H_idler = HTE;
                                    else
                                        H_idler = HTM;
                                    end
                                    if pola_p==1
                                        H_signal = HTE;
                                    else
                                        H_signal = HTM;
                                    end
                                    
                                    ex_pump = exist(H_pump(ii).W(jj).O(ord_p).fit_neff);
                                    ex_idler = exist(H_idler(ii).W(jj).O(ord_i).fit_neff);
                                    ex_signal = exist(H_signal(ii).W(jj).O(ord_s).fit_neff);
                                    ind_i_max = max(find(H_idler(ii).W(jj).O(kk_i).neff>0));
                                    ind_i_min = min(find(H_idler(ii).W(jj).O(kk_i).neff>0));
                                    ind_s_max = max(find(H_signal(ii).W(jj).O(kk_s).neff>0));
                                    ind_s_min = min(find(H_signal(ii).W(jj).O(kk_s).neff>0));
                                    
                                    if not(isempty(ind_i_max)) && not(isempty(ind_i_min))
                                        if lam(ind_i_max)>=lambda_pump_half+sigma(end) && lam(ind_i_min)<= lambda_pump_half-sigma(end)
                                            ok_idler = 1;
                                        else
                                            ok_idler = 0;
                                        end
                                    else
                                        ok_idler = 0;
                                    end
                                    
                                    if not(isempty(ind_s_max)) && not(isempty(ind_s_min))
                                        if lam(ind_s_max)>=lambda_pump_half+sigma(end) && lam(ind_s_min)<= lambda_pump_half-sigma(end)
                                            ok_signal = 1;
                                        else
                                            ok_signal = 0;
                                        end
                                    else
                                        ok_signal = 0;
                                    end
                                    if (ex_pump*ex_idler*ex_signal*ok_idler*ok_signal)
                                        %Now we are sure that the effective index has a meaning
                                        n_idler = H_idler(ii).W(jj).O(ord_i).fit_neff(lambda_pump_half+sigma);
                                        n_signal = H_signal(ii).W(jj).O(ord_s).fit_neff(lambda_pump_half-sigma);
                                        n_pump = H_pump(ii).W(jj).O(ord_p).fit_neff(lambda_pump);
                                        
                                        for pp=1:length(P0)
                                            
                                        phase_m = f_pump*n_pump-f_idler.*n_idler-f_signal.*n_signal;%+dkNL;
                                        zr = fnzeros(spline(sigma,phase_m),[sigma(1) sigma(end)]);
                                        if not(isempty(zr))
                                            for j=1:size(zr,2)
                                                if (zr(1,j)<delta_sigma_min)
                                                    riempi(j) = 1;
                                                else
                                                    riempi(j) = 0;
                                                end
                                            end
                                        end
                                        if (not(isempty(zr)) && sum(riempi)==length(riempi))
                                            %Walk off distance
                                            if combo_found==0
                                                combo_found = 1;
                                                comb = comb+1;
                                            end
                                            n_pumps(comb) = n_pumps(comb)+1;
                                            if pp==round(length(P0)/2)
                                            vg_p = 1/H_pump(ii).W(jj).O(ord_p).fit_beta1(lambda_pump);
                                            vg_i = 1/H_idler(ii).W(jj).O(ord_i).fit_beta1(lambda_pump_half+zr(1,:));
                                            vg_s = 1/H_signal(ii).W(jj).O(ord_s).fit_beta1(lambda_pump-zr(1,:));
                                            %Check the minimum distance
                                            walkoff_i(find(vg_p>=vg_i)) = vg_p*vg_p*dt./(abs(vg_p-vg_i(find(vg_p>=vg_i))));
                                            walkoff_i(find(vg_p<vg_i)) = vg_i(find(vg_p<vg_i)).*vg_p*dt./(abs(vg_p-vg_i(find(vg_p<vg_i))));
                                            walkoff_s(find(vg_p>=vg_s)) = vg_p*vg_p*dt./(abs(vg_p-vg_s(find(vg_p>=vg_s))));
                                            walkoff_s(find(vg_p<vg_s)) = vg_s(find(vg_p<vg_s)).*vg_p*dt./(abs(vg_p-vg_s(find(vg_p<vg_s))));
                                            %Store solution
                                            sol_down_conv.data = [sol_down_conv.data,[h(ii),w(jj),pola_p,ord_p,lambda_pump,pola_i, ...
                                                                  ord_i,zr(1,:)+lambda_pump_half,pola_s,ord_s,lambda_pump_half-zr(1,:)]'];
                                            sol_down_conv.walkoff_idler = [sol_down_conv.walkoff_idler,walkoff_i];
                                            sol_down_conv.walkoff_signal = [sol_down_conv.walkoff_signal,walkoff_s];
                                            sol_down_conv.walkoff = [sol_down_conv.walkoff,min(walkoff_i,walkoff_s)];
                                            clear walkoff_i;
                                            clear walkoff_s;
                                            end
                                        end
                                      end
                                    end
                                    if mod(count,1e4)==0
                                        fprintf('Combination: %d/%d\n',count,n_para_3*n_para_2*2*2*2*n_modi^3);
                                    end
                                end %for ord_s
                            end %for ord_i
                        end %for ord_p
                    end %for pola_s
                end %for pola_i
            end %for pola_p
        end %for n_para_2
    end %for n_para_3
    
n_pumps = n_pumps./max(n_pumps);
%Make the longest walkoff combinations the favourite ones
[~,pos] = sort(sol_down_conv.walkoff,'descend');
n_pumps = n_pumps(pos);
sol_down_conv.data = sol_down_conv.data(:,pos);
sol_down_conv.walkoff = sol_down_conv.walkoff(pos);
sol_down_conv.walkoff_idler = sol_down_conv.walkoff_idler(:,pos);
sol_down_conv.walkoff_signal = sol_down_conv.walkoff_signal(:,pos);
%printSFWM(sol_sFWM.data,sol_sFWM.walkoff_idler,sol_sFWM.walkoff_signal,sol_sFWM.cross,P0(round(length(P0)/2)),n_pumps);
    
toc
%}

%% SPONTANEOUS FOUR WAVE MIXING ZONE
%{
tic;

%phase mathing sFWM
lam_pump=1.55;% lambda pompa um
delta_lam_min=0.01e-3; %evita di trovare SFWM tra stessi modi a stessa lambda

%Chi3 tensor. Alpha, beta and gamma are the Euler's angle of rotation
%Set all equal to 0 if the waveguide is aligned with the z crystal.
%axis
alpha = 0;
beta = 0;
gamma = 0;
chi3 = chi3tensor(alpha,beta,gamma);

%limite ricerca generate
l1=1.501;%um
l2=1.599;%um
n_sample=100;


%silicon n2 chi3 real part in um^2/W
n2=(0.45*10^-17)*10^12;

%potenza picco di ingresso (W)
P0_max=500;
P0_min=0;
P0_step = 50;
P0 = [P0_min:P0_step:P0_max];
%P0 = 250;
vcc=VC*1e6; %um/s

o_p=2*pi*vcc/(lam_pump); %frequenza della pompa

lam_i=linspace(l1,l2,n_sample);

dt = 40; %pulse width in time (s)
o_i=2*pi*vcc./lam_i; %frequenza dell' idler
o_s = 2*o_p-o_i;
lam_s = 2*pi*vcc./o_s;

cost=n2*o_p; % self phase modulation contribution to dk


min_dk=10;
%mat_sol=zeros(n_sample,n_sample);
sol_sFWM.data=zeros(20,1);
sol_sFWM.walkoff_signal = zeros(5,1);
sol_sFWM.walkoff_idler = zeros(5,1);
sol_sFWM.walkoff = 0;
sol_sFWM.cross = zeros(3,1);

%  for j=1:length(P0)
%      for c=1:100
%      sFWM(j).combo(c).pol_p=0;
%      end
%  end

comb = 0;
count = 0;
n_pumps = zeros(1000,1);

%Modal phase matching
for ii=1:n_wg_hgt
    for jj=1:n_wg_wid
        
        for pol_p=2:2
            for pol_i=2:2
                for pol_s=2:2
                    
                    %for kk_p=1:n_modi %cicla sui possibili modi della pompa
                        %for kk_i=1:n_modi %cicla sui possibili modi dove puo finire l'idler
                            for kk_s=1:n_modi %cicla sui possibii modi dove puo finire il signal
                                %Check mode existence before
                                count = count+1;
                                combo_found = 0;
                                kk_p = kk_s;
                                kk_i = kk_s;
                                
                                if pol_p==1
                                    H_pump = HTE;
                                else
                                    H_pump = HTM;
                                end
                                if pol_i==1
                                    H_idler = HTE;
                                else
                                    H_idler = HTM;
                                end
                                if pol_s==1
                                    H_signal = HTE;
                                else
                                    H_signal = HTM;
                                end
                                
                                ex_pump = exist(H_pump(ii).W(jj).O(kk_p).fit_neff);
                                ex_signal = exist(H_signal(ii).W(jj).O(kk_s).fit_neff);
                                ex_idler = exist(H_idler(ii).W(jj).O(kk_i).fit_neff);
                                ind_p_max = max(find(H_pump(ii).W(jj).O(kk_p).neff>0));
                                ind_p_min = min(find(H_pump(ii).W(jj).O(kk_p).neff>0));
                                ind_i_max = max(find(H_idler(ii).W(jj).O(kk_i).neff>0));
                                ind_i_min = min(find(H_idler(ii).W(jj).O(kk_i).neff>0));
                                ind_s_max = max(find(H_signal(ii).W(jj).O(kk_s).neff>0));
                                ind_s_min = min(find(H_signal(ii).W(jj).O(kk_s).neff>0));
                                
                                if not(isempty(ind_p_max)) && not(isempty(ind_p_min))
                                    if lam(ind_p_max)>=l2 && lam(ind_p_min)<= l1
                                        ok_pump = 1;
                                    else
                                        ok_pump = 0;
                                    end
                                else
                                    ok_pump = 0;
                                end
                                
                                if not(isempty(ind_i_max)) && not(isempty(ind_i_min))
                                    if lam(ind_i_max)>=l2 && lam(ind_i_min)<= l1
                                        ok_idler = 1;
                                    else
                                        ok_idler = 0;
                                    end
                                else
                                    ok_idler = 0;
                                end
                                
                                if not(isempty(ind_s_max)) && not(isempty(ind_s_min))
                                    if lam(ind_s_max)>=l2 && lam(ind_s_min)<= l1
                                        ok_signal = 1;
                                    else
                                        ok_signal = 0;
                                    end
                                else
                                    ok_signal = 0;
                                end
                                
                                if (ex_pump*ex_signal*ex_idler*ok_pump*ok_idler*ok_signal)
                                    
                                    %Phase matching condition for idler wave
                                    for pp=1:length(P0)
                                    mat_sol= 2*cost*P0(pp)/H_pump(ii).W(jj).O(kk_p).fit_Aeff(lam_pump)-2*o_p*H_pump(ii).W(jj).O(kk_p).fit_neff(lam_pump)+o_i'.*H_idler(ii).W(jj).O(kk_i).fit_neff(lam_i) +(o_s)'.*H_signal(ii).W(jj).O(kk_s).fit_neff((lam_s));
                                    
                                         if (combo_found==0)
                                                 combo_found = 1;
                                                 comb = comb+1;
                                         end
                                        %n_pumps(comb) = n_pumps(comb)+1;  
                                        dk_integrated(comb,pp) = 1/(lam_i(end)-lam_i(1))*trapz(lam_i,abs(mat_sol));
                                        w_info(comb) = w(jj);
                                        ord(comb) = kk_s;
                                    end
                                end %close if existence
                                if mod(count,1e4)==0
                                    ffprintf('Combination: %d/%d\n',count,n_wg_hgt*n_wg_wid*2*2*2*n_modi^3);
                                end
                            %end %signal for
                        %end %idler for
                    end %pump for
                end %pol_s for
            end %pol_i for
        end %pol_p for
    end %n_wg_wid for
end %n_wg_hgt for

[vals,pos] = sort(mean(dk_integrated,2),'ascend');
w_info = w_info(pos);
ord = ord(pos);
fid = fopen('sFWMSelfTM.txt','w');
for j=1:length(vals)
    ffprintf(fid,'Integrated dk value: %e\n',vals(j));
    ffprintf(fid,'Mode order: %f\n',ord(j));
    ffprintf(fid,'Waveguide width: %f\n',w_info(j));
    ffprintf(fid,'--------------------------------------\n\n');
end
toc;

%}

%%%%
% saves data
%{
TEv=lam;
for ii=1:n_wg_hgt
    TEv=[0,lam']';
    for jj=1:n_wg_wid
        for kk=1:n_modi
            if(HTE(ii).W(jj).O(kk).neff(1)>0)
                TEv=[TEv,[kk,HTE(ii).W(jj).O(kk).neff']'];
            end
        end
    end
    eval(['TEv',num2str(ii),'=TEv;'])
end

TMv=lam;
for ii=1:n_wg_hgt
    TMv=[0,lam']';
    for jj=1:n_wg_wid
        for kk=1:n_modi
            if(HTM(ii).W(jj).O(kk).neff(1)>0)
                TMv=[TMv,[kk,HTM(ii).W(jj).O(kk).neff']'];
            end
        end
    end
    eval(['TMv',num2str(ii),'=TMv;'])
end
%}

%% COUPLING
% it selects the section and the order of pump guide and collector guide
% (selezione sezione guida di pompa e raccolta e relativo ordine)
%{
% ATTENZIONE - CONTROLLARE QUI IL PHASE MATCHING
% fnzeros(spilines)

%ordine pompa (SEMPRE 1)
o_p=1;
lambda_pump = 1.58-0.0001;
pos= find(lam>lambda_pump);
lamb = pos(1);

[o_p_x,o_p_y,o_p_fit]  = ord_para( lamb, HTE,'w', 1,o_p,'neff'); %Polarizzazione di entrata TE

%\martino I wish to have the interaction between all the modes. This
%means I have to make a row of fits...mh...maybe I should change
%o_a_fit into a vector of fits.
resample_w = 50;
ww=linspace(w(1),w(end),resample_w);
for(o_a=1:ord_max)
    [o_a_x,o_a_y,o_a_fit]  = ord_para( lamb, HTE,'w', 1,o_a,'neff');
    modes(o_a).TE.fit=o_a_fit;
    modes(o_a).TE.dom=o_a_x;
    modes(o_a).TE.dat=o_a_y;
    hold on;
    grid on;
    plot(modes(o_a).TE.dom,modes(o_a).TE.dat,'o-b');
end
xlim([min(modes(1).TM.dom) max(modes(1).TM.dom)]);
ylim([1.5 3]);
%}
    
%% interpolazione w
%{
    resample_w=1500;
    w_pump = 3.263158-0.001;
    pos_o = find(w>=w_pump);
    
    arrive_index = HTE(1).W(pos_o(1)).O(5).fit_neff(lambda_pump);
    
    ww=linspace(w(1),w(end),resample_w);
        
        %TE section
        yy_p_TE=modes(1).TE.fit(ww);
        yy_p_TE(find(yy_a_TE<1.5))=5; %Evita i valori dove il fit non ha piu senso
        yy_p_TE(find(yy_a_TE>3.5))=5; %Evita i valori dove il fit non ha piu senso
        [val,pos] = min(abs(yy_p_TE-arrive_index));
        ffprintf('\nDelta n: %f \t Width: %f',val,ww(pos));
        ffprintf('\nArrive width: %f\n',w(pos_o(1)));
%}

%% COUPLED MODES DISPERSION
% Find the possible coupled modes and compute dispersion (for DC) (IT NEEDS DISPERSE DISPERSION!)
%{
ord_ind_wg1 = 1;
ord_ind_wg2 = 1;
lambda_pump = 1.55;
for j=1:n_modi
    if exist(HTE(1).W(1).O(j).fit_neff)
        coupled_ord_wg1(ord_ind_wg1) = j;
        neff_WG1_temp{ord_ind_wg1} = coeffvalues(HTE(1).W(1).O(j).fit_neff);
        ord_ind_wg1 = ord_ind_wg1+1;
    end
    if exist(HTE(1).W(2).O(j).fit_neff)
        coupled_ord_wg2(ord_ind_wg2) = j;
        neff_WG2_temp{ord_ind_wg2} = coeffvalues(HTE(1).W(2).O(j).fit_neff);
        ord_ind_wg2 = ord_ind_wg2+1;
    end
end

ffprintf('\n\n----DIRECTIONAL COUPLER----\n\n');
ffprintf('Pump waveguide width (WG1) : %1.3f \t Slave waveguide width (WG2) : %1.3f\n',w(1),w(2));
ffprintf('WG1 mode order --- Effective index at %1.6f um\n',lambda_pump);

for j=1:length(coupled_ord_wg1)
    ffprintf('%d\t\t%1.5f\n',coupled_ord_wg1(j),HTE(1).W(1).O(coupled_ord_wg1(j)).fit_neff(lambda_pump));
end

ffprintf('\n\nWG2 mode order --- Effective index at %1.6f um\n',lambda_pump);

for j=1:length(coupled_ord_wg2)
    fprintf('%d\t\t%1.5f\n',coupled_ord_wg2(j),HTE(1).W(2).O(coupled_ord_wg2(j)).fit_neff(lambda_pump));
end
%}

%% MODE SELECTION
% it selects modes
%{
% if SELECT_COUPLED_MANUALLY
    
    wg1_cm_select = [1];
    wg2_cm_select = [1 2 3 4];
    neff_WG1{1} = HTE(1).W(1).O(1).neff(1);
    neff_WG2{1} = HTE(1).W(2).O(1).neff(1);
    neff_WG2{2} = HTE(1).W(2).O(2).neff(1);
    neff_WG2{3} = HTE(1).W(2).O(3).neff(1);
    neff_WG2{4} = HTE(1).W(2).O(4).neff(1);
%}

%% CICLA CICLI ?
%{
cicla_cicli %Execute scalar product (compute F matrix)

cicla_cicli_vettoriale %Compute M matrix
ffprintf('\nDONE\n\n');
display_FM(fc,mc,1.55);
%}

%% DIRECTIONAL COUPLER
%{
%function DC( M,F,neffWG1,neffWG2 )
M = {[38.660145,0],[0,0];[0,0],[45.724496 0]};
F = {[0,0],[11.044841 0];[11.947274 0],[0 0]};
neff_WG1 = mat2cell(2.6419);
neff_WG2 = mat2cell(2.6452);
DC(M,F,neff_WG1,neff_WG2,wg1_cm_select,wg2_cm_select);
%}

%% PRINT SECTIONS
%{
if(PRINT_SECTIONS)
    %run print_section.m
end
%}
