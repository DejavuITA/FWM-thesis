function [ output ] = wavenumber( TE, wid, hgt, order, wlen, temp, vec, HTE, HTM, TOC )
%wavenumber: It calculates the wavenumber with dependence of wavelength, order
%of the mode and temperature, interpolating data from HTE or HTM structures.
%   wid, hgt and order must be integer (they are considered indexes).
%   wlen, temp can be any real number included the interval considered.
%   other parameters are vectors and structures required to calculate

c     = 299792458;  % [m/s]
cps   = c*1e-12;	% [m/ps]

[value.wlen, index.wlen] = min( abs(vec.wlen - wlen) );
[value.temp, index.temp] = min( abs(vec.temp - temp) );

wlen0 = vec.wlen(index.wlen);
Omega  = wlen - wlen0;
Fmega  = cps/wlen - cps/wlen0;  % [1/ps]

% physical units are correct ??
% betas are in [ps^x/cm] for beta_x, therefore beta_x are in 1e-10*[s^x/m]
% c/Omega is in [m/s/m]=[1/s] --> (c/Omega)^x*beta_x = [s^x/m/s^x] = [1/m]
% alpha is in [1/T] --> alpha*T = [1]
if ( TE == 1 )
    fprintf('%d,\t', index.wlen );
    tmp =       2*pi/wlen0     * HTE.DIM(wid, hgt).T(1).O(order).neff(index.wlen);           % [1/m] ~1e7
    tmp = tmp + (Fmega)        * HTE.DIM(wid, hgt).T(1).O(order).beta1(index.wlen)*100;      % [1/m] ~1e6
    tmp = tmp + (Fmega)^2 / 2  * HTE.DIM(wid, hgt).T(1).O(order).beta2(index.wlen)*100;      % [1/m] ~1e4
    tmp = tmp + (Fmega)^3 / 6  * HTE.DIM(wid, hgt).T(1).O(order).beta3(index.wlen)*100;      % [1/m] ~1e3
    tmp = tmp + (Fmega)^4 / 24 * HTE.DIM(wid, hgt).T(1).O(order).beta4(index.wlen)*100;      % [1/m] ~1e2
    tmp = tmp + 2*pi/wlen0     * TOC.TE.DIM(wid, hgt).O(order).WL(index.wlen).alpha * temp;  % [1/m] ~1e5
    output = tmp;
end
if ( TE == 0 )
    tmp =       2*pi/wlen0     * HTM.DIM(wid, hgt).T(1).O(order).neff(index.wlen);           % [1/m] ~1e7
    tmp = tmp + (Fmega)        * HTM.DIM(wid, hgt).T(1).O(order).beta1(index.wlen)*100;      % [1/m] ~1e6
    tmp = tmp + (Fmega)^2 / 2  * HTM.DIM(wid, hgt).T(1).O(order).beta2(index.wlen)*100;      % [1/m] ~1e4
    tmp = tmp + (Fmega)^3 / 6  * HTM.DIM(wid, hgt).T(1).O(order).beta3(index.wlen)*100;      % [1/m] ~1e3
    tmp = tmp + (Fmega)^4 / 24 * HTM.DIM(wid, hgt).T(1).O(order).beta4(index.wlen)*100;      % [1/m] ~1e2
    tmp = tmp + 2*pi/wlen0     * TOC.TM.DIM(wid, hgt).O(order).WL(index.wlen).alpha * temp;  % [1/m] ~1e5
    
    output = tmp;
end

