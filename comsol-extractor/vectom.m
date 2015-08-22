function [ M ] = vectom( V,n,m )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% This function should convert a n*m-long vector into a n time m matrix.
if(numel(V)==n*m)
    %printf('Checked dimensions'); 
    for(i=1:numel(V))
        j=idivide(i-1, uint16(n), 'floor')+1;
        M(mod(i-1,n)+1,j)=V(i);
    end
    M;
else
fprintf('Wrong dimensions in vectom!');    
end

end

