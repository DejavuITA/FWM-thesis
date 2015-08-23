function [ output_args ] = get_par( model, parameter, data_set)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

tmp = mphglobal(model,{parameter},'dataset', data_set,'outersolnum','all');
tmp = tmp(1,:);
tmp = sort(tmp);
tmp = unique(tmp);
output_args = tmp;

clear ans tmp
end