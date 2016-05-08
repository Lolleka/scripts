function [ val ] = Gauss2D(x,y,x0,y0,sigx,sigy)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    val = gaussmf(x,[sigx x0])*gaussmf(y,[sigy y0]);
end

