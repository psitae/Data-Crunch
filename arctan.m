function [ angle ] = arctan( a, b )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
first = atan(b/a);
if a<0
    angle = first + pi;
else
    angle = first;    
end

end

