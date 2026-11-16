function [list] = pickMitagationWindows(se,skipSome, endSkip)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
arguments
se;
skipSome = 0; % add indeces of this number of elements from the beginning without checking them in forceFit
endSkip = 0; % add indeces of this number of elements from the end without checking them in forceFit
end

ff = se.forceFit;

len = size(ff,1);
index = 1;
list(index) = 1; % always include the index of the first element

for i = 2:1+skipSome
       index = index+1;
       list(index) = i;
end

for i = 1+skipSome:len-endSkip-1
   if ff(i)==1
       index = index+1;
       list(index) = i;
   end
end

for i = len-endSkip:len-1
      index = index+1;
     list(index) = i;
end

list(index+1) = len;

end 
