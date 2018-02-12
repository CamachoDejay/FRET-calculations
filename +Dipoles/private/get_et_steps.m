function [ et_steps ] = get_et_steps( no_ET_prob )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

test = true;
et_steps = 0;
while test
    
    n = 100;
    tmp = rand(n,1);
    tt = find(tmp < no_ET_prob,1);
   
    if isempty(tt)
        et_steps = et_steps + n;
    else
        et_steps = et_steps + tt;
        test = false; 

    end
    
   if et_steps >= 1500
       test = false;
       disp('manual stop')
   end

end
end

