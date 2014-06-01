function current_target= f_whichisit(randnum, sum_prob)

%% A simple mapping function to map the probablity drawn to the target list
% by M. Baghaie July 2012

sum_prob = sum_prob.*1000;
randnum= randnum*1000;

current_target = 0;



if 0 <= randnum && randnum <= sum_prob(1) 
    
    current_target = 1;
    
else

for ii = 2:1:size(sum_prob,1)
    
    if sum_prob(ii-1) < randnum && randnum <= sum_prob(ii) 

        current_target = ii;
        break;

    end
    
end

if current_target == 0
    
    display('error in the f_whichisit function')
    %randnum;
    
end

end



