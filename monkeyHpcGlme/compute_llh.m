function [log_llh,log_llh_model] = compute_llh(A,n,param)
%original function that came with the LN Model code

r = exp(A * param); 
n = reshape(n,numel(n),1); 
meanFR = mean(n);

log_llh_model = sum(r-n.*log(r)+log(factorial(n)))/sum(n);
log_llh_mean = sum(meanFR-n.*log(meanFR)+log(factorial(n)))/sum(n);
log_llh = (-log_llh_model + log_llh_mean);
log_llh = log_llh/log(2); %# convert from nats to bits

end
