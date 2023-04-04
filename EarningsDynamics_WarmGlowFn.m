function warmglow=EarningsDynamics_WarmGlowFn(aprime,a,z,upsilon,epsilon,alpha,kappabeta,agej,Jr,warmglow1,warmglow2,warmglow3,beta,sj)
% Earning dynamics process following GKOS2021

warmglow=0;
% add the warm glow to the return, but only near end of life
if agej>=Jr+10
    % Warm glow of bequests
    warmglow=warmglow1*((aprime-warmglow2)^(1-warmglow3))/(1-warmglow3);
    % Modify for beta and sj (get the warm glow next period if die)
    warmglow=beta*(1-sj)*warmglow;
    % add the warm glow to the return
%     if c>0 % I don't think this should be needed, but added to be sure
%         F=F+warmglow;
%     end
end

end
