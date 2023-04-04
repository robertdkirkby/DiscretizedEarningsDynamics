function F=EarningsDynamics_ReturnFn(aprime,a,z,upsilon,epsilon,alpha,kappabeta,w,sigma,agej,Jr,pension,incomefloor,r,kappa_j,warmglow1,warmglow2,warmglow3,beta,sj,eta1,eta2,CEV)
% Earning dynamics process following GKOS2021

F=-Inf;
if agej<Jr % If working age
    Income=w*(1-upsilon)*exp(kappa_j+alpha+kappabeta+z+epsilon)+r*a;
    if Income>0
        IncomeTax=eta1+eta2*log(Income)*Income;
    else
        IncomeTax=0;
    end
    AfterTaxIncome=Income-IncomeTax;
    if AfterTaxIncome<incomefloor
        AfterTaxIncome=incomefloor;
    end
    c=AfterTaxIncome+a-aprime;
else % Retirement
    Income=r*a;
    if Income>0
        IncomeTax=eta1+eta2*log(Income)*Income;
    else
        IncomeTax=0;
    end
    % Income floor is not relevant as all get pension (and pension>incomefloor)
    c=pension+(Income-IncomeTax)+a-aprime;
end

if c>0
    c=c*(1+CEV); % CEV is typically zero, is used to compute consumption equivalent variation
    F=(c^(1-sigma))/(1-sigma); % The utility function
end

% add the warm glow to the return, but only near end of life
if agej>=Jr+10
    % Warm glow of bequests
    warmglow=warmglow1*((aprime-warmglow2)^(1-warmglow3))/(1-warmglow3);
    % Modify for beta and sj (get the warm glow next period if die)
    warmglow=beta*(1-sj)*warmglow;
    % add the warm glow to the return
    if c>0 % I don't think this should be needed, but added to be sure
        F=F+warmglow;
    end
end

end
