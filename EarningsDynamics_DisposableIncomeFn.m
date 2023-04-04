function DisposableIncome=EarningsDynamics_DisposableIncomeFn(aprime,a,z,upsilon,epsilon,alpha,kappabeta,w,agej,Jr,pension,incomefloor,r,kappa_j,eta1,eta2)
% Earning dynamics process of Model 6 of GKOS2021

DisposableIncome=0;
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
    DisposableIncome=AfterTaxIncome;
else % Retirement
    Income=r*a;
    if Income>0
        IncomeTax=eta1+eta2*log(Income)*Income;
    else
        IncomeTax=0;
    end
    % Income floor is not relevant as all get pension (and pension>incomefloor)
    DisposableIncome=pension+Income-IncomeTax;
end

end
