function f = package_envelope_s0_3drug(c0,surv1,dose1,surv2,dose2,surv3,dose3)

% c = [lam1, h1, s0, lam2, h2]
function S = sub3(c)
y1 = (1-c(3))./((dose1/c(1)).^c(2) + 1) + c(3);
y2 = (1-c(3))./((dose2/c(4)).^c(5) + 1) + c(3);
y3 = (1-c(3))./((dose3/c(6)).^c(7) + 1) + c(3);
S = ((norm(y1-surv1))^2+(norm(y2-surv2))^2+(norm(y3-surv3))^2);
end

options = optimset('TolX',1e-8,'TolFun',1e-8,'MaxIter',1e4,'MaxFunEval',1e4,'Display','off');
f = fmincon(@sub3,c0,[],[],[],[],[0 0 0 0 0 0 0],[Inf 4 1 Inf 4 Inf 4],[],options); %Hill coefficients and EC50's cannot be negative

end