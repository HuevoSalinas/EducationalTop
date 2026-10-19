function [xnew,opt] = TOPupdateDV(x,opt)

tol = 1e-6;

opt.dJ = OCshift_dJdx(opt.dJ); % Check if sensitiv are all negative

l2 = 1;
l1 = 0;

xnew = OCupdate_embedded(x,opt,l2);
g = opt.v_func(xnew, opt) - opt.volfrac;
while (g > 0)
    l1 = l2;
    l2 = 10 * l2;
    xnew = OCupdate_embedded(x,opt,l2);
    g = opt.v_func(xnew, opt) - opt.volfrac;
end
while ((l2-l1)/(l1+l2) > tol)
    lmid = (l2+l1)/2;
    xnew = OCupdate_embedded(x,opt,lmid);
    g = opt.v_func(xnew, opt) - opt.volfrac;
    if (g > 0)
        l1 = lmid;
    else
        l2 = lmid;
    end
end
lmid = (l2+l1)/2;
xnew = OCupdate_embedded(x,opt,lmid);
     
end


function [xnew] = OCupdate_embedded(x,opt,lambda)
B = (-opt.dJ./opt.dv/lambda).^opt.eta;
xnew = max(opt.xLB,min(opt.xUB,x.*B));
xnew(opt.pass_solid) = 1;
xnew(opt.pass_void) = opt.xmin;
end

function [dJ] = OCshift_dJdx(dJ)
max_dJ = max(dJ);
if (max_dJ>0)
    if (max_dJ>1e-12)
        fprintf('INFO -- OC found a positive sensitivity. max_dJ = %g\n (shifting scale).\n',max_dJ)
    end
    dJ = dJ - max_dJ;
end
end
