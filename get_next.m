function [gap] = get_next(Mu,dst)
r = -log(rand());
if r > 11
    gap = 11/Mu;
    return;
end
a = 0;
b = 11/Mu;
while b-a > 1e-04
    t = Mu*(a+b)/2 + dst - dst*exp(-(a+b)/2);
    if t>r
        b=(a+b)/2;
    else
        a = (a+b)/2;
    end
end
gap = a;
end

