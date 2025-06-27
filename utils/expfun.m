function R = expfun(hatr)
%EXPFUNC 
r = veemap(hatr);
if norm(r) == 0
    R = eye(3);
else
    k = (1-cos(norm(r)))/(norm(r)*norm(r));
    R = eye(3)+sin(norm(r))*hatr/norm(r)+k*hatr*hatr;
end
end

