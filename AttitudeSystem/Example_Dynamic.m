function [Next_r,Next_Ang_Vel] = Example_Dynamic(Cur_r,Cur_Ang_Vel,time_step,A)
% vectors (3,1)
norm_r = norm(Cur_r);
if norm_r == 0
    dotr = Cur_Ang_Vel;
else
    
    dotr = ((eye(3) ...
            +skewfun(Cur_r)/2 ...
            +skewfun(Cur_r)*skewfun(Cur_r)/(norm_r*norm_r)...
            -sin(norm_r)*skewfun(Cur_r)*skewfun(Cur_r)/(2*norm_r*(1-cos(norm_r))))...
            )*Cur_Ang_Vel;
end

Next_r = Cur_r + dotr*time_step;
dotw = A*Cur_Ang_Vel;
Next_Ang_Vel = Cur_Ang_Vel + dotw*time_step;

end
