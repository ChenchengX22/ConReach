function R = exp_coord2rot(r,num_chart)
%GET_ROT 
Lambda = get_transition_matrix(num_chart);
R = expfun(skewfun(r));
R = Lambda*R;
end

