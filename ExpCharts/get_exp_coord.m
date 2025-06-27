function [flag_valid,cur_chart_num,r] = get_exp_coord(R,chart_num,flag_sw_chart)
% GET_EXP_COORD 
% when flag_sw_chart == true it helps to find the best chart
% nearest to the origin points
if flag_sw_chart
    r_map = zeros(3,4);
    dist = zeros(1,4);
    for i=1:4
        Lambda = get_transition_matrix(i-1);
        R1 = Lambda*R;
        % norm_r = acos((trace(R1)-1)/2);
        [hat_r,~] = logm(R1);
        r = veemap(hat_r);
        r_map(:,i) = r;
        dist(i) = norm(r);
    end
    num = find(dist == min(dist,[],'all'));
    cur_chart_num = num-1;
    r = r_map(:,num);
    flag_valid = true;
else
    % One may choose to display it in a chart, even if the chart is not valid.
    Lambda = get_transition_matrix(chart_num);
    R1 = Lambda*R;
    % norm_r = acos((trace(R1)-1)/2);
    [hat_r,flag_sw] = logm(R1);
    flag_valid = ~flag_sw;
    r = veemap(hat_r);
    cur_chart_num = chart_num;
end
end





