function Lambda = get_transition_matrix(chart_num)
%GET_TRANSITION_MATRIX 
switch chart_num
    case 0
        Lambda = [1,0,0;0,1,0;0,0,1];
    case 1
        Lambda = [1,0,0;0,-1,0;0,0,-1];
    case 2
        Lambda = [-1,0,0;0,1,0;0,0,-1];
    case 3
        Lambda = [-1,0,0;0,-1,0;0,0,1];
    otherwise
        fprintf("Err_wrong_chart_num.\n");
end

end

