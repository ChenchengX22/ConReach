% function matrix_type = Check_Matrix_Definiteness(A)
function if_pd = Check_Matrix_Definiteness(A)
    if ~isequal(A, A')
        error('matrix should be symmetric');
    end

    eigenvalues = eig(A);
    if_pd = false;

    num_positive = sum(eigenvalues > 0);
    num_negative = sum(eigenvalues < 0);
    num_zero = sum(eigenvalues == 0);
    n = length(eigenvalues);

    if num_positive == n
        matrix_type = 'positive_definite';
        if_pd = true;
    elseif num_positive + num_zero == n
        matrix_type = 'positive_semidefinite';
        if_pd = false;
    elseif num_negative == n
        matrix_type = 'negative_definite';
        if_pd = false;
    else
        matrix_type = 'indefinite';
        if_pd = false;
    end
end
