function [realMask,check_fitness,Complex_index,DATA_OUT] = Find_real_points(DATA_IN)

    %rowsWithComplex = find( any( imag(DATA_IN)~=0, 2 ) );
    realMask = ~any( imag(DATA_IN)~=0, 2 );
    DATA_OUT = DATA_IN(realMask, :); 
    Complex_index = double( any( imag(DATA_IN)~=0, 'all' ) );
    % if isempty(DATA_OUT)
    %     check_fitness = 0;
    % else
    %     check_fitness = 1;
    % end  \
    if Complex_index
        check_fitness = 0;
    else
        check_fitness = 1;
    end

end

