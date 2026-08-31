function info = stiffness_pattern_info(pat)
% STIFFNESS_PATTERN_INFO
% One-dimensional transfer kinds, local factor ordering, and required
% derivative-direction indices for one numerator pattern.

    switch pat
        case 'p12'
            info.kinds = {'DDBB','DDBB','BBBB'};
            info.orders = {[1 2 3 4], [3 4 1 2], [1 2 3 4]};
            info.needed_q = 9;

        case 'p13'
            info.kinds = {'DDBB','BBBB','DDBB'};
            info.orders = {[1 2 3 4], [1 2 3 4], [3 4 1 2]};
            info.needed_q = 5;

        case 'p23'
            info.kinds = {'BBBB','DDBB','DDBB'};
            info.orders = {[1 2 3 4], [1 2 3 4], [3 4 1 2]};
            info.needed_q = 1;

        case 'm1'
            info.kinds = {'DDBB','DBBB','DBBB'};
            info.orders = {[1 2 3 4], [3 1 2 4], [4 1 2 3]};
            info.needed_q = [6 8];

        case 'm2'
            info.kinds = {'DBBB','DDBB','DBBB'};
            info.orders = {[1 2 3 4], [2 3 1 4], [4 1 2 3]};
            info.needed_q = [3 7];

        case 'm3'
            info.kinds = {'DBBB','DBBB','DDBB'};
            info.orders = {[1 2 3 4], [2 1 3 4], [3 4 1 2]};
            info.needed_q = [2 4];

        otherwise
            error('Unknown stiffness numerator pattern "%s".', pat);
    end
end