function [T, slot] = variable_to_tt(varName, C_1_tt, C_2_tt, C_3_tt)

    switch varName
        case 'a', T = C_1_tt; slot = 1;
        case 'b', T = C_1_tt; slot = 2;
        case 'c', T = C_1_tt; slot = 3;
        case 'd', T = C_2_tt; slot = 1;
        case 'e', T = C_2_tt; slot = 2;
        case 'f', T = C_2_tt; slot = 3;
        case 'g', T = C_3_tt; slot = 1;
        case 'h', T = C_3_tt; slot = 2;
        case 'i', T = C_3_tt; slot = 3;
        otherwise
            error('Unknown variable name.');
    end
end


