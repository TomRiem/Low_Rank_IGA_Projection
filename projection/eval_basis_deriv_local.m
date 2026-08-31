function [I, dB] = eval_basis_deriv_local(U, p, xa, xb, xq)
% EVAL_BASIS_DERIV_LOCAL  Evaluate first derivatives of the (p+1) active
% B-spline basis functions on one element span [xa, xb].
%
% Output
%   I   : 1 x (p+1) global indices of the active basis functions
%   dB  : nq x (p+1) derivatives

    xmid = 0.5*(xa + xb);
    span = find_span(U, p, xmid);
    I    = (span-p):span;

    nq = numel(xq);
    dB = zeros(nq, p+1);

    if p == 0
        return
    end

    % basis_funs(span, x, p-1, U) returns p values for
    % N_{span-(p-1)} ... N_{span}, i.e. global indices span-p+1 .. span.
    % The degree-p active window starts one earlier at span-p, so
    % N_{span-p, p-1} is missing. It is identically zero on this span
    % (its support ends exactly at U(span-p+p) = U(span)), so we
    % place it as a zero in column 1 and shift the rest into columns 2..p+1.
    %
    % After the shift:
    %   B_prev(:, k)  corresponds to global index span-p + k - 1  (k = 1..p+1)
    %
    % The derivative formula for the k-th degree-p function
    % (global i = span-p+k, 1-based k = 1..p+1):
    %   first  term: N_{i, p-1}  = B_prev(:, k)     valid for k = 1..p+1
    %   second term: N_{i+1, p-1}  = B_prev(:, k+1)   valid for k = 1..p

    B_prev = zeros(nq, p+1);          % col 1 stays zero (missing left function)
    for q = 1:nq
        B_prev(q, 2:p+1) = basis_funs(span, xq(q), p-1, U);
    end

    for k = 1:p+1
        i = span - p + k - 1;   % Correct: i = I(k)
    
        % First term:
        % + p/(U(i+p)-U(i)) * N_{i,p-1}
        denom1 = U(i+p) - U(i);
        if denom1 > 1e-14
            dB(:,k) = dB(:,k) + (p / denom1) * B_prev(:,k);
        end
    
        % Second term:
        % - p/(U(i+p+1)-U(i+1)) * N_{i+1,p-1}
        if k <= p
            denom2 = U(i+p+1) - U(i+1);
            if denom2 > 1e-14
                dB(:,k) = dB(:,k) - (p / denom2) * B_prev(:,k+1);
            end
        end
    end
end










