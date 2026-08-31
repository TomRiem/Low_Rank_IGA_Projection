function span = find_span(U, p, x)
% Locate the active knot span using 1-based MATLAB indexing.

    n = numel(U) - p - 1;

    if x >= U(n+1)
        span = n;
        return;
    end

    low = p+1;
    high = n+1;
    span = floor((low+high)/2);

    while x < U(span) || x >= U(span+1)
        if x < U(span)
            high = span;
        else
            low = span;
        end
        span = floor((low+high)/2);
    end
end