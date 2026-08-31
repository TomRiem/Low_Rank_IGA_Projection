function Num_pat = build_Q_entry_num_pattern_tt(entryName, C_1_tt, C_2_tt, C_3_tt, ops_num, tol)
% BUILD_Q_ENTRY_NUM_PATTERN_TT
%
% Directly builds the pushed numerator coefficient fields.
%
% Output:
%   Num_pat.(pat) is already in the numerator spline coefficient space,
%   i.e. after applying the corresponding L4 pushforward matrices.
%
% This avoids constructing the raw n^4 x n^4 x n^4 fourfold coefficient TT.

    minors = make_metric_minors();

    switch entryName
        case 'Q11'
            pairs = {'A','A'; 'D','D'; 'G','G'};
        case 'Q12'
            pairs = {'A','B'; 'D','E'; 'G','H'};
        case 'Q13'
            pairs = {'A','C'; 'D','F'; 'G','I'};
        case 'Q22'
            pairs = {'B','B'; 'E','E'; 'H','H'};
        case 'Q23'
            pairs = {'B','C'; 'E','F'; 'H','I'};
        case 'Q33'
            pairs = {'C','C'; 'F','F'; 'I','I'};
        otherwise
            error('Unknown entryName.');
    end

    Num_pat = struct('p12', [], 'p13', [], 'p23', [], 'm1', [], 'm2', [], 'm3', []);

    for a = 1:size(pairs,1)
        M1 = minors.(pairs{a,1});
        M2 = minors.(pairs{a,2});

        for i = 1:numel(M1)
            for j = 1:numel(M2)

                coeff = M1(i).coeff * M2(j).coeff;
                vars  = [M1(i).vars, M2(j).vars];

                [pat, T] = monomial_to_pushed_num_tt3(vars, C_1_tt, C_2_tt, C_3_tt, ops_num, tol);

                T = round(coeff * T, tol);

                if isempty(Num_pat.(pat))
                    Num_pat.(pat) = T;
                else
                    Num_pat.(pat) = round(Num_pat.(pat) + T, tol);
                end
            end
        end
    end
end


function [pattern_name, T] = monomial_to_pushed_num_tt3(vars, C_1_tt, C_2_tt, C_3_tt, ops_num, tol)
% MONOMIAL_TO_PUSHED_NUM_TT3
%
% Builds one monomial contribution directly in the pushed numerator space.
%
% It does NOT build the raw fourfold product TT.
% It directly constructs the TT cores after applying L4 in each dimension.

    slot_fac = { {}, {}, {} };

    for q = 1:numel(vars)
        [Tq, slot] = variable_to_tt(vars{q}, C_1_tt, C_2_tt, C_3_tt);
        slot_fac{slot}{end+1} = Tq;
    end

    counts = cellfun(@numel, slot_fac);

    switch sprintf('%d%d%d', counts(1), counts(2), counts(3))
        case '220'
            pattern_name = 'p12';
            fac = {slot_fac{1}{1}, slot_fac{1}{2}, slot_fac{2}{1}, slot_fac{2}{2}};

        case '202'
            pattern_name = 'p13';
            fac = {slot_fac{1}{1}, slot_fac{1}{2}, slot_fac{3}{1}, slot_fac{3}{2}};

        case '022'
            pattern_name = 'p23';
            fac = {slot_fac{2}{1}, slot_fac{2}{2}, slot_fac{3}{1}, slot_fac{3}{2}};

        case '211'
            pattern_name = 'm1';
            fac = {slot_fac{1}{1}, slot_fac{1}{2}, slot_fac{2}{1}, slot_fac{3}{1}};

        case '121'
            pattern_name = 'm2';
            fac = {slot_fac{1}{1}, slot_fac{2}{1}, slot_fac{2}{2}, slot_fac{3}{1}};

        case '112'
            pattern_name = 'm3';
            fac = {slot_fac{1}{1}, slot_fac{2}{1}, slot_fac{3}{1}, slot_fac{3}{2}};

        otherwise
            error('monomial_to_pushed_num_tt3: unsupported slot multiplicities [%d %d %d].', counts(1), counts(2), counts(3));
    end

    info = stiffness_pattern_info(pattern_name);
    ord_dim = info.orders;

    C = cell(3,1);

    for dim = 1:3
        ops1d = ops_num.(pattern_name){dim};
    
        C{dim} = fused_fourproduct_push_core_staged(fac, dim, ord_dim{dim}, ops1d);
    end

    T = cell2core(tt_tensor, C);
    T = round(T, tol);
end


function H = fused_fourproduct_push_core_staged(fac, dim, ord, ops1d)
% FUSED_FOURPRODUCT_PUSH_CORE_STAGED
%
% Builds one TT core of the pushed fourfold product directly, using the
% staged one-dimensional pushforward operators.
%
% This avoids:
%   1. building the raw n^4 physical mode,
%   2. applying the preassembled L4 matrix,
%   3. looping over all n^4 coefficient tuples.
%
% The TT rank order remains the original factor order fac{1},...,fac{4}.
% Only the physical factor order is changed by ord.

    Cfac = cell(4,1);
    G    = cell(4,1);

    lvec = zeros(1,4);
    rvec = zeros(1,4);
    nvec = zeros(1,4);

    A = cell(4,1);

    for k = 1:4
        Cfac{k} = core2cell(fac{k});
        G{k}    = Cfac{k}{dim};

        [lvec(k), nvec(k), rvec(k)] = size(G{k});

        % A{k} has rows = physical coefficient index,
        % columns = local TT-rank pair (left,right).
        %
        % Column order inside A{k}:
        %   p_k = sub2ind([lvec(k), rvec(k)], left_k, right_k)
        A{k} = reshape(permute(G{k}, [2 1 3]), nvec(k), lvec(k)*rvec(k));
    end

    if any(nvec ~= nvec(1))
        error('fused_fourproduct_push_core_staged: physical sizes do not match.');
    end

    % Reorder only the physical factor order expected by the 1D operator.
    % The rank order will be restored afterwards.
    Aord = A(ord);

    switch ops1d.kind
        case 'BBBB'
            % Physical order: B B B B.
            % First push factors 1 and 2 into S_{2p}, push factors 3 and 4
            % into S_{2p}, and then multiply the two S_{2p} factors.
            % No ordered n^2 or n_2^2 transfer matrices are assembled.
            Y12  = contract_same_space_pair_factorized(Aord{1}, Aord{2}, ops1d.Pbb, ops1d.IJbb);
            Y34  = contract_same_space_pair_factorized(Aord{3}, Aord{4}, ops1d.Pbb, ops1d.IJbb);
            Htmp = contract_same_space_pair_factorized(Y12, Y34, ops1d.P22, ops1d.IJ22);

        case 'DDBB'
            % Physical order: D D B B.
            % Apply the derivative map to the first two factors, push their
            % product into S_{2p-2}, push the last two basis factors into S_{2p},
            % and finally multiply the two intermediate spaces.
            Ydd  = contract_derivative_derivative_pair_factorized(Aord{1}, Aord{2}, ops1d.D, ops1d.Pdd, ops1d.IJdd);
            Ybb  = contract_same_space_pair_factorized(Aord{3}, Aord{4}, ops1d.Pbb, ops1d.IJbb);
            Htmp = contract_mixed_pair_factorized(Ydd, Ybb, ops1d.Rcoef, ops1d.AM);

        case 'DBBB'
            % Physical order: D B B B.
            % Apply the derivative map to factor 1, multiply it with factor 2,
            % push factors 3 and 4 into S_{2p}, and finally multiply the two
            % intermediate spaces.
            Ydb  = contract_derivative_basis_pair_factorized(Aord{1}, Aord{2}, ops1d.D, ops1d.Rdb, ops1d.AMdb);
            Ybb  = contract_same_space_pair_factorized(Aord{3}, Aord{4}, ops1d.Pbb, ops1d.IJbb);
            Htmp = contract_mixed_pair_factorized(Ydb, Ybb, ops1d.Rcoef, ops1d.AM);

        otherwise
            error('Unknown fourfold pushforward kind "%s".', ops1d.kind);
    end

    % Htmp has columns ordered by rank pairs in physical order ord.
    % Convert back to TT rank order fac{1},fac{2},fac{3},fac{4}, with
    % all left ranks first and all right ranks second.
    Hmat = reorder_staged_rank_columns_to_original_lr(Htmp, lvec, rvec, ord);

    nout = size(Hmat,1);
    rl   = prod(lvec);
    rr   = prod(rvec);

    H = reshape(Hmat, [nout, rl, rr]);
    H = permute(H, [2 1 3]);
end



function Y = contract_same_space_pair_factorized(A1, A2, Pcoef, IJ)
% CONTRACT_SAME_SPACE_PAIR_FACTORIZED
%
% Matrix-free application of the ordered same-space product operator.
%
% Inputs:
%   A1 : n x c1 coefficient/rank-slice matrix for factor 1
%   A2 : n x c2 coefficient/rank-slice matrix for factor 2
%   Pcoef(:,q) are product coefficients for the unordered overlap pair
%       IJ(q,:) = [i,j], i <= j.
%
% Output:
%   Y : nout x (c1*c2), with column ordering sub2ind([c1,c2],s1,s2).
%
% Algebraically this equals
%
%   build_ordered_same_space_operator(Pcoef,IJ,n) * kron(A2,A1),
%
% but avoids assembling the ordered n^2-column matrix.

    [n1, c1] = size(A1);
    [n2, c2] = size(A2);

    if n1 ~= n2
        error('contract_same_space_pair_factorized: A1 and A2 have incompatible row sizes.');
    end

    nout = size(Pcoef, 1);
    Y = zeros(nout, c1*c2);

    for q = 1:size(IJ,1)
        i = IJ(q,1);
        j = IJ(q,2);

        % Row of kron(A2,A1) corresponding to sub2ind([n,n],i,j).
        xij = kron(A2(j,:), A1(i,:));

        if i ~= j
            % Same product coefficients apply to the mirrored ordered pair.
            xij = xij + kron(A2(i,:), A1(j,:));
        end

        pcol = Pcoef(:,q);
        [rr, ~, vv] = find(pcol);
        if isempty(rr)
            continue;
        end

        Y(rr,:) = Y(rr,:) + full(vv) * xij;
    end
end


function Y = contract_derivative_derivative_pair_factorized(A1, A2, D, Pdd, IJdd)
% CONTRACT_DERIVATIVE_DERIVATIVE_PAIR_FACTORIZED
%
% Matrix-free application of the D-D pair transfer
%
%   beta_i' beta_j' -> S_{2p-2}.
%
% This equals the action of the formerly assembled Pdd_ord on kron(A2,A1).

    DA1 = D * A1;
    DA2 = D * A2;
    Y   = contract_same_space_pair_factorized(DA1, DA2, Pdd, IJdd);
end


function Y = contract_derivative_basis_pair_factorized(A1, A2, D, Rdb, AMdb)
% CONTRACT_DERIVATIVE_BASIS_PAIR_FACTORIZED
%
% Matrix-free application of the D-B product transfer
%
%   beta_i' beta_j -> S_{2p-1}.
%
% This equals the action of the formerly assembled Rdb_ord on kron(A2,A1).

    DA1 = D * A1;
    Y   = contract_mixed_pair_factorized(DA1, A2, Rdb, AMdb);
end


function Y = contract_mixed_pair_factorized(A1, A2, Rcoef, AM)
% CONTRACT_MIXED_PAIR_FACTORIZED
%
% Matrix-free application of a mixed product transfer between two different
% spline spaces.
%
% Inputs:
%   A1 : nA x c1
%   A2 : nB x c2
%   Rcoef(:,q) are product coefficients for AM(q,:) = [a,b].
%
% Output:
%   Y : nout x (c1*c2), with column ordering sub2ind([c1,c2],s1,s2).
%
% Algebraically this equals
%
%   build_ordered_mixed_operator(Rcoef,AM,nA,nB) * kron(A2,A1),
%
% but avoids assembling the ordered nA*nB-column matrix.

    [nA, c1] = size(A1);
    [nB, c2] = size(A2);

    nout = size(Rcoef, 1);
    Y = zeros(nout, c1*c2);

    for q = 1:size(AM,1)
        a = AM(q,1);
        b = AM(q,2);

        if a < 1 || a > nA || b < 1 || b > nB
            error('contract_mixed_pair_factorized: AM contains an out-of-range index.');
        end

        xab = kron(A2(b,:), A1(a,:));

        rcol = Rcoef(:,q);
        [rr, ~, vv] = find(rcol);
        if isempty(rr)
            continue;
        end

        Y(rr,:) = Y(rr,:) + full(vv) * xab;
    end
end


function Hmat = reorder_staged_rank_columns_to_original_lr(Htmp, lvec, rvec, ord)
% REORDER_STAGED_RANK_COLUMNS_TO_ORIGINAL_LR
%
% Htmp columns are ordered by rank-pairs in the physical order ord:
%
%   (left_ord1, right_ord1,
%    left_ord2, right_ord2,
%    left_ord3, right_ord3,
%    left_ord4, right_ord4)
%
% We need columns ordered as a TT core:
%
%   (left_1,left_2,left_3,left_4,
%    right_1,right_2,right_3,right_4)
%
% in the original factor order fac{1},...,fac{4}.

    dims_ord = zeros(1,8);

    for t = 1:4
        k = ord(t);
        dims_ord(2*t-1) = lvec(k);
        dims_ord(2*t)   = rvec(k);
    end

    idx = reshape(1:prod(lvec .* rvec), dims_ord);

    perm = zeros(1,8);

    % Positions of left ranks of original factors 1,2,3,4
    for k = 1:4
        t = find(ord == k, 1);
        perm(k) = 2*t - 1;
    end

    % Positions of right ranks of original factors 1,2,3,4
    for k = 1:4
        t = find(ord == k, 1);
        perm(4+k) = 2*t;
    end

    old_of_new = permute(idx, perm);
    old_of_new = old_of_new(:);

    Hmat = Htmp(:, old_of_new);
end