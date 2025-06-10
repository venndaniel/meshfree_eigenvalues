function be = cond_solve(F, f, iter)
    % solves F*b = f directly with a symmetric, diagonal pre-conditioner
    % rescales rows/columns to account for differences in magnitude;
    % results in a more stable solve
        rd = 1./sqrt(vecnorm(F, 2, 2));
        
        R = sparse(diag(rd)); % simple pre-conditioner
        
        Fp = R*F*R';
        Fp = 1/2*real(Fp + Fp'); % ensure symmetric
        if iter == 1
            be = R*(Fp\(R*f));
        else
            be = R*cond_solve(Fp, R*f, iter - 1);
        end
end