function psi_value = get_psi(R, Rd)
    psi_value = 0.5*trace(eye(3) - Rd'*R);
    % % psi__Rd = eye(3) - Rd'*Rd;
    % fprintf("psi__Rd = %.6f \n", norm(psi__Rd,2));
    % 
    % psi_R = eye(3) - R'*R;
    % fprintf("psi_R = %.6f \n", norm(psi_R,2));    
    % 
    % psi_RRd = Rd - R;
    % fprintf("psi_RRd = %.6f \n", norm(psi_RRd,2));    

end