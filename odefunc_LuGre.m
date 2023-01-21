function[xdot, zdot, F_f] = odefunc_LuGre (t, q, M,K,v_b, Fs, Fc, sigma_0, sigma_1, sigma_2, vs)
    y = v_b*t;              % 0.1 m/s
    u = K * (y - q(1,:));   % force by the spring
   
    zdot = q(2,:) - ( (q(3,:).*abs(q(2,:))*sigma_0) ./ ...
           (Fc+(Fs-Fc)*exp(-(q(2,:)/vs).^2)) );
       
    F_f = sigma_0*q(3,:) + sigma_1 * zdot + sigma_2*q(2,:);

    qdot_1 = q(2,:);
    qdot_2 = (u - F_f) / M;
    qdot_3 = zdot;
    xdot = [qdot_1 ; qdot_2; qdot_3 ];

end