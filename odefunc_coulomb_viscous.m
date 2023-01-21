function dydt=odefunc_coulomb_viscous(t,x,M,miu_c,miu_v,K,v_b)
dydt=zeros(2,1);
x_act=v_b*t;
F_n=M*9.81; % normal force in newton
%F_f=miu_c*F_n* sign(x(2)) +(miu_v*x(2));
a=0.1
F_f=miu_c*F_n* tanh(a*x(2)) +(miu_v*x(2));
dydt(1)=x(2);
dydt(2)=(1/M)*((K*x_act)-(K*x(1))-F_f);
end
