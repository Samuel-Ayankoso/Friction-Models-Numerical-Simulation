function dydt=odefunc_coulomb(t,x,M,miu_c,K,v_b)
dydt=zeros(2,1);

% I think this is not correct as the relative velocity
%v_rel=x(2)-v_b ;% the relative velocity
%F_f=miu_c*(K*x(1))*sign(v_rel)
%F_f=miu_c*F_n* sign(v_rel)
%F_c=1;F_f=F_c* sign(v_rel);
%dydt(2)=(1/M)*(F_f-(K*x(1)));

x_act=v_b*t;
F_n=M*9.81; % normal force in newton
a=0.1
F_f=miu_c*F_n* tanh(a*x(2));
dydt(1)=x(2);
dydt(2)=(1/M)*((K*x_act)-(K*x(1))-F_f);
end
