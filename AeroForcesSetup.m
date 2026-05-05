function [cl0,cl1,cd0,cd1,cd2] = AeroForcesSetup(v_stall, v_ldmax, ld_max, weight, S, rho)

% lift calculation
cl_stall = 2*weight/(rho*S*v_stall^2);

cl0 = 0;

cl1 = cl_stall;

% Drag Calculation
cl_ldmax = 2*weight/(rho*S*v_ldmax^2);

a_ldmax = (cl_ldmax - cl0)/cl1;

k = (cl1*cl_ldmax)/(ld_max*(cl1^3*a_ldmax^2 + 2*cl0*cl1^2*a_ldmax + cl0^2*cl1 +  cl1*cl_ldmax^2));

Cd0 = (cl_ldmax/ld_max) - k*cl_ldmax^2;
cd0 = Cd0 + k*cl0^2;
cd1 = 2*k*cl0*cl1;
cd2 = k*cl1^2;

end