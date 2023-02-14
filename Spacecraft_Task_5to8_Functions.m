
% function for spacecraft tasks course 
format long;
%% TASK 5
computeRcN = RcN_t(330);
% disp(computeRcN);

wRcN = omegaRcN(330);
% disp(wRcN);

%% TASK 6
%   sigmaBN: MRP spacecraft’s initial attitude 
%   b_omegaBN: initial body angular velocity
%   N_omegaRN: calculated angular velocity Rx-frame wrt N
%   dcmRN: Calculated DCM of Rx-frame wrt N
t = 0; % initial time t0 
sigmaBN = [(0.3) (-0.4) (0.5)].';
b_omegaBN = [deg2rad(1.00) deg2rad(1.75) deg2rad(-2.20)].';
dcmRN_Power = [-1 0 0; 0 0 1; 0 1 0];%RsN
n_omegaRN_Power = [0 0 0].';
dcmRN_Nadir = [-0.2133    -0.8758   -0.4330;
   -0.9619    0.1107    0.2500;
   -0.1710    0.4698   -0.8660]; %RnN
n_omegaRN_Nadir = [0.0001513 -0.0004157 0.0007663].';
dcmRN_Comm = RcN_t(0); %RcN
n_omegaRN_Comm = omegaRcN(0);


[OrientError_PowerMode, AngularVelocityError_PowerMode] = Errors(t, sigmaBN, b_omegaBN, dcmRN_Power, n_omegaRN_Power);
disp("OrientError_PowerMode");
disp(OrientError_PowerMode);
disp("AngularVelocityError_PowerMode");
disp(AngularVelocityError_PowerMode);

[OrientError_NadirPointing, AngularVelocityError_NadirPointing] = Errors(t, sigmaBN, b_omegaBN, dcmRN_Nadir, n_omegaRN_Nadir);
disp("OrientErrors_NadirPointing");
disp(OrientError_NadirPointing);
disp("AngularVelocityError_NadirPointing");
disp(AngularVelocityError_NadirPointing);

[OrientErrors_CommMode, AngularVelocityError_CommMode]  = Errors(t, sigmaBN, b_omegaBN, dcmRN_Comm, n_omegaRN_Comm);
disp("OrientErrors_CommMode");
disp(OrientErrors_CommMode);
disp("AngularVelocityError_CommMode");
disp(AngularVelocityError_CommMode);

%% Task 7 -- UPDATE

% Integrate State vector X forward 500s with u = 0
% Provide H = [I]*b_omegaBN at 500s expressed in the B frame

[tout, xout] = rk4(@sc_dynamics, b_omegaBN, [0;0;0], 1, 500, 0);
b_I = [10 0 0; 0 5 0; 0 0 7.5]; % [kg-m^2]
disp("H in body frame at 500s:")
b_H = b_I * xout(:, end);
disp(b_H);
% Provide Rotational Kinetic Energy T = 1/2 b_omegaBN' * I .* b_omegaBN at 500 seconds
disp("KE at 500s:")
T = 1/2* xout(:,end).' * b_I * xout(:,end); % first row of matrix
disp(T);
% Provide MRP attitude sigmaBN (500s)
x0 = [sigmaBN; b_omegaBN*pi/180];
[tout, xout] = rk4(@sc_dynamics_full, x0, [0;0;0], 1, 500, 1);
disp("MRP Attitude at 500s:")
disp(xout(1:3, end))
plot(tout, xout)

% Provide angular momentum vector n_H(500s) in inertial frame components
% Get DCM [BN] from the MRP
BN = mrp2dcm(xout(1:3,end));
disp("Angular Momentum at 500s:")
n_H = BN' * b_H;
disp(n_H);


% If you apply a fixed control torque b_u = (0.01, -0.01, 0.02) Nm provide
% the attitude sigmaBN (t = 100s)
[tout, xout] = rk4(@sc_dynamics_full, x0, [0.01;-0.01;0.02], 1, 100, 1);
disp("MRP Attitude with Fixed Control Torque at 100s:")
disp("u = [0.01; -0.01; 0.02]")
disp(xout(1:3,end))
plot(tout,xout)

function [t, yout] = rk4(fcn, x0, u_control, dt, tend, mrp_use)
% Runge-Kutta 4 Integrator
% Function should be setup as y_dot = f(x, u, t)
% Hold control vector u piece-wise constant over the RK4 integration to the
% next time step
% Can update the control u at ever time step in advance to the RK4
% integration step

t = 0:dt:tend;
yout = zeros(length(x0), length(t));
yout(:, 1) = x0;
for i = 1:(length(t)-1)
    % Update control vector
    u = u_control;
    
    
    k1 = fcn(yout(:, i),u, i);
    k2 = fcn(yout(:, i) + k1 * dt/2, u, i + dt/2);
    k3 = fcn(yout(:, i) + k2 * dt/2, u, i + dt/2);
    k4 = fcn(yout(:, i) + k3 * dt,   u, i + dt);
    
    % Calculate y(i + dt)
    yout(:, i+1) = yout(:, i) + (k1 + 2*k2 + 2*k3 + k4)*dt/6;
    
    % Check norm of MRP if flag set
    if (mrp_use)
       % If Using MRPs check norm and see if we're approaching singularity
       norm_o_b_n = norm(yout(1:3, i+1));
       if (norm_o_b_n > 1)
           yout(1:3, i+1) = -(yout(1:3,i+1) ./(norm_o_b_n^2));
           
       end
    end
end
end


function x_dot = sc_dynamics(x, u, t)
% x = [w_b_n(1) w_b_n(2) w_b_n(3)]
% Calculate x_dot at time t
b_I = [10 0 0; 0 5 0; 0 0 7.5]; % [kg-m^2]


b_w_b_n = x;

% Skew Matrix b_w_b_n_skew
b_w_b_n_skew = [0          -b_w_b_n(3) b_w_b_n(2); 
                b_w_b_n(3)  0         -b_w_b_n(1); 
               -b_w_b_n(2)  b_w_b_n(1) 0];


b_w_dot_b_n = inv(b_I)* (-b_w_b_n_skew * b_I * b_w_b_n + u);

x_dot = b_w_dot_b_n;



end

function x_dot = sc_dynamics_full(x, u, t)
% x = [o_b_n b_w_b_n]
% Calculate x_dot at time t
b_I = [10 0 0; 0 5 0; 0 0 7.5]; % [kg-m^2]

% Get States
o_b_n = x(1:3); 
b_w_b_n = x(4:6);

% Angular acceleration differential equation
b_w_b_n_skew = [0          -b_w_b_n(3) b_w_b_n(2); 
                b_w_b_n(3)  0         -b_w_b_n(1); 
               -b_w_b_n(2)  b_w_b_n(1) 0];


b_w_dot_b_n = b_I \ ((-b_w_b_n_skew * (b_I * b_w_b_n)) + u);


% MRP Differential Equation
skew_o_b_n = [0        -o_b_n(3) o_b_n(2); 
              o_b_n(3)  0       -o_b_n(1); 
             -o_b_n(2)  o_b_n(1) 0];


o_dot_b_n = 1/4 * ( (1-(o_b_n.'*o_b_n)) * eye(3) + 2 * skew_o_b_n + 2*(o_b_n*o_b_n.')) * b_w_b_n;

% Output x_dot
x_dot = [o_dot_b_n; b_w_dot_b_n];

end











% To get the tracking error sigma_BR:
% 1. convert sigma_BN to euler parameters, 
% 2. calculate DCMs of sigmaBN and RN, 
% 3. multiply the two DCMs (dot(DCM for sigma_BN, DCM for sigma_RN.T)), 
% 4. convert the result back to MRPs

% Attitude errors for all frames wrt spacecraft body frame
function [BR_error_attitute_in_b, BR_error_omega_in_b] = Errors(t, sigmaBN, b_omegaBN, dcmRN, n_omegaRN)

%   sigmaBN: MRP spacecraft’s initial attitude 
%   b_omegaBN: initial body angular velocity
%   N_omegaRN: calculated angular velocity Rx-frame wrt N
%   dcmRN: Calculated DCM of Rx-frame wrt N
%   t: time in seconds 

    dcmBN = mrp2dcm(sigmaBN);
    dcmRB = dcmRN*dcmBN.';
    BR_error_attitute_in_b = dcm2mrp(dcmRB.').';
    
    n_omegaBR = dcmBN.'*b_omegaBN - n_omegaRN;    
    BR_error_omega_in_b = dcmBN*n_omegaBR;

end

% Conversion MRP to DCM
function M = mrp2dcm(Q)
    Qn = norm(Q);
    Qt = tilde(Q);
    M = eye(3) + (8*Qt^2 - 4*(1-Qn^2)*Qt)/(1 + Qn^2)^2;
end

% Conversion DCM to MRP
function s = dcm2mrp(Q)
    b = dcm2quat(Q); % Always convert DCMs to Euler Parameters and then convert that to MRPs to avoid singularity
    s(1) = b(2)/(1+b(1));
    s(2) = b(3)/(1+b(1));
    s(3) = b(4)/(1+b(1));

end

% Tilde Matrix
function M_tilde = tilde(q)
    M_tilde = [0 -q(3) q(2); q(3) 0 -q(1); -q(2) q(1) 0];
end


% N-frame to Rc-frame - Comm mode (with GMO sat)
function angularvelocityRcN = omegaRcN(t)
    RcN_1 = RcN_t(t);
    RcN_2 = RcN_t(t+0.01); % finite difference to be considered
    omega = RcN_1.'*(RcN_2-RcN_1)/0.01; % DCM multiplied
    angularvelocityRcN = [-omega(3,2); -omega(1,3); -omega(2,1)];
end

function dcmRcN = RcN_t(t)
    [r_GMO] = EAtoDCMtoORBIT(20424.2,0,0,250,t);
    [r_LMO] = EAtoDCMtoORBIT(3.7962e+03,20,30,60,t);
    
    delta_r = r_GMO - r_LMO;
    delta_r = delta_r/norm(delta_r);
    r1 = -(delta_r);
    r2 = cross(delta_r,[0;0;1])/norm(cross(delta_r,[0;0;1]));
    r3 = cross(r1,r2);
    NRc = [r1,r2,r3];
    dcmRcN = NRc';
end



function [Rn] = EAtoDCMtoORBIT(radius, gamma, i, theta_t0, t)

%   radius: radius of the circular orbit
%   gamma: right ascension angle
%   i: inclination angle
%   theta: true latitude angle
%   t: time in seconds 

    mu = 42828.3;   % km^3/s^2 
    theta_dot = rad2deg(sqrt(mu/radius^3)); % rad/sec
    gamma = deg2rad(gamma);
    i = deg2rad(i);
    theta_at_t = theta_t0 + theta_dot*t;
    theta = deg2rad(theta_at_t);

% argument of perigee
    w = 0; % for circular orbit
    x = radius*cos(theta); % e = 0
    y = radius*sin(theta);

    % r_n = DCM * [x,y,0].'
    cos_g = cos(gamma);
    sin_g = sin(gamma);
    cos_i = cos(i);
    sin_i = sin(i);
    cos_w = cos(w);
    sin_w = sin(w);

% intertial position
    Rn = [ (cos_g*cos_w-sin_g*sin_w*cos_i)*x + (-cos_g*sin_w-sin_g*cos_w*cos_i)*y;
        (sin_g*cos_w+cos_g*sin_w*cos_i)*x + (-sin_g*sin_w+cos_g*cos_w*cos_i)*y;
        (sin_w*sin_i)*x + (cos_w*sin_i)*y;
    ];
end


% intertial velocity 
%x_dot = -sqrt(mu/radius)*sin(theta);
% y_dot = sqrt(mu/radius)*cos(theta);

% RdotN = [ (cos_g*cos_w-sin_g*sin_w*cos_i)*x_dot + (-cos_g*sin_w-sin_g*cos_w*cos_i)*y_dot;
%     (sin_g*cos_w+cos_g*sin_w*cos_i)*x_dot + (-sin_g*sin_w+cos_g*cos_w*cos_i)*y_dot;
%     (sin_w*sin_i)*x_dot + (cos_w*sin_i)*y_dot;
% ];
