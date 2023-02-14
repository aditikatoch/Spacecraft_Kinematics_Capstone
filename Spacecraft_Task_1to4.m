clear;
format short; 

R = 3396.19; %km
h = 400; %km
r_LMO = R + h; %km
r_GMO = 20424.2; %km
mu = 42828.3; %km^3 /s^2

% pv = r * ir = r * DCM * ir
% dpv = wbn X r = DCM * dtheta X r * ir

EAtoDCM = @computeEAtoDCM;

% LMO - Hill
EA_LMO = [deg2rad(20) deg2rad(30) deg2rad(60)].';
rate_LMO = sqrt(mu/(r_LMO.^3)); %0.000884797 rad/s
dEA_LMO = [0, 0, rate_LMO].';
b_LMO = [r_LMO, 0, 0].'; % nano-sat
n_LMO = EAtoDCM(EA_LMO).'* b_LMO; % LMO in intertial frame

% GMO
EA_GMO = [deg2rad(0) deg2rad(0) deg2rad(250)].';
rate_GMO = sqrt(mu/(r_GMO.^3)); %0.0000709003 rad/s
dEA_GMO = [0, 0, rate_GMO].';
b_GMO = [r_GMO, 0, 0].';
n_GMO = EAtoDCM(EA_GMO).'* b_GMO; % GMO in intertial frame

% POWER MODE - Sun pointing
EA_SUN = [deg2rad(180), deg2rad(0), deg2rad(90)].';
dEA_SUN = [0 0 0].';


% SCIENCE MODE - Nadir pointing (center of Mars)
EA_MARS_LMO = [deg2rad(180) deg2rad(0) deg2rad(0)].';
dEA_MARS = [0 0 rate_LMO].';



h = 0.001;
time = 1200;
prev_t = 0;

for ti = linspace(0, time, (time/h + 1))
    dt = ti - prev_t;
    prev_t = ti;
    EA_LMO_t = EA_LMO + dEA_LMO * dt;
    EA_GMO_t = EA_GMO + dEA_GMO * dt;
    
    n_LMO_t = n_LMO + dEA_LMO * dt;
    n_GMO_t = n_GMO + dEA_GMO * dt;
    % COMM MODE - GMO facing
    delta_r = n_GMO_t - n_LMO_t;
    Rc_1 = -delta_r/norm(delta_r);
    Rc_2 = cross(delta_r, [0 0 1].')/norm(cross(delta_r, [0 0 1].'));  
    Rc_3 = cross(Rc_1, Rc_2);
    RcN = [Rc_1 Rc_2 Rc_3].';
    
    % Internial Position and velocity
    Nr_LMO_t = EAtoDCM(EA_LMO_t).' * b_LMO;
    Nr_GMO_t = EAtoDCM(EA_GMO_t).' * b_GMO;

    Nrdot_LMO = cross(dEA_LMO, b_LMO);
    Nrdot_LMO = EAtoDCM(EA_LMO_t).'* Nrdot_LMO;
    Nrdot_GMO = cross(dEA_GMO, Nr_GMO_t);
  
    %% TASK 3
    if ti == 0
        disp("DCM [RsN]");
        disp(angle2dcm(deg2rad(180), deg2rad(0), deg2rad(90)));
        disp("Angular velocity RsN");
        disp(cross(EA_SUN, dEA_SUN));
    end

    %% TASK 2
    if ti == 300
        disp("DCM [HN]");
        disp(EAtoDCM(EA_LMO_t));
    end

    if ti == 0
       %% TASK 5
        disp("DCM [RcN]");
        disp(RcN);

       %% TASK 4 
       RnN_MARS = angle2dcm(deg2rad(180), deg2rad(0), deg2rad(180)) * (EAtoDCM(EA_LMO_t)); %DCMs are multiplied
       disp("DCM [RnN]")
       disp(RnN_MARS);
       omega_RnN_MARS = RnN_MARS.'*(-dEA_MARS);
       disp("Angular velocity RnN ")
       disp(omega_RnN_MARS); 

    end
    
    %% TASK 1
    if ti == 450
        disp("Position LMO");
        disp(Nr_LMO_t);
        disp("Velocity LMO");
        disp(Nrdot_LMO);

    end

    if ti == 1150
        disp("Position GMO");
        disp(Nr_GMO_t);
        disp("Velocity GMO");
        disp(Nrdot_GMO);

    end
   

        
    EA_LMO = EA_LMO_t;
    EA_GMO = EA_GMO_t;

end



function DCM = computeEAtoDCM(EA)
    s1 = sin(EA(1));
    c1 = cos(EA(1));
    s2 = sin(EA(2));
    c2 = cos(EA(2));
    s3 = sin(EA(3));
    c3 = cos(EA(3));
    
    DCM = [[c3*c1-s3*c2*s1, c3*s1+s3*c2*c1, s3*s2];
                    [-s3*c1-c3*c2*s1, -s3*s1+c3*c2*c1, c3*s2]; 
                    [s2*s1, -s2*c1, c2]];
end
