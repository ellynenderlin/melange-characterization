function [Ft] = calculate_iceberg_torque(H,ang,E,rho_i,rho_w,calve_mode)

% friction coefficients for terminus and fjord bottom
muT=0; muN=0;

% calculate gammas (iceberg length above water)
if contains(calve_mode,'float')
    gammas(3)= H*(1 - rho_i/rho_w - (E/2)*tand(ang));
else %grounded
    gammas(3)=H*(1-(rho_i/rho_w)*secd(ang));
end
gammas(2) = NaN;
gammas(1) = gammas(3)+E*H*tand(ang);

%calculate torques
taus(1) = (1/12)*rho_w*9.81*cosd(ang)*(H-gammas(1))^2*(H+2*gammas(1));
taus(2) = (1/12)*rho_w*9.81*sind(ang)*E^3*H^3;
taus(3) = -(1/12)*rho_w*9.81*cosd(ang)*(H-gammas(3))^2*(H+2*gammas(3));
tauB=sum(taus);

% force due to gravity
Fg=rho_i*9.81*H^2*E;

% force due to buoyancy
Fb=rho_w*9.81*E*H/2*(2*H-gammas(3)-gammas(1));

% normal force on lowest corner of the iceberg
Fn=Fg-Fb;

% torque from normal force on bottom of iceberg
tauN=Fn*H/2*(E*cosd(ang)-sind(ang));

%solve separate conditions for floating vs lightly grounded termini
if contains(calve_mode,'float')
    Ft = -tauB/(gammas(3)*cosd(ang)); %bottom-out rotation
    % Ft = -tauB/((H*(cosd(ang) - E*sind(ang))) - gammas(3)*cosd(ang)); %top-out rotation
else
    Ft=(-tauB-tauN-muN*Fn*H/2*(cosd(ang)+E*sind(ang)))/ ... 
        (H*(cosd(ang)-(rho_i/rho_w)-muT/2*(sind(ang)+E*cosd(ang))));
end


end
