function Shock_wave=shock_wave(M,gamma)
 Shock_wave=sqrt((M.^2.*(gamma-1)+2)./(2.*gamma.*M.^2-(gamma-1))); 