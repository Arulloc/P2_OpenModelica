package Sattrak
  model Satellite
    constant Real pi = Modelica.Constants.pi;
    constant Real d2r = Modelica.Constants.D2R "Degrees to radians pi/180";
    parameter Real ecc "Eccentricity";
    parameter Real M0 "Mean anomaly at Epoch (deg)";
    parameter Real N0 "Mean motion at Epoch (rev/d)";
    parameter Real Ndot2 "1st der Mean Motion /2 rev/d^2";
    parameter Real Nddot6 "2nd der Mean Motion /6 rev/d^3";
    parameter Real tstart "Simulation start time, seconds since Epoch (s)";
    parameter Real incl "inclination angle (deg)";
    parameter Real argper0 "argument of perigee at start time (deg)";
    parameter Real RAAN0 "RAAN at start time (deg)";
    Real J2 = 1.081874e-3 "Earth's Gravity Field";
    Real Re = 6378.135 "Earth's Radius (km)";
    Real M "Mean Anomaly (deg)";
    Real N "Mean Motion (rev/d)";
    Real E "Eccentric Anomaly (deg)";
    Real theta "True Anomaly (deg)";
    Real a "Semi-major axis (km)";
    Real r "satellite radial distance (km)";
    Real p_sat_pf[3] "Position, Perifocal coords";
    Real v_sat_pf[3] "Velocity Perifocal coords";
    Real RAAN "RAAN (deg)";
    Real argper "argument of perigee (deg)";
    Real angles[3] "3 angles between perifocal and ECI (rad)";
  initial equation
  N = N0 + 2*Ndot2*tstart/86400 +3*Nddot6*tstart^2/86400^2;
  M = M0 + (N*360.)/86400.*tstart + Ndot2*tstart^2*360/86400^2 + Nddot6*tstart^3*360/86400^3;
  RAAN = RAAN0 + (-3*J2*Re^2*cos(incl))/(2*a^2*(1-ecc^2)^2)*(N*360.)/86400.*tstart;
  argper = argper0 + 3*J2*Re^2*(5*(cos(incl))^2-1)/(2*a^2*(1-ecc^2)^2)*(N*360.)/86400.*tstart;
  
  equation
    M*d2r = E*d2r - ecc*sin(E*d2r);
    tan(theta*d2r/2.) = sqrt((1 + ecc)/(1 - ecc))*tan(E*d2r/2.);
    a^3*(N*2*pi/86400)^2 = 398600.4;
    der(M) = (N*360.)/86400. + 2.*Ndot2*time*360/86400^2 + 3.*Nddot6*time^2*360/86400^3;
    der(N) = 2*Ndot2/86400 + 6*Nddot6*time/86400^2;
    r = a*(1 - ecc^2)/(1 + ecc*cos(theta*d2r));
    der(RAAN) = ((-3*J2*Re^2*cos(incl))/(2*a^2*(1-ecc^2)^2))*(N*360.)/86400.;
    der(argper) = ((3*J2*Re^2*(5*(cos(incl))^2-1))/(2*a^2*(1-ecc^2)^2))*(N*360.)/86400.;
    p_sat_pf[1] = r*cos(theta*d2r);
    p_sat_pf[2] = r*sin(theta*d2r);
    p_sat_pf[3] = 0;
    der(p_sat_pf[1]) = v_sat_pf[1];
    der(p_sat_pf[2]) = v_sat_pf[2];
    der(p_sat_pf[3]) = v_sat_pf[3];
    angles[1] = -argper*d2r;
    angles[2] = -incl*d2r;
    angles[3] = -RAAN*d2r;
  end Satellite;


model SatTest
  Sattrak.Satellite GPS(tstart=64800., M0=308.6693 , N0=2.00564286187898 , ecc=0.0066028 , Ndot2=-.00000001, Nddot6=0., incl=55.5495, RAAN0=145.4799, argper0=51.9711);
  
  Real r "Satellite radius(km)";
  Real theta "Satellite true anomaly (deg)";
  Real E "eccentric Anomaly (deg)";
  Real M "Mean Anomaly (deg)";
  Real p_sat_ECI[3] "position of sat in ECI (km)";
  Real v_sat_ECI[3] "velocity of sat in ECI (km/s)";
  Real p_sat_ECF[3] "position of sat in ECF (km)";
  Real v_sat_ECF[3] "velocity of sat in ECF (km/s)";
equation
  E = mod(GPS.E, 360.);
  M = mod(GPS.M, 360.);
  r = mod(GPS.r, 360.);
  theta = mod(GPS.theta, 360.);
  (p_sat_ECI,v_sat_ECI) =sat_ECI(ang=GPS.angles,p_pf=GPS.p_sat_pf,v_pf=GPS.v_sat_pf);
  (p_sat_ECF,v_sat_ECF) = sat_ECF(ang=fill(Sattrak.theta_d.GMST, 3),p_ECI=sat_ECI.p_ECI,v_ECI=sat_ECI.v_ECI);
  

  
  annotation(
  Documentation(info = "GPS BIIR-2  (PRN 13)    
1 24876U 97035A   23064.61816216 -.00000001  00000+0  00000+0 0  9997
2 24876  55.5495 145.4799 0066028  51.9711 308.6693  2.00564286187898"), Icon);

end SatTest;

  function sat_ECI
    import Modelica.Mechanics.MultiBody.Frames.TransformationMatrices.axesRotations;
    import Modelica.Mechanics.MultiBody.Frames.TransformationMatrices.resolve2;
    
    input Real ang[3] "angles from perifocal to ECI(rad)";
    input Real p_pf[3] "Posn vector in Perifocal coords (km)";
    input Real v_pf[3] "Velocity vector in Perifocal coords (km/s)";
    
    output Real p_ECI[3] "Posn vector in ECI coord (km)";
    output Real v_ECI[3] "Velocity vector in ECI coord (km/s)";
  
  protected
    Real TM[3,3] = axesRotations(sequence = {3,1,3},angles = ang);
   
  algorithm
    
    p_ECI :=resolve2(TM,p_pf);
    v_ECI :=resolve2(TM,v_pf);
  end sat_ECI;

  model GndStn
    import Modelica.Mechanics.MultiBody.Frames.TransformationMatrices.axesRotations;
    
    parameter Real stn_long "Station longitude (degE)"; //from python
    parameter Real stn_lat "Station latitude (degN)"; //from python
    parameter Real stn_elev "Station elevation (m)"; //from python
    parameter Real f "Earth reference ellipsoid";
    Real p_stn_ECF[3] "Station coordinates in ECF (km)";
    Real TM[3,3] "Transform matrix from ECF to topo";
    
    Real Re = 6378.137 "Earth equitorial radius (km)";
  
  equation
    f = 1/298.25223563 "Earth reference ellipsoid flattening";
    e = sqrt(2*f-f^2) "ellipsoidal eccentricity";
    Nlat = Re / sqrt(1-e^2*sin(lat)^2) "Earth ellipsoidal radius of curvature of the meridian";
    
    u[3] = {-sin(stn_long/180*Modelica.Constants.pi), cos(stn_long/180*Modelica.Constants.pi), 0};
    v[3] = {-sin(stn_lat/180*Modelica.Constants.pi)*cos(stn_long/180*Modelica.Constants.pi),
                 -sin(stn_lat/180*Modelica.Constants.pi)*sin(stn_long/180*Modelica.Constants.pi),
                 cos(stn_lat/180*Modelica.Constants.pi)};
    w[3] = {cos(stn_lat/180*Modelica.Constants.pi)*cos(stn_long/180*Modelica.Constants.pi),
                 cos(stn_lat/180*Modelica.Constants.pi)*sin(stn_long/180*Modelica.Constants.pi),
                 sin(stn_lat/180*Modelica.Constants.pi)};
    TM[1,1:3] = u;
    TM[2,1:3] = v;
    TM[3,1:3] = w;
  end GndStn;

  function sat_ECF
    import Modelica.Mechanics.MultiBody.Frames.TransformationMatrices.axesRotations;
    import Modelica.Mechanics.MultiBody.Frames.TransformationMatrices.resolve2;
    
    
    input Real ang[3] "GMST angle (deg)";
    input Real p_ECI[3] "Posn vector in ECI coords (km)";
    input Real v_ECI[3] "Velocity vector in ECI coords (km/s)";
    
    output Real p_ECF[3] "Posn vector in ECF coord (km)";
    output Real v_ECF[3] "Velocity vector in ECF coord (km/s)";
    
  protected
    constant Real d2r = Modelica.Constants.D2R "Degrees to radians pi/180";
    Real TM[3,3] = {{cos(ang[1]*d2r),sin(ang[1]*d2r),0}, {-sin(ang[2]*d2r),cos(ang[2]*d2r),0},{0,0,1}} "rotation matrix from ECI to ECF";
  algorithm
  
    p_ECF :=resolve2(TM,p_ECI);
    v_ECF :=resolve2(TM,v_ECI);
  end sat_ECF;

  function range_ECF2topo
    import Modelica.Mechanics.MultiBody.Frames.TransformationMatrices.resolve2;
    
    input Real p_stn_ECF[3] "Position of station in ECF coords";
    input Real p_sat_ECF[3] "Position of satellite in ECF coords";
    input Real v_sat_ECF[3] "Relative Velocity of satellite in ECF coords";
    input Real TM[3, 3] "Transform matrix from ECF to topo";
    output Real p_sat_topo[3] "Position of satellite relative to station, topo coords (km)";
    output Real v_sat_topo[3] "Velocity of satellite relative to station, topo coords(km/s)";
  
  protected
    Real p_sat_topo_ECF[3] "store pos of sat relative to gs in topoECF coords";
    Real v_sat_topo_ECF[3] "store vel of sat relative to gs in topoECF coords";
    
  algorithm
    p_sat_topo_ECF := resolve2(TM,p_sat_ECF - p_stn_ECF) "computes pos of sat realtive to gs in topoECF";
    p_sat_topo := p_sat_topo_ECF;
    
    v_sat_topo_ECF := resolve2(TM, v_sat_ECF) "computes vel of sat relative to gs in topECF";
    v_sat_topo := v_sat_topo_ECF;
    
  end range_ECF2topo;

  function range_topo2look_angles
    import Modelica.Mechanics.MultiBody.Frames.TransformationMatrices.axesRotations;
    import Modelica.Mechanics.MultiBody.Frames.TransformationMatrices.resolve2;
    
    input Real p_sat_topo[3] "Position of satellite in topo coords (km)";
    input Real v_sat_topo[3] "Velocity of satellite in topo coords (km/s)";
    
    output Real Azimuth "Azimuth look angle (deg)";
    output Real Elevation "Elevation look angle (deg)";
    output Real Azrate "Azimuth look angle (deg/min)";
    output Real Elrate "Elevation look angle (deg/min)";
    output Real Rrate "Range rate to dish (km/s)";
  
  protected
    Real Re = 6378.137 "Earth equitorial radius (km)";
    Real f = 1/298.25223563 "Earth reference ellipsoid flattening";
    Real e = sqrt(2*f-f^2) "ellipsoidal eccentricity";
    
    Real Nlat = Re / sqrt(1-e^2*sin(lat)^2) "Earth ellipsoidal radius of curvature of the meridian";
    Real T[3,3] = { (Nlat + h) * cos(lat) * cos(long), (Nlat + h) * cos(lat) * sin(long), (Nlat * (1 - e^2) + h) * sin(lat) } "ECF Cartesian coords of tracking station"; 
    Real TM[3,3] = {p_sat_topo[1] + T[1], p_sat_topo[2] + T[2], p_sat_topo[3] +T[3]} "ECF to topo coordinates";
    
    Real R = sqrt(pow(TM[1], 2) + pow(TM[2], 2) + pow(TM[3], 2)) "calculating range";
    Real Rrate = (TM[1] * v_sat_topo[1] + TM[2]* v_sat_topo[2] + TM[3] * v_sat_topo[3])/R "range rate to dish";
    
    Real Az = atan2(TM[1], TM[2])*180/Modelica.Constants.pi "Azimuth look angle (deg)";
    Real Azrate = (-TM[1] * v_sat_topo[2] + TM[2] * v_sat_topo[1]) / (R^2 * cos(e)) * 60 "Azimuth rate (deg/min)";
    
    Real El = atan2(TM[3], sqrt(TM[1]^2 + TM[2]^2)) * 180 / Modelica.Constants.pi "Elevation look angle (deg)";
    Real Elrate = (TM[3] / R - Rrate * tan(e)) * 60 "Elevation look angle (deg/min)";
  algorithm
    Azimuth :=Az;
    Elevation := El;
    Azrate := Azrate;
    Elrate := Elrate;
    Rrate :=Rrate;
  
  end range_topo2look_angles;

  function theta_d
  
    input Real days "Number of days from J2000 to start of day in question";
    input Real hours "hours from midnight of the day in question to time in question";
    output Real GMST "GMST angle (deg)";
  
    Real Du = hours / 24.-0.5 "calculating Du time preceeding midnight in 0.5 intervals";
    Real Tu = Du /36525. "Tu expressed in Julian Centuries (days)";
    Real GMST_00h = 24110.5484 + 8640184.*Tu + 0.093104*Tu^2 - 6.2e-6 *Tu^3 "GMST at start of observations (s)";
    Real GMST_00h_mod = mod(GMST_00h,86400) "modulo of 86400 to remove integer multiples of days";
    Real theta_mid = 360*GMST_00h_mod/86400 "converting to radians";
    Real rTu = 1.002737909350795 + 5.9006e-11*Tu - 5.9e-15*Tu^2 "ratio of sideral seconds to mean solar seconds";
    
  algorithm
    GMST := theta_mid + 360*rTu* (hours/24.-0.5-days)*360/86400 "Greenwich Sidereal Time";
  end theta_d;
end Sattrak;
