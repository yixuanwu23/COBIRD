function delta = compute_tidal_linearization_rate(torque_magnitude, inertia_z, density, G)
%COMPUTE_TIDAL_LINEARIZATION_RATE Jacobson near-synchronous transition rate.
% delta = (|Gamma|/I_z) * sqrt(6/(pi*G*rho)).

arguments
    torque_magnitude (1,1) double {mustBeNonnegative, mustBeFinite}
    inertia_z (1,1) double {mustBePositive, mustBeFinite}
    density (1,1) double {mustBePositive, mustBeFinite}
    G (1,1) double {mustBePositive, mustBeFinite}
end

delta = (torque_magnitude / inertia_z) * sqrt(6 / (pi * G * density));
end
