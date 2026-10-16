function Q_ref = compute_reflected_light_flux(bond_albedo, Q_incident)
%COMPUTE_REFLECTED_LIGHT_FLUX Build the Lambert reflected-light recoil flux.
% Q_incident is the projected direct solar flux after shadow and eclipse
% [W/m^2]. Absorbed-light incident momentum is deliberately excluded.
arguments
    bond_albedo (1,1) double {mustBeFinite, mustBeNonnegative}
    Q_incident (:,:) double {mustBeFinite, mustBeNonnegative}
end

assert(bond_albedo <= 1, ...
    'cobird:thermal:compute_reflected_light_flux:badAlbedo', ...
    'Bond albedo must not exceed one.');

Q_ref = bond_albedo .* Q_incident;
end
