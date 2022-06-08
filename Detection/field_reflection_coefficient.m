function r = field_reflection_coefficient(dielectric_constant1, dielectric_constant2)
    r = (dielectric_constant1 - dielectric_constant2)./(dielectric_constant1 + dielectric_constant2);
end