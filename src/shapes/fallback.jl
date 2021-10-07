has_symmetric_plane(scatterer::AbstractScatterer) = false

function volume_equivalent_radius(scatterer::AbstractScatterer)
    throw(ArgumentError("volume_equivalent_radius() has not been implemented for $(typeof(scatterer))."))
    return
end
