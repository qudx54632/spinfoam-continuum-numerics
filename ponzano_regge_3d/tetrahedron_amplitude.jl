using WignerSymbols

"""
    tetrahedron_amplitude((a, b, c), (d, e, f))

Return the Ponzano-Regge amplitude of one tetrahedron,

                  { a  b  c }
    A_t(j_e) =    { d  e  f }.

The arguments are SU(2) spins `j`, not edge lengths `j + 1/2`.
Half-integer spins should be written exactly, for example `1//2`.
"""
function tetrahedron_amplitude(
    top::Tuple{Real,Real,Real},
    bottom::Tuple{Real,Real,Real},
)
    a, b, c = top
    d, e, f = bottom
    return wigner6j(a, b, c, d, e, f)
end
