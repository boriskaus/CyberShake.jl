"""
    citation(io::IO = stdout)

Print how to cite the CyberShake software that this package runs: the reference paper, a suggested
sentence for the body of a paper, the acknowledgement requested by the Southern California Earthquake
Center (SCEC), and the repository.
"""
function citation(io::IO = stdout)
    print(io, """
    CyberShake.jl is only an interface; the computations are done by the CyberShake software of the
    Southern California Earthquake Center (SCEC): https://github.com/SCECcode/cybershake-core (BSD 3-Clause).

    Please cite:
      Graves, R., Jordan, T.H., Callaghan, S. et al. CyberShake: A Physics-Based Seismic Hazard Model for
      Southern California. Pure Appl. Geophys. 168, 367-381 (2011). https://doi.org/10.1007/s00024-010-0161-6
      (SCEC Contribution 1354)

    Suggested text:
      The research described in this article used CyberShake software (Graves et al., 2011) published under
      the BSD-3 license.

    Acknowledgement:
      We would like to acknowledge use of the CyberShake software provided by the Southern California
      Earthquake Center (http://scec.org) which is funded by NSF Cooperative Agreement EAR-1600087 and
      USGS Cooperative Agreement G17AC00047.
    """)
    return nothing
end
