# CR3BP Parameters
struct CRTBPParams 

    # Primary and Secondary Body Mass
    m1::Float64
    m2::Float64

    # Radius of Primary and Secondary Bodies
    R1::Float64
    R2::Float64

    # Total Mass
    mtot::Float64

    # Gravitational Parameters
    μ::Float64

    # Length Unit
    LU::Float64

    # Time Unit
    TU::Float64

    # Speed Unit
    VU::Float64

    # Mass Unit
    MU::Float64

    # Gravitational Constant
    G::Float64

    function CRTBPParams(m1,m2,R1,R2,MU,r12) 
        G       = 6.673e-20
        mtot    = m1 + m2
        μ       = m2 / (m1 + m2)
        LU      = r12 
        TU      = 1.0 / sqrt(G * (m1 + m2) / (LU^3))
        VU      = LU / TU 
        MU      = MU 
        new(m1,m2,R1,R2,mtot,μ,LU,TU,VU,MU,G)
    end
end