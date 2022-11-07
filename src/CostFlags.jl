
abstract type CostFlag end

struct MinimumFuel <: CostFlag end
struct MinimumFuelMayer <: CostFlag end

struct MinimumEnergyToFuel <: CostFlag end

struct MinimumTime <: CostFlag end