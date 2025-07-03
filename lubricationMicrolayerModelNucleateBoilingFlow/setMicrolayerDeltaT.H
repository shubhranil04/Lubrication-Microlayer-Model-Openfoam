forAll(fluidRegions, i)
{
    Info << "\nSetting time step using microlayer model for fluid region "
         << fluidRegions[i].name() << endl;
#include "setRegionFluidFields.H"

    const scalar currDeltaT = runTime.deltaT().value();

    const scalar microlayerDeltaT = max(1e-9, 0.1 * microlayerContext.ds / (mag(microlayerContext.ub) + 1e-10));

    if (microlayerDeltaT < currDeltaT)
    {
        Info << "Microlayer model time step = " << microlayerDeltaT << " s" << nl;
        runTime.setDeltaT(microlayerDeltaT);
    }
    else
    {
        Info << "Microlayer model time step is larger than current time step, using current time step." << nl;
    }
}