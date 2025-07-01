// FOAM objects

PtrList<volScalarField> microlayerHeatSourceFluid(fluidRegions.size());

// Populate microlayer field pointer lists
forAll(fluidRegions, i)
{
    Info << "*** Initializing microlayer heat source for region "
         << solidRegions[i].name() << nl << endl;

    microlayerHeatSourceFluid.set(
        i,
        new volScalarField(
            IOobject(
                "microlayerHeatSource",
                runTime.timeName(),
                solidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE),
            solidRegions[i],
            dimensionedScalar("microlayerHeatSource", dimPower / (dimLength * dimLength * dimLength), 0)));
}