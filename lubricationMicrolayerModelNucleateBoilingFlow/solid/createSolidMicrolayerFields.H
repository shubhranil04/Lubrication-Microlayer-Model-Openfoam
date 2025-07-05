// FOAM objects

PtrList<volScalarField> sharpMicrolayerHeatSourceSolid(solidRegions.size());
PtrList<volScalarField> spreadMicrolayerHeatSourceSolid(solidRegions.size());

// Populate microlayer field pointer lists
forAll(solidRegions, i)
{
    Info << "*** Initializing microlayer heat source for region "
         << solidRegions[i].name() << nl << endl;

    sharpMicrolayerHeatSourceSolid.set(
        i,
        new volScalarField(
            IOobject(
                "sharpMicrolayerHeatSource",
                runTime.timeName(),
                solidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE),
            solidRegions[i],
            dimensionedScalar("sharpMicrolayerHeatSource", dimPower / (dimLength * dimLength * dimLength), 0)));

    spreadMicrolayerHeatSourceSolid.set(
        i,
        new volScalarField(
            IOobject(
                "spreadMicrolayerHeatSource",
                runTime.timeName(),
                solidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE),
            solidRegions[i],
            dimensionedScalar("spreadMicrolayerHeatSource", dimPower / (dimLength * dimLength * dimLength), 0),
            "zeroGradient"));
}