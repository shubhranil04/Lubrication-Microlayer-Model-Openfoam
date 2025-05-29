// Bubble Objects

PtrList<scalar> bubbleRadius1Fluid(fluidRegions.size());

PtrList<scalar> bubbleRadius2Fluid(fluidRegions.size());

PtrList<scalar> dt12Fluid(fluidRegions.size());

// Microlayer Objects

// FOAM objects

PtrList<volScalarField> microlayerMassSourceFluid(fluidRegions.size());

PtrList<IOdictionary> microlayerDictFluid(fluidRegions.size());

PtrList<scalar> cSlopeFluid(fluidRegions.size());

// PETSC objects
PtrList<Vec> microlayerGridFluid(fluidRegions.size());

PtrList<Vec> microlayerThicknessFluid(fluidRegions.size());

// Populate microlayer field pointer lists
forAll(fluidRegions, i)
{
    Info << "*** Initializing microlayer mass source for region "
         << fluidRegions[i].name() << nl << endl;

    microlayerMassSourceFluid.set(
        i,
        new volScalarField(
            IOobject(
                "microlayerMassSource",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE),
            fluidRegions[i],
            dimensionedScalar("microlayerMassSource", dimMass / (dimLength * dimLength * dimLength * dimTime), 0)));

    Info << "*** Reading microlayer properties for region "
         << fluidRegions[i].name() << nl << endl;

    microlayerDictFluid.set(
        i,
        new IOdictionary(
            IOobject(
                "microlayerProperties",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::AUTO_WRITE ))

    );

    cSlopeFluid.set(
        i,
        new scalar(microlayerDictFluid[i].lookupOrDefault<scalar>("cSlope", 0)));

    bubbleRadius1Fluid.set(
        i,
        new scalar(0.0));

    bubbleRadius2Fluid.set(
        i,
        new scalar(0.0));

    dt12Fluid.set(
        i,
        new scalar(0.0));
}
