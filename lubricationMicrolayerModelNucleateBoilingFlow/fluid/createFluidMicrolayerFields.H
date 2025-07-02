// Bubble Objects

PtrList<scalar> bubbleRadius1Fluid(fluidRegions.size());

PtrList<scalar> bubbleRadius2Fluid(fluidRegions.size());

PtrList<scalar> dt12Fluid(fluidRegions.size());

// Microlayer Objects

// FOAM objects

PtrList<volScalarField> sharpMicrolayerMassSourceFluid(fluidRegions.size());

PtrList<volScalarField> spreadMicrolayerMassSourceFluid(fluidRegions.size());

PtrList<IOdictionary> microlayerDictFluid(fluidRegions.size());

// PETSC objects

PtrList<MicrolayerContext> microlayerContextFluid(fluidRegions.size());

// Populate microlayer field pointer lists
forAll(fluidRegions, i)
{
    Info << "*** Initializing microlayer mass source for region "
         << fluidRegions[i].name() << nl << endl;

    sharpMicrolayerMassSourceFluid.set(
        i,
        new volScalarField(
            IOobject(
                "sharpMicrolayerMassSource",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE),
            fluidRegions[i],
            dimensionedScalar("sharpMicrolayerMassSource", dimMass / (dimLength * dimLength * dimLength * dimTime), 0)));

    spreadMicrolayerMassSourceFluid.set(
        i,
        new volScalarField(
            IOobject(
                "spreadMicrolayerMassSource",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE),
            fluidRegions[i],
            dimensionedScalar("spreadMicrolayerMassSource", dimMass / (dimLength * dimLength * dimLength * dimTime), 0),
            "zeroGradient"));

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
                IOobject::AUTO_WRITE))

    );

    bubbleRadius1Fluid.set(
        i,
        new scalar(0.0));

    bubbleRadius2Fluid.set(
        i,
        new scalar(0.0));

    dt12Fluid.set(
        i,
        new scalar(0.0));

    Info << "*** Setting microlayer PETSC objects for region "
         << fluidRegions[i].name() << nl << endl;

    microlayerContextFluid.set(
        i,
        new MicrolayerContext(microlayerDictFluid[i]));
}
// ************************************************************************* //
