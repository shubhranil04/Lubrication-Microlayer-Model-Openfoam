// Bubble Objects

scalar& bubbleRadius1 = bubbleRadius1Fluid[i];

scalar& bubbleRadius2 = bubbleRadius2Fluid[i];

scalar& dt12 = dt12Fluid[i];

// Microlayer Objects

// FOAM objects

volScalarField& microlayerMassSource = microlayerMassSourceFluid[i];

IOdictionary& microlayerDict = microlayerDictFluid[i];

// PETSC objects

MicrolayerContext& microlayerContext = microlayerContextFluid[i];


// ************************************************************************* //
