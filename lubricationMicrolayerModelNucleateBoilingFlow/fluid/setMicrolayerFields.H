// Bubble Objects

scalar& bubbleRadius1 = bubbleRadius1Fluid[i];

scalar& bubbleRadius2 = bubbleRadius2Fluid[i];

scalar& dt12 = dt12Fluid[i];

// Microlayer Objects

// FOAM objects

volScalarField& microlayerMassSource = microlayerMassSourceFluid[i];

IOdictionary& microlayerDict = microlayerDictFluid[i];

const scalar cSlope =  cSlopeFluid[i];

Info << "cSlope = " << cSlope << nl;

// Petsc Objects 

// Not initialized in createMicrolayerFields yet

//Vec& microlayerGrid = microlayerGridFluid[i];

//Vec& microlayerThickness = microlayerThicknessFluid[i];

//Info << "Successfully set PETSC fields" << nl;
