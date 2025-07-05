volScalarField& spreadMicrolayerHeatSource = spreadMicrolayerHeatSourceSolid[i];

volScalarField& sharpMicrolayerHeatSource = sharpMicrolayerHeatSourceSolid[i];

// PETSC objects (for now, use the same context as fluid)

MicrolayerContext& microlayerContext = microlayerContextFluid[i];

IOdictionary& microlayerDict = microlayerDictFluid[i];