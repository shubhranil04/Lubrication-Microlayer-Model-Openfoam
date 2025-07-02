volScalarField& microlayerHeatSource = microlayerHeatSourceSolid[i];

// PETSC objects (for now, use the same context as fluid)

MicrolayerContext& microlayerContext = microlayerContextFluid[i];