// Get boundary grid points

{
    Info << "Calculating microlayer source terms ... " << nl;


    const label patchIndex = mesh.boundaryMesh().findPatchID("fluid_to_solid");
    const fvPatch &patch = mesh.boundary()[patchIndex];
    const PetscInt nBoundaryFaces = patch.size();
    PetscScalar *boundaryPoints = new PetscScalar[nBoundaryFaces + 1];

    // Loop over faces in the patch
    forAll(patch, faceI)
    {
        label meshFaceI = patch.start() + faceI;
        const face &f = mesh.faces()[meshFaceI];

        scalar xMin = 1e20;
        scalar xMax = -1e20;
        // Loop over points on this face
        forAll(f, vertI)
        {
            label pointI = f[vertI];
            const point &pt = mesh.points()[pointI];
            scalar xCoord = pt.x();

            xMin = Foam::min(xMin, xCoord);
            xMax = Foam::max(xMax, xCoord);
        }

        boundaryPoints[faceI] = xMin;
        boundaryPoints[faceI + 1] = xMax;
    }

    scalar spread_ = microlayerDict.lookupOrDefault<scalar>("spread", 3);

    dimensionedScalar DPsi
    (
        "DPsi",
        dimensionSet(0,2,0,0,0,0,0),
        spread_/sqr(gAverage(mesh.nonOrthDeltaCoeffs())
    ));

    forAll(patch, faceI)
    {
        label meshFaceI = patch.start() + faceI;

        // Face area magnitude
        scalar faceAreaMag = mesh.magSf()[meshFaceI]; // or: mesh.Sf()[meshFaceI].mag()

        // Owner cell of the face
        label ownerCellI = mesh.faceOwner()[meshFaceI];

        // Volume of the owner cell
        scalar cellVolume = mesh.V()[ownerCellI];

        PetscScalar cellHeatFlux = computeMeanHeatFlux(boundaryPoints[faceI], boundaryPoints[faceI + 1], &microlayerContext);

	if (alpha1[ownerCellI] < 1 - 1e-3)
        {
            sharpMicrolayerMassSource[ownerCellI] = cellHeatFlux * faceAreaMag / cellVolume / microlayerContext.hfg;
        }
    }

    delete[] boundaryPoints;

    // smear microlayer mass source terms
    dimensionedScalar intPsi0 = fvc::domainIntegrate(sharpMicrolayerMassSource);

    //- Smearing of source term field
    fvScalarMatrix psiEqn
    (
        fvm::Sp(scalar(1),spreadMicrolayerMassSource) - fvm::laplacian(DPsi,spreadMicrolayerMassSource) == sharpMicrolayerMassSource
    );

    psiEqn.solve();

    // Cut cells with cutoff < alpha1 < 1-cutoff
    // and rescale remaining source term field
    dimensionedScalar intPsiLiquid
    (
        "intPsiLiquid",
        dimensionSet(1,0,-1,0,0,0,0),
        0.0
    );
    dimensionedScalar intPsiVapor
    (
        "intPsiLiquid",
        dimensionSet(1,0,-1,0,0,0,0),
        0.0
    );

    forAll(mesh.C(),celli)
    {
        if (alpha1[celli] < 1e-3)
        {
            intPsiVapor.value() +=
                (1.0-alpha1[celli])*spreadMicrolayerMassSource[celli]*mesh.V()[celli];
        }
        else if (alpha1[celli] > 1.0-1e-3)
        {
            intPsiLiquid.value() += alpha1[celli]*spreadMicrolayerMassSource[celli]*mesh.V()[celli];
        }
    }

    //- Calculate Nl and Nv
    dimensionedScalar Nl ("Nl", dimless, 2.0);
    dimensionedScalar Nv ("Nv", dimless, 2.0);

    reduce(intPsiLiquid.value(),sumOp<scalar>());
    reduce(intPsiVapor.value(),sumOp<scalar>());

     if (intPsiLiquid.value() > 1e-99)
    {
        Nl = intPsi0/intPsiLiquid;
    }
    if (intPsiVapor.value() > 1e-99)
    {
        Nv = intPsi0/intPsiVapor;
    }

    //- Set source terms in cells with alpha1 < cutoff or alpha1 > 1-cutoff
    forAll(mesh.C(),celli)
    {
        if (alpha1[celli] < 1e-3)
        {
            spreadMicrolayerMassSource[celli] = Nv.value()*(1.0-alpha1[celli])*spreadMicrolayerMassSource[celli];
        }
        else if (alpha1[celli] > 1.0-1e-3)
        {
            spreadMicrolayerMassSource[celli] = 0;
        }
        else
        {
            spreadMicrolayerMassSource[celli] = 0.0;
        }
    }
}
