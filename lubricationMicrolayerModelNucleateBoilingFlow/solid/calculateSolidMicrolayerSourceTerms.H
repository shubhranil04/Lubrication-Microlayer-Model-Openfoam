
// Get boundary grid points

{
    Info << "Calculating solid region microlayer source terms ... " << nl;

    const label patchIndex = mesh.boundaryMesh().findPatchID("solid_to_fluid");
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

    forAll(patch, faceI)
    {
        label meshFaceI = patch.start() + faceI;
	
        // Face area magnitude
        scalar faceAreaMag = mag(mesh.Sf()[meshFaceI]); //mesh.magSf()[meshFaceI]; // or: mesh.Sf()[meshFaceI].mag()
	
	scalar x1 = boundaryPoints[faceI];
        scalar x2 = Foam::min(boundaryPoints[faceI + 1], microlayerContext.rbCalc);

        // Calculate face area assuming wedge angle = 5 degrees 
        scalar calcFaceArea = microlayerContext.wedgeAngle * Foam::constant::mathematical::pi/360*(x2*x2-x1*x1);

        faceAreaMag = calcFaceArea;
        // Owner cell of the face
        label ownerCellI = mesh.faceOwner()[meshFaceI];

        // Volume of the owner cell
        scalar cellVolume = mesh.V()[ownerCellI];

        PetscScalar cellHeatFlux = computeMeanHeatFlux(x1, x2, &microlayerContext);
	
	if (boundaryPoints[faceI] < microlayerContext.rbCalc)
	{	
        	sharpMicrolayerHeatSource[ownerCellI] = - cellHeatFlux * faceAreaMag / cellVolume ;
	}
    }

    delete[] boundaryPoints;
	
    scalar maxHeatSource = -gMin(sharpMicrolayerHeatSource.internalField());

    Info << "Max Heat Source = " << maxHeatSource << nl;

    scalar spread_ = microlayerDict.lookupOrDefault<scalar>("spreadSolid", 3);

    dimensionedScalar DPsi
    (
        "DPsi",
        dimensionSet(0,2,0,0,0,0,0),
        spread_/sqr(gAverage(mesh.nonOrthDeltaCoeffs())
    ));

    //- Smearing of source term field
    fvScalarMatrix psiEqn
    (
        fvm::Sp(scalar(1),spreadMicrolayerHeatSource) - fvm::laplacian(DPsi,spreadMicrolayerHeatSource) == sharpMicrolayerHeatSource
    );

    psiEqn.solve();
	
    // smear microlayer mass source terms
    dimensionedScalar intPsi0 = fvc::domainIntegrate(sharpMicrolayerHeatSource);

    // smear microlayer mass source terms
    dimensionedScalar intPsi = fvc::domainIntegrate(spreadMicrolayerHeatSource) + dimensionedScalar("smallFix", spreadMicrolayerHeatSource.dimensions()*dimLength*dimLength*dimLength, SMALL);;

    dimensionedScalar Ns ("Ns", dimless, 1.0);

    Ns = intPsi0 / intPsi;

    forAll(mesh.C(),celli)
    {
	spreadMicrolayerHeatSource[celli] = Ns.value() * spreadMicrolayerHeatSource[celli];
    }
}
