// Get boundary grid points

{
    Info << "Calculating microlayer source terms ... " << nl;

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
        scalar faceAreaMag = mesh.magSf()[meshFaceI]; // or: mesh.Sf()[meshFaceI].mag()

        // Owner cell of the face
        label ownerCellI = mesh.faceOwner()[meshFaceI];

        // Volume of the owner cell
        scalar cellVolume = mesh.V()[ownerCellI];

        PetscScalar cellHeatFlux = computeMeanHeatFlux(boundaryPoints[faceI], boundaryPoints[faceI + 1], &microlayerContext);

        Info << "Cell " << ownerCellI << ": Heat flux = " << cellHeatFlux << nl;

        microlayerHeatSource[ownerCellI] = - cellHeatFlux * faceAreaMag / cellVolume ;

    }

    delete[] boundaryPoints;
}
