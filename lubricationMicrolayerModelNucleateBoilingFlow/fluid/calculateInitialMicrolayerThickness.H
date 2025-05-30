// Look up PLIC normals and positions.

const auto &meshPoints = mesh.points();
const auto &meshFaces = mesh.faces();

const auto &db = alpha1.db();

const auto normalsName = IOobject::groupName(
    "interfaceNormal",
    alpha1.internalField().group());
const auto centresName = IOobject::groupName(
    "interfaceCentre",
    alpha1.internalField().group());

bool hasNormals = db.found(normalsName);
bool hasCentres = db.found(centresName);

// Get patch fields for interface normals and centers
const label patchIndex = mesh.boundaryMesh().findPatchID("fluid_to_solid");
const fvPatch &patch = mesh.boundary()[patchIndex];

bool isBubble = false;

if (hasNormals && hasCentres)
{
    // Info << "Interpolating the PLIC" << nl;
    const volVectorField &interfaceCentre =
        db.lookupObject<volVectorField>(centresName);

    const volVectorField &interfaceNormal =
        db.lookupObject<volVectorField>(normalsName);

    const auto &pInterfaceCentres = interfaceCentre.boundaryField()[patchIndex];
    const auto &pInterfaceNormals = interfaceNormal.boundaryField()[patchIndex];

    const auto pInternalCentresTmp = pInterfaceCentres.patchInternalField();
    const auto pInternalNormalsTmp = pInterfaceNormals.patchInternalField();

    const auto &pInternalCentres = pInternalCentresTmp.cref();
    const auto &pInternalNormals = pInternalNormalsTmp.cref();

    scalar numContactLineCells = 0.0;
    scalar bubbleRadius = 0.0;

    forAll(pInternalCentres, faceI)
    {
        bool hasContactLine = false;

        const vector &cellInterfaceNormal = pInternalNormals[faceI];
        const vector &cellInterfaceCentre = pInternalCentres[faceI];

        const auto &thisFace = meshFaces[patch.start() + faceI];

        // Get face points.
        for (auto pointI = 0; pointI < (thisFace.size() - 1); ++pointI)
        {
            // Compute the signed distance of the first point.
            const point &firstFacePoint = meshPoints[thisFace[pointI]];
            const scalar firstDist = (firstFacePoint - cellInterfaceCentre) & cellInterfaceNormal;

            // Compute the signed distance of the second point.
            const point &secondFacePoint = meshPoints[thisFace[pointI + 1]];
            const scalar secondDist = (secondFacePoint - cellInterfaceCentre) & cellInterfaceNormal;

            if (firstDist * secondDist < 0)
            {
                hasContactLine = true;
                break;
            }
        }

        if (hasContactLine)
        {
            numContactLineCells += 1;
            bubbleRadius += pInternalCentres[faceI].x();
        }
    }

    if (numContactLineCells > 0)
    {
        bubbleRadius /= numContactLineCells;

        scalar dt01 = runTime.deltaTValue();

        scalar bubbleExpansionVelocity = (bubbleRadius - bubbleRadius1) / (dt01 + 1e-20);

        scalar bubbleExpansionAcceleration = 2 / ((dt01 + 1e-20) * (dt12 + 1e-20) * (dt01 + dt12 + 1e-20)) * (bubbleRadius * dt12 * dt12 - bubbleRadius1 * (dt01 + dt12) * (dt01 + dt12) + bubbleRadius2 * dt01 * dt01);

        Info << "" << nl << "Bubble radius = " << bubbleRadius << " m" << nl;
        Info << "" << nl << "Bubble expansion velocity = " << bubbleExpansionVelocity << " m" << nl;
        Info << "Bubble expansion acceleration = " << bubbleExpansionAcceleration << " m/s" << nl << "" << nl;

        microlayerDict.set("bubbleRadius", bubbleRadius);
        microlayerDict.set("bubbleExpansionVelocity", bubbleExpansionVelocity);
        microlayerDict.set("bubbleExpansionAcceleration", bubbleExpansionAcceleration);

        bubbleRadius2 = bubbleRadius1;
        bubbleRadius1 = bubbleRadius;

        dt12 = dt01;
    }
    else
    {
        Info << "No contact line cells found" << nl << "" << nl;
    }
}

else
{
    Info << "No normals and centres found" << nl;
}
