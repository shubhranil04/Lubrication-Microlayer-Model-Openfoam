
{
    // Update microlayerContext

    Info << "Updating microlayer model parameters ..." << nl;

    microlayerContext.dt = runTime.deltaT().value();

    microlayerContext.ub = max(microlayerDict.lookupOrDefault<scalar>("bubbleExpansionVelocity", 0), 0);

    microlayerContext.xCl += microlayerContext.ucl * microlayerContext.dt;

    microlayerContext.l = microlayerDict.lookupOrDefault<scalar>("bubbleRadius", 0);

    if (microlayerContext.l < microlayerContext.dsReq)
    {
        microlayerContext.l = 0;
        microlayerContext.isMl = PETSC_FALSE;
    }

    microlayerContext.del0 = microlayerDict.lookupOrDefault<Foam::scalar>("initialMicrolayerThickness", 1e-6);

    microlayerContext.rc = microlayerDict.lookupOrDefault<Foam::scalar>("bubbleFootRadiusOfCurvature", 1e-6);

    microlayerContext.delM = 0.3 * microlayerContext.rc;

    // Update mesh

    if (microlayerContext.isMl)
    {

        Info << "Remeshing ..." << nl;

        microlayerContext.n = max(int(microlayerContext.l / microlayerContext.dsReq) + 1, 11); // Have at least 11 nodes

        microlayerContext.ds = microlayerContext.l / (microlayerContext.n - 1);

        // Declare vectors to store old grid data and new grid data

        // Declare vectors to store old grid data and new grid data
        Vec sOld;
        Vec sNew;
        Vec xNew;
        Vec hOld;
        Vec hNew;

        // Create vectors
        PetscCall(VecDuplicate(microlayerContext.s, &sOld));
        PetscCall(VecDuplicate(microlayerContext.h, &hOld));

        PetscCall(VecCreate(PETSC_COMM_WORLD, &sNew));
        PetscCall(VecSetSizes(sNew, PETSC_DECIDE, microlayerContext.n));
        PetscCall(VecSetFromOptions(sNew));

        PetscCall(VecDuplicate(sNew, &xNew));
        PetscCall(VecDuplicate(sNew, &hNew));

        // Copy data into vectors
        PetscCall(VecCopy(microlayerContext.s, sOld));
        PetscCall(VecCopy(microlayerContext.h, hOld));

        // Get pointers to PETSc vectors
        const PetscScalar *hOldp, *sOldp;
        PetscInt nOld;
        PetscScalar *hNewp, *sNewp, *xNewp;

        PetscCall(VecGetArrayRead(sOld, &sOldp));
        PetscCall(VecGetArrayRead(hOld, &hOldp));
        PetscCall(VecGetSize(sOld, &nOld));
        PetscCall(VecGetArray(sNew, &sNewp));
        PetscCall(VecGetArray(hNew, &hNewp));
        PetscCall(VecGetArray(xNew, &xNewp));

        // Update grid points
        for (PetscInt i = 0; i < microlayerContext.n; i++)
        {
            sNewp[i] = i * microlayerContext.ds;
            xNewp[i] = microlayerContext.xCl + i * microlayerContext.ds;
        }

        // Update microlayer thickness

        Info << "Interpolating microlayer profile onto new mesh ..." << nl;

        for (PetscInt j = 0; j < microlayerContext.n; j++)
        {
            // for new nodes in between old grid end points
            if ((sNewp[j] > sOldp[0]) && (sNewp[j] < sOldp[nOld - 1]))
            {
                // Find k such that sOldp[k] <= s < sOldp[k+1]
                PetscInt k = 0;
                while ((k < nOld - 2) && (sOldp[k + 1] < sNewp[j]))
                {
                    ++k;
                }

                // Linear interpolation
                hNewp[j] = ((sOldp[k + 1] - sNewp[j]) * hOldp[k] + (sNewp[j] - sOldp[k]) * hOldp[k + 1]) / (sOldp[k + 1] - sOldp[k]);
            }

            // Boundary conditions
            if (sNewp[j] <= sOldp[0])
            {
                hNewp[j] = 0;
            }

            if (sNewp[j] >= sOldp[nOld - 1])
            {
                hNewp[j] = 3 * hNewp[j - 1] - 3 * hNewp[j - 2] + hNewp[j - 3];
            }
        }

        // Restore vectors
        PetscCall(VecRestoreArrayRead(hOld, &hOldp));
        PetscCall(VecRestoreArrayRead(sOld, &sOldp));
        PetscCall(VecRestoreArray(sNew, &sNewp));
        PetscCall(VecRestoreArray(xNew, &xNewp));
        PetscCall(VecRestoreArray(hNew, &hNewp));

        // Store new mesh date in the struct
        PetscCall(VecDuplicate(sNew, &microlayerContext.s));
        PetscCall(VecCopy(sNew, microlayerContext.s));

        PetscCall(VecDuplicate(xNew, &microlayerContext.x));
        PetscCall(VecCopy(xNew, microlayerContext.x));

        PetscCall(VecDuplicate(hNew, &microlayerContext.h));
        PetscCall(VecCopy(hNew, microlayerContext.h));

        // Destroy vectors
        PetscCall(VecDestroy(&sNew));
        PetscCall(VecDestroy(&xNew));
        PetscCall(VecDestroy(&hNew));
        PetscCall(VecDestroy(&sOld));
        PetscCall(VecDestroy(&hOld));

        // Initialize SNES

        SNES snes;
        Vec h, r;
        Mat J;
        PetscInt niter, maxit, maxf;

        PetscCall(VecDuplicate(microlayerContext.h, &h));
        PetscCall(VecCopy(microlayerContext.h, h));
        PetscCall(VecDuplicate(h, &r));

        PetscCall(MatCreate(PETSC_COMM_WORLD, &J));
        PetscCall(MatSetSizes(J, PETSC_DECIDE, PETSC_DECIDE, microlayerContext.n, microlayerContext.n));
        PetscCall(MatSetFromOptions(J));
        PetscCall(MatSeqAIJSetPreallocation(J, 7, NULL));

        PetscCall(SNESCreate(PETSC_COMM_WORLD, &snes));
        PetscCall(SNESSetFromOptions(snes));
        PetscCall(SNESGetTolerances(snes, &abstol, &rtol, &stol, &maxit, &maxf));
        PetscCall(SNESSetFunction(snes, r, FormFunction, &microlayerContext));
        PetscCall(SNESSetJacobian(snes, J, J, FormJacobian, &microlayerContext));

        Info << "Solving microlayer model ..." << nl;
        PetscCall(SNESSolve(snes, NULL, h));
        PetscCall(SNESGetIterationNumber(snes, &niter));
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "number of SNES iterations = %" PetscInt_FMT "\n\n", niter));

        PetscCall(VecCopy(h, microlayerContext.h));

        // Write files
        if (runTime.writeTime())
        {
            Info << "Writing microlayer data at time " << runTime.timeName() << nl;

            fileName writeDir = runTime.path() / runTime.timeName();

            mkDir(writeDir);

            fileName xml_file = writeDir / "xml.csv";

            PetscCall(PetscViewerASCIIOpen(PETSC_COMM_WORLD, xml_file.c_str(), &petscviewer));
            PetscCall(VecView(microlayerContext.x, petscviewer));
            PetscCall(PetscViewerDestroy(&petscviewer));

            fileName hml_file = writeDir / "hml.csv";

            PetscCall(PetscViewerASCIIOpen(PETSC_COMM_WORLD, hml_file.c_str(), &petscviewer));
            PetscCall(VecView(microlayerContext.h, petscviewer));
            PetscCall(PetscViewerDestroy(&petscviewer));
        }

        PetscCall(VecDestroy(&h));
        PetscCall(VecDestroy(&r));
        PetscCall(MatDestroy(&J));
        PetscCall(SNESDestroy(&snes));
    }
}