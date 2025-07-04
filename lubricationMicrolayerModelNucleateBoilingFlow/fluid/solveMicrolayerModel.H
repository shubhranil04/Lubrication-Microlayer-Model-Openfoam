
{
    // Update microlayerContext

    Info << "Updating microlayer model parameters ..." << nl;

    if (mag(microlayerDict.lookupOrDefault<scalar>("bubbleExpansionVelocity", microlayerContext.ub)) < microlayerContext.umax)
    {
        microlayerContext.ub = microlayerDict.lookupOrDefault<scalar>("bubbleExpansionVelocity", microlayerContext.ub);
    }

    //microlayerContext.del0 = microlayerDict.get<Foam::scalar>("initialMicrolayerThickness");

    //microlayerContext.rc = microlayerDict.get<Foam::scalar>("bubbleFootRadiusOfCurvature");

    //microlayerContext.delM = 10 * microlayerContext.del0; // min(0.3 * microlayerContext.rc, 20 * microlayerContext.del0);

    // microlayerContext.dt = runTime.deltaT().value();

    // microlayerContext.xCl += microlayerContext.ucl * microlayerContext.dt;

    // microlayerContext.l += (microlayerContext.ub - microlayerContext.ucl) * microlayerContext.dt; //= microlayerDict.lookupOrDefault<scalar>("bubbleRadius", 0);

    // Divide foam time step into microlayer time step

    scalar dtFoam = runTime.deltaTValue();
    scalar dtReq = min(microlayerContext.dtMax, 0.5 * microlayerContext.ds / (mag(microlayerContext.ub) + 1e-10));
    label n_time_steps = int(dtFoam / dtReq);

    if (n_time_steps < 1)
    {
        n_time_steps = 1;
    }

    microlayerContext.dt = dtFoam / n_time_steps;

    Info << "Divided foam time step into " << n_time_steps << " microlayer time steps." << nl;

    for (label nt = 0; nt < n_time_steps; ++nt)
    {   
        //const scalar B = microlayerDict.get<scalar>("growthConstant");

        //scalar Rb = microlayerContext.xCl + microlayerContext.l;

        //scalar Rbdot = 2 * B * B * microlayerContext.alphal / Rb;

        //microlayerContext.ub = 2 * B * B * microlayerContext.alphal / Rb * Foam::sin(microlayerContext.thetaApp);

        microlayerContext.xCl += microlayerContext.ucl * microlayerContext.dt;
        microlayerContext.l += (microlayerContext.ub - microlayerContext.ucl) * microlayerContext.dt;

        if (microlayerContext.l < microlayerContext.dsReq)
        {
            microlayerContext.l = 0;
            microlayerContext.isMl = PETSC_FALSE;
        }
        // Update mesh

        if (microlayerContext.isMl)
        {
            // Compute initial microlayer thickness

            const word initialMicrolayerThicknessModel = microlayerDict.lookupOrDefault<word>("initialMicrolayerThicknessModel", "linear");

            if (initialMicrolayerThicknessModel == "linear")
            {
                Info << "Using linear initial microlayer thickness model" << nl;

                const scalar cSlope = microlayerDict.lookupOrDefault<scalar>("cSlope", 4.46e-3); // Default value from Utaka et al. 2013

                Info << "cSlope = " << cSlope << nl;

                microlayerContext.del0 = cSlope * (microlayerContext.xCl + microlayerContext.l);

                Info << "Initial microlayer thickness = " << microlayerContext.del0 * 1e6 << " micrometers" << nl << "" << nl;

                // Calculate meniscus curvature based on the initial microlayer thickness
                microlayerContext.rc = microlayerContext.del0 / (1.34 * Foam::pow(microlayerContext.mul * microlayerContext.ub / microlayerContext.sigma, 2.0 / 3.0));

                microlayerContext.curv = 1.0 / microlayerContext.rc;

                microlayerContext.delM = 0.3 * microlayerContext.rc; //10 * microlayerContext.del0;
            }

            else if (initialMicrolayerThicknessModel == "landauLevichScriven")
            {
                Info << "Using Landau-Levich initial microlayer thickness model (Zhang et al. 2024) with Scriven solution" << nl;

                const scalar B = microlayerDict.get<scalar>("growthConstant");

                // scalar Rb = microlayerDict.lookupOrDefault<scalar>("bubbleRadius", 0) + 1e-20;

                scalar Rb = microlayerContext.xCl + microlayerContext.l;

                scalar Rbdot = 2 * B * B * microlayerContext.alphal / Rb;

                const scalar Rbddot = -4 * Foam::pow(B, 4) * Foam::pow(microlayerContext.alphal, 2) / Foam::pow(Rb, 3);

                // Coefficients

                const scalar a3 = microlayerContext.rhol / microlayerContext.sigma * (Foam::pow(Rbdot / Rb, 2.0) - Rbddot / (3 * Rb));

                const scalar a2 = -microlayerContext.rhol * Rbddot / (2 * microlayerContext.sigma);

                const scalar a1 = 1 / Rb;

                const scalar a0 = -1;

                struct CubicFunction
                {
                    Foam::scalar a, b, c, d;

                    std::pair<Foam::scalar, Foam::scalar> operator()(Foam::scalar x) const
                    {
                        Foam::scalar fx = a * x * x * x + b * x * x + c * x + d;
                        Foam::scalar dfx = 3 * a * x * x + 2 * b * x + c;
                        return std::make_pair(fx, dfx);
                    }
                };

                CubicFunction f{a3, a2, a1, a0};

                const scalar guess = 10e-6;
                const scalar min = -100.0;
                const scalar max = 100.0;
                const int digits = std::numeric_limits<Foam::scalar>::digits;

                scalar root = boost::math::tools::newton_raphson_iterate(
                    f, guess, min, max, digits);

                const scalar xbar = root;

                //const scalar temp = a3 * Foam::pow(Rb, 3) + a2 * Foam::pow(Rb, 2) + a1 * Rb + a0;

                microlayerContext.curv = 3 * a3 * Foam::pow(xbar, 2) + 2 * a2 * xbar + a1;

                microlayerContext.rc = 1.0 / microlayerContext.curv;

                microlayerContext.del0 = 1.34 * microlayerContext.rc * Foam::pow(microlayerContext.mul * Rbdot / microlayerContext.sigma, 2.0 / 3.0);

                microlayerContext.delM = microlayerContext.del0; //0.3*microlayerContext.rc; //10 * microlayerContext.del0;

                Info << "Initial microlayer thickness = " << microlayerContext.del0 * 1e6 << " micrometers" << nl;

            }

            microlayerDict.set("meniscusCurvature", microlayerContext.curv);
            microlayerDict.set("bubbleFootRadiusOfCurvature", microlayerContext.rc);
            microlayerDict.set("initialMicrolayerThickness", microlayerContext.del0);

            Info << "Remeshing ..." << nl;

            microlayerContext.n = max(int(microlayerContext.l / microlayerContext.dsReq) + 1, microlayerContext.nMin); // Have at least 11 nodes

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
                if (j == 0)
                {
                    hNewp[j] = 0;
                }

                if (j == microlayerContext.n - 1)
                {
                    // hNewp[j] = 2 * hNewp[j - 1] - 5 * hNewp[j - 2] + 4 * hNewp[j - 3] - hNewp[j - 4] - microlayerContext.ds * microlayerContext.ds * microlayerContext.curv;
                    hNewp[j] = microlayerContext.delM; // Matching thickness at end of microlayer
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
            PetscCall(SNESSetTolerances(snes, 1e-10, 1e-8, 1e-12, 100, 1000));
            PetscCall(SNESGetTolerances(snes, &abstol, &rtol, &stol, &maxit, &maxf));
            PetscCall(SNESSetFunction(snes, r, FormFunction, &microlayerContext));
            PetscCall(SNESSetJacobian(snes, J, J, FormJacobian, &microlayerContext));

            Info << "Solving microlayer model ..." << nl;
            PetscCall(SNESSolve(snes, NULL, h));
            PetscCall(SNESGetIterationNumber(snes, &niter));
            PetscCall(PetscPrintf(PETSC_COMM_WORLD, "number of SNES iterations = %" PetscInt_FMT "\n\n", niter));

            PetscCall(VecCopy(h, microlayerContext.h));

            PetscCall(VecDestroy(&h));
            PetscCall(VecDestroy(&r));
            PetscCall(MatDestroy(&J));
            PetscCall(SNESDestroy(&snes));
        }
    }

    computeHeatFlux(&microlayerContext);

    // Write files
    if (runTime.writeTime())
    {
        PetscViewer petscviewer;

        Info << "Writing microlayer model results at time " << runTime.timeName() << nl;

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

        fileName qml_file = writeDir / "qml.csv";

        PetscCall(PetscViewerASCIIOpen(PETSC_COMM_WORLD, qml_file.c_str(), &petscviewer));
        PetscCall(VecView(microlayerContext.q, petscviewer));
        PetscCall(PetscViewerDestroy(&petscviewer));
    }
}