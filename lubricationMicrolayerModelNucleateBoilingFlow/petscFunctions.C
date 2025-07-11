#include "petscFunctions.H"

PetscErrorCode FormFunction(SNES snes, Vec h, Vec f, void *ctx_)
{
    MicrolayerContext *ctx = (MicrolayerContext *)ctx_;

    // Pointers to PETSc vectors
    const PetscScalar *hp, *hnp; // hp : pointer to next interface profile (variable of nonlinear eqns), hnp : pointer to current interface profile (stored in struct)
    PetscScalar *fp;             // ff : pointer to residual vector
    PetscInt i, n;
    PetscCall(VecGetSize(h, &n));

    // Fluid properties
    PetscScalar sigma = ctx->sigma;   // Surface tension
    PetscScalar mul = ctx->mul;       // Viscosity
    PetscScalar rhol = ctx->rhol;     // Liquid density
    PetscScalar Ri = ctx->Ri;         // Interface Heat Transfer Resistance
    PetscScalar kappal = ctx->kappal; // Thermal conductivity
    PetscScalar Tsat = ctx->Tsat;     // Saturation temperature
    PetscScalar hfg = ctx->hfg;       // Heat of vaporization

    // For now, assume constant wall temperature
    PetscScalar Tw = ctx->Tw; // Wall temperature
    // Model parameters
    PetscScalar ds = ctx->ds;          // Grid size
    PetscScalar dt = ctx->dt;          // Time step
    PetscScalar delm = ctx->delm;      // Microregion height
    PetscScalar delM = ctx->delM;      // Matching thickness at end of microlayer
    PetscScalar rc = ctx->rc;          // Bubble foot radius of curvature
    PetscScalar ls = ctx->ls;          // Slip length
    PetscScalar ucl = ctx->ucl;        // Contact line speed
    PetscScalar theta = ctx->thetaApp; // Apparent contact angle
    PetscScalar curv = ctx->curv;      // Meniscus curvature
    // Current State
    Vec hn = ctx->h; // Current interface profile

    PetscFunctionBeginUser;
    PetscCall(VecGetArrayRead(h, &hp));
    PetscCall(VecGetArrayRead(hn, &hnp));
    PetscCall(VecGetArray(f, &fp));

    // Compute function

    // Boundary conditions at Contact Line
    /*
    fp[0] = hp[0] - delm;                             // microlayer height = microregion height
    fp[1] = hp[1] - hp[0] - ds * tan(theta);          // apparent contact angle condition
    fp[2] = (-hp[0] + 3 * hp[1] - 3 * hp[2] + hp[3]); // dp/ds = 0 at contact line
    */

    // Higher order boundary conditions at contact line
    fp[0] = hp[0] - delm;                                                    // microlayer height = microregion height
    fp[1] = -1.5 * hp[0] + 2 * hp[1] - 0.5 * hp[2] - ds * tan(theta);        // apparent contact angle condition
    fp[2] = -2.5 * hp[0] + 9 * hp[1] - 12 * hp[2] + 7 * hp[3] - 1.5 * hp[4]; // dp/ds = 0 at contact line

    // Boundary conditions at Bubble foot
    /*
    fp[n - 1] = hp[n - 1] - delM;                                         // microlayer height = matching thickness at end of microlayer
    fp[n - 2] = (hp[n - 1] - 2 * hp[n - 2] + hp[n - 3]) - (ds * ds) * curv; // curvature condition at bubble foot
    fp[n - 3] = (hp[n - 1] - 3 * hp[n - 2] + 3 * hp[n - 3] - hp[n - 4]);  // dp/ds = 0 at bubble foot
    */

    // Higher order boundary conditions at bubble foot
    fp[n - 1] = hp[n - 1] - delM;                                                                   // microlayer height = matching thickness at end of microlayer
    fp[n - 2] = 2 * hp[n - 1] - 5 * hp[n - 2] + 4 * hp[n - 3] - hp[n - 4] - ds * ds * curv;         // curvature condition at bubble foot
    fp[n - 3] = 2.5 * hp[n - 1] - 9 * hp[n - 2] + 12 * hp[n - 3] - 7 * hp[n - 4] + 1.5 * hp[n - 5]; // dp/ds = 0 at bubble foot

    // Interior points
    for (i = 3; i < n - 3; i++)
    {
        // fp[i] = hp[i] - hnp[i] + dt / (6 * mu * ds) * (hp[i + 1] * hp[i + 1] * (hp[i + 1] + 3 * Ls) * (sigma * (hp[i + 3] - 2 * hp[i + 2] + 2 * hp[i] - hp[i - 1]) / (2 * ds * ds * ds) + rhol * g) - hp[i - 1] * hp[i - 1] * (hp[i - 1] + 3 * Ls) * (sigma * (hp[i + 1] - 2 * hp[i] + 2 * hp[i - 2] - hp[i - 3]) / (2 * ds * ds * ds) + rhol * g)) - dt * Ucl * (hp[i + 1] - hp[i - 1]) / (2 * ds) + dt * (Tw - Tsat) / (rhol * hfg * (Ri + hp[i] / k)) - dt * sigma * Tsat / (pow(rhol * hfg, 2) * (Ri + hp[i] / k)) * (hp[i + 1] - 2 * hp[i] + hp[i - 1]) / (ds * ds);

        // Central difference for first derivative

        //- dt * sigma * Tsat / (PetscPowScalar(rhol * hfg, 2) * (Ri + hp[i] / kappal)) * (hp[i + 1] - 2 * hp[i] + hp[i - 1]) / (ds * ds);

        if (hnp[i] <= 1e-10)
        {
            fp[i] = hp[i];
        }

        else
        {
            fp[i] = hp[i] - hnp[i] + dt / (2 * ds) * (sigma / (3 * mul) * hp[i + 1] * hp[i + 1] * (hp[i + 1] + 3 * ls) * (hp[i + 3] - 2 * hp[i + 2] + 2 * hp[i] - hp[i - 1]) / (2 * ds * ds * ds) - ucl * hp[i + 1]) - dt / (2 * ds) * (sigma / (3 * mul) * hp[i - 1] * hp[i - 1] * (hp[i - 1] + 3 * ls) * (hp[i + 1] - 2 * hp[i] + 2 * hp[i - 2] - hp[i - 3]) / (2 * ds * ds * ds) - ucl * hp[i - 1]) + dt * (Tw - Tsat) / (rhol * hfg * (Ri + hp[i] / kappal));
        }

        // Upwinded first derivative
        /*
        fp[i] = hp[i] - hnp[i]
                + dt / (2 * ds) * (sigma / (3 * mul) * hp[i + 1] * hp[i + 1] * (hp[i + 1] + 3 * ls) * (hp[i + 3] - 2 * hp[i + 2] + 2 * hp[i] - hp[i - 1]) / (2 * ds * ds * ds))
                - dt / (2 * ds) * (sigma / (3 * mul) * hp[i - 1] * hp[i - 1] * (hp[i - 1] + 3 * ls) * (hp[i + 1] - 2 * hp[i] + 2 * hp[i - 2] - hp[i - 3]) / (2 * ds * ds * ds))
                - dt * ucl * (hp[i + 1] - hp[i]) / ds
                + dt * (Tw - Tsat) / (rhol * hfg * (Ri + hp[i] / kappal)) ;
                //- dt * sigma * Tsat / (PetscPowScalar(rhol * hfg, 2) * (Ri + hp[i] / kappal)) / (ds * ds);
        */
    }

    // Restore vectors
    PetscCall(VecRestoreArrayRead(h, &hp));
    PetscCall(VecRestoreArrayRead(hn, &hnp));
    PetscCall(VecRestoreArray(f, &fp));
    PetscFunctionReturn(PETSC_SUCCESS);
}

// Function to calculate the Jacobian
PetscErrorCode FormJacobian(SNES snes, Vec h, Mat jac, Mat B, void *ctx_)
{
    MicrolayerContext *ctx = (MicrolayerContext *)ctx_;

    // Pointers to PETSc vectors
    const PetscScalar *hp, *hnp; // hp : pointer to next interface profile (variable of nonlinear eqns)

    PetscInt i, j[7], n; // i : Row Index, j : Column Index
    PetscCall(VecGetSize(h, &n));
    PetscScalar A[7]; // Jacobian entries

    PetscCall(VecGetSize(h, &n));

    // Fluid properties
    PetscScalar sigma = ctx->sigma;   // Surface tension
    PetscScalar mul = ctx->mul;       // Viscosity
    PetscScalar rhol = ctx->rhol;     // Liquid density
    PetscScalar Ri = ctx->Ri;         // Interface Heat Transfer Resistance
    PetscScalar kappal = ctx->kappal; // Thermal conductivity
    PetscScalar Tsat = ctx->Tsat;     // Saturation temperature
    PetscScalar hfg = ctx->hfg;       // Heat of vaporization

    // For now, assume constant wall temperature
    PetscScalar Tw = ctx->Tw; // Wall temperature
    // Model parameters
    PetscScalar ds = ctx->ds; // Grid size
    PetscScalar dt = ctx->dt; // Time step

    PetscScalar ls = ctx->ls;   // Slip length
    PetscScalar ucl = ctx->ucl; // Contact line speed

    PetscFunctionBeginUser;
    PetscCall(VecGetArrayRead(h, &hp));
    PetscCall(VecGetArrayRead(ctx->h, &hnp));

    // Interior Grid Points
    for (i = 3; i < n - 3; i++)
    {

        // Central difference for first derivative

        if (hnp[i] > 1e-10)
        {
	    j[0] = i - 3;
            j[1] = i - 2;
            j[2] = i - 1;
            j[3] = i;
            j[4] = i + 1;
            j[5] = i + 2;
            j[6] = i + 3;

            A[0] = dt / (12 * ds * ds * ds * ds) * sigma / mul * hp[i - 1] * hp[i - 1] * (hp[i - 1] + 3 * ls);

            A[1] = -dt / (6 * ds * ds * ds * ds) * sigma / mul * hp[i - 1] * hp[i - 1] * (hp[i - 1] + 3 * ls);

            A[2] = -dt / (4 * ds * ds * ds * ds) * sigma / mul * (hp[i - 1] * (hp[i - 1] + 2 * ls) * (hp[i + 1] - 2 * hp[i] + 2 * hp[i - 2] - hp[i - 3]) + hp[i + 1] * hp[i + 1] / 3 * (hp[i + 1] + 3 * ls)) + dt * ucl / (2 * ds);
            //-dt * sigma * Tsat / (PetscPowScalar(rhol * hfg, 2) * (Ri + hp[i] / kappal)) / (ds * ds);

            A[3] = 1 + dt / (6 * ds * ds * ds * ds) * sigma / mul * (hp[i - 1] * hp[i - 1] * (hp[i - 1] + 3 * ls) + hp[i + 1] * hp[i + 1] * (hp[i + 1] + 3 * ls)) - dt * (Tw - Tsat) / (rhol * hfg) * (1 / (Ri + hp[i] / kappal)) * (1 / (Ri + hp[i] / kappal)) * 1 / kappal;
            //+ dt * Tsat * sigma / (rhol * rhol * hfg * hfg) * (1 / (Ri + hp[i] / kappal)) * (1 / (Ri + hp[i] / kappal)) / kappal  * (hp[i + 1] - 2 * hp[i] + hp[i - 1]) / (ds * ds);

            A[4] = dt / (4 * ds * ds * ds * ds) * sigma / mul * (hp[i + 1] * (hp[i + 1] + 2 * ls) * (hp[i + 3] - 2 * hp[i + 2] + 2 * hp[i] - hp[i - 1]) - hp[i - 1] * hp[i - 1] / 3 * (hp[i - 1] + 3 * ls)) - dt * ucl / (2 * ds);
            //-dt * sigma * Tsat / (PetscPowScalar(rhol * hfg, 2) * (Ri + hp[i] / kappal)) / (ds * ds);

            A[5] = -dt / (6 * ds * ds * ds * ds) * sigma / mul * hp[i + 1] * hp[i + 1] * (hp[i + 1] + 3 * ls);

            A[6] = dt / (12 * ds * ds * ds * ds) * sigma / mul * hp[i + 1] * hp[i + 1] * (hp[i + 1] + 3 * ls);

            // Upwinded first derivative
            /*
            A[0] = dt / (12 * ds * ds * ds * ds) * sigma / mul * hp[i - 1] * hp[i - 1] * (hp[i - 1] + 3 * ls);

            A[1] = -dt / (6 * ds * ds * ds * ds) * sigma / mul * hp[i - 1] * hp[i - 1] * (hp[i - 1] + 3 * ls);

            A[2] = -dt / (4 * ds * ds * ds * ds) * sigma / mul * (hp[i - 1] * (hp[i - 1] + 2 * ls) * (hp[i + 1] - 2 * hp[i] + 2 * hp[i - 2] - hp[i - 3]) + hp[i + 1] * hp[i + 1] / 3 * (hp[i + 1] + 3 * ls));
                   //-dt * sigma * Tsat / (PetscPowScalar(rhol * hfg, 2) * (Ri + hp[i] / kappal)) / (ds * ds);

            A[3] = 1 + dt / (6 * ds * ds * ds * ds) * sigma / mul * (hp[i - 1] * hp[i - 1] * (hp[i - 1] + 3 * ls) + hp[i + 1] * hp[i + 1] * (hp[i + 1] + 3 * ls)) - dt * (Tw - Tsat) / (rhol * hfg) * (1 / (Ri + hp[i] / kappal)) * (1 / (Ri + hp[i] / kappal)) * 1 / kappal + dt * ucl / ds;
                   //+ dt * Tsat * sigma / (rhol * rhol * hfg * hfg) * (1 / (Ri + hp[i] / kappal)) * (1 / (Ri + hp[i] / kappal)) / kappal  * (hp[i + 1] - 2 * hp[i] + hp[i - 1]) / (ds * ds);

            A[4] = dt / (4 * ds * ds * ds * ds) * sigma / mul * (hp[i +1] * (hp[i + 1] + 2 * ls) * (hp[i + 3] - 2 * hp[i + 2] + 2 * hp[i] - hp[i - 1]) - hp[i - 1] * hp[i - 1] / 3 * (hp[i - 1] + 3 * ls)) - dt * ucl / ds;
                   //-dt * sigma * Tsat / (PetscPowScalar(rhol * hfg, 2) * (Ri + hp[i] / kappal)) / (ds * ds);

            A[5] = -dt / (6 * ds * ds * ds * ds) * sigma / mul * hp[i + 1] * hp[i + 1] * (hp[i + 1] + 3 * ls);

            A[6] = dt / (12 * ds * ds * ds * ds) * sigma / mul * hp[i + 1] * hp[i + 1] * (hp[i + 1] + 3 * ls);
            */

            PetscCall(MatSetValues(B, 1, &i, 7, j, A, INSERT_VALUES));
        }

        else
        {
            j[0] = i;

            A[0] = 1.0;

            PetscCall(MatSetValues(B, 1, &i, 1, j, A, INSERT_VALUES));
        }
    }

    // Boundary Conditions at Contact Line

    /*
    i = 0;
    A[0] = 1.0;

    PetscCall(MatSetValues(B, 1, &i, 1, &i, A, INSERT_VALUES));

    i = 1;
    j[0] = 0;
    j[1] = 1;

    A[0] = -1.0;
    A[1] = 1.0;

    PetscCall(MatSetValues(B, 1, &i, 2, j, A, INSERT_VALUES));

    i = 2;
    j[0] = 0;
    j[1] = 1;
    j[2] = 2;
    j[3] = 3;

    A[0] = -1.0;
    A[1] = 3.0;
    A[2] = -3.0;
    A[3] = 1.0;

    PetscCall(MatSetValues(B, 1, &i, 4, j, A, INSERT_VALUES));
    */

    // Higher order boundary conditions at contact line
    i = 0;
    A[0] = 1.0;

    PetscCall(MatSetValues(B, 1, &i, 1, &i, A, INSERT_VALUES));

    i = 1;
    j[0] = 0;
    j[1] = 1;
    j[2] = 2;

    A[0] = -3.0 / 2.0;
    A[1] = 2.0;
    A[2] = -1.0 / 2.0;

    PetscCall(MatSetValues(B, 1, &i, 3, j, A, INSERT_VALUES));

    i = 2;
    j[0] = 0;
    j[1] = 1;
    j[2] = 2;
    j[3] = 3;
    j[4] = 4;

    A[0] = -5.0 / 2.0;
    A[1] = 9.0;
    A[2] = -12.0;
    A[3] = 7.0;
    A[4] = -3.0 / 2.0;

    PetscCall(MatSetValues(B, 1, &i, 5, j, A, INSERT_VALUES));

    // Boundary Conditions at Bubble Foot
    /*
    i = n - 1;
    A[0] = 1.0;

    PetscCall(MatSetValues(B, 1, &i, 1, &i, A, INSERT_VALUES));

    i = n - 2;
    j[0] = n - 3;
    j[1] = n - 2;
    j[2] = n - 1;

    A[0] = 1.0;
    A[1] = -2.0;
    A[2] = 1.0;

    PetscCall(MatSetValues(B, 1, &i, 3, j, A, INSERT_VALUES));

    i = n - 3;
    j[0] = n - 4;
    j[1] = n - 3;
    j[2] = n - 2;
    j[3] = n - 1;

    A[0] = -1.0;
    A[1] = 3.0;
    A[2] = -3.0;
    A[3] = 1.0;

    PetscCall(MatSetValues(B, 1, &i, 4, j, A, INSERT_VALUES));
    */

    // Higher order boundary conditions at bubble foot
    i = n - 1;
    A[0] = 1.0;

    PetscCall(MatSetValues(B, 1, &i, 1, &i, A, INSERT_VALUES));

    i = n - 2;
    j[0] = n - 4;
    j[1] = n - 3;
    j[2] = n - 2;
    j[3] = n - 1;

    A[0] = -1.0;
    A[1] = 4.0;
    A[2] = -5.0;
    A[3] = 2.0;

    PetscCall(MatSetValues(B, 1, &i, 4, j, A, INSERT_VALUES));

    i = n - 3;
    j[0] = n - 5;
    j[1] = n - 4;
    j[2] = n - 3;
    j[3] = n - 2;
    j[4] = n - 1;

    A[0] = 3.0 / 2.0;
    A[1] = -7.0;
    A[2] = 12.0;
    A[3] = -9.0;
    A[4] = 5.0 / 2.0;

    PetscCall(MatSetValues(B, 1, &i, 5, j, A, INSERT_VALUES));

    // Restore vectors
    PetscCall(VecRestoreArrayRead(h, &hp));

    // Assemble Matrix
    PetscCall(MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY));
    if (jac != B)
    {
        PetscCall(MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY));
        PetscCall(MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY));
    }
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode boundThickness(void *ctx_, PetscScalar lim_)
{
    MicrolayerContext *ctx = (MicrolayerContext *)ctx_;
    // Fluid properties

    // Current State
    Vec h_ = ctx->h; // Current interface profile

    // Get access to the array data
    PetscScalar *hp;
    PetscInt n;

    PetscCall(VecGetSize(h_, &n));
    PetscCall(VecGetArray(h_, &hp));

    for (PetscInt i = 1; i < n; ++i)
    {
        if (hp[i] < lim_)
        {
            // If the microlayer thickness is below the limit, set it to the limit
            hp[i] = 0;
        }
    }

    VecRestoreArray(h_, &hp);

    return PETSC_SUCCESS;
}

PetscErrorCode computeHeatFlux(void *ctx_)
{
    MicrolayerContext *ctx = (MicrolayerContext *)ctx_;
    // Fluid properties
    PetscScalar Ri = ctx->Ri;         // Interface Heat Transfer Resistance
    PetscScalar kappal = ctx->kappal; // Thermal conductivity
    PetscScalar Tsat = ctx->Tsat;     // Saturation temperature

    // For now, assume constant wall temperature
    PetscScalar Tw = ctx->Tw; // Wall temperature

    // Current State
    Vec h_ = ctx->h; // Current interface profile

    VecDuplicate(h_, &ctx->q);

    Vec q_ = ctx->q; // Heat flux vector

    // Get access to the array data
    const PetscScalar *hp;
    PetscScalar *qp;
    PetscInt n;

    PetscCall(VecGetSize(h_, &n));
    PetscCall(VecGetArrayRead(h_, &hp));
    PetscCall(VecGetArray(q_, &qp));

    for (PetscInt i = 0; i < n; ++i)
    {
        
        if (hp[i] <= 1e-10)
        {
            // If the microlayer thickness is too small, set the heat flux to zero
            qp[i] = 0.0;
        }

        else
        {
            qp[i] = (Tw - Tsat) / (Ri + hp[i] / kappal);
        }
    }

    VecRestoreArrayRead(h_, &hp);
    VecRestoreArray(q_, &qp);

    return PETSC_SUCCESS;
}

PetscScalar computeMeanHeatFlux(PetscScalar xMin, PetscScalar xMax, void *ctx_)
{

    MicrolayerContext *ctx = (MicrolayerContext *)ctx_;

    if (!ctx->isMl)
    {
        // If microlayer is not present, return 0
        return 0.0;
    }

    const PetscScalar *xp, *qp;
    Vec q = ctx->q;
    Vec x = ctx->x;
    PetscInt n;

    PetscCall(VecGetSize(q, &n));
    PetscCall(VecGetArrayRead(x, &xp));
    PetscCall(VecGetArrayRead(q, &qp));

    PetscScalar integral = 0.0;
    PetscScalar length = 0.0;
    PetscBool first = PETSC_TRUE;

    if (xMax <= xp[0] || xMin >= xp[n - 1])
    {
        // If the range is outside the grid, return 0
        PetscCall(VecRestoreArrayRead(q, &qp));
        PetscCall(VecRestoreArrayRead(x, &xp));
        return 0.0;
    }

    for (PetscInt i = 0; i < n; ++i)
    {
        if (xp[i] >= xMin && xp[i] <= xMax)
        {
            if (first)
            {
                PetscScalar dx = xp[i] - xMin;
                PetscScalar avgQ = qp[i];
                PetscScalar avgx = xp[i];
                integral += avgQ * avgx * dx;
                length += dx;
                first = PETSC_FALSE;
            }

            else
            {
                PetscScalar dx = xp[i] - xp[i - 1];
                PetscScalar avgQ = (qp[i] + qp[i - 1]) / 2.0; // Trapezoidal rule
                PetscScalar avgx = (xp[i] + xp[i - 1]) / 2.0; // Midpoint for average
                integral += 2 * avgQ * avgx * dx;             // Integral over circular arc segment
                length += dx;
            }
        }
    }

    PetscScalar denom = xMax * xMax - xMin * xMin;

    PetscCall(VecRestoreArrayRead(q, &qp));
    PetscCall(VecRestoreArrayRead(x, &xp));

    return integral / denom; // Mean flux
}
// ************************************************************************* //
