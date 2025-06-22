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
    PetscScalar sigma = ctx->sigma; // Surface tension
    PetscScalar mul = ctx->mul;     // Viscosity
    PetscScalar rhol = ctx->rhol;   // Liquid density
    PetscScalar Ri = ctx->Ri;       // Interface Heat Transfer Resistance
    PetscScalar kappal = ctx->kappal; // Thermal conductivity
    PetscScalar Tsat = ctx->Tsat;   // Saturation temperature
    PetscScalar hfg = ctx->hfg;     // Heat of vaporization

    // For now, assume constant wall temperature
    PetscScalar Tw = ctx->Tw; // Wall temperature
    // Model parameters
    PetscScalar ds = ctx->ds;            // Grid size
    PetscScalar dt = ctx->dt;            // Time step
    PetscScalar delm = ctx->delm;        // Microregion height
    PetscScalar delM = ctx->delM; // Matching thickness at end of microlayer
    PetscScalar rc = ctx->rc;            // Bubble foot radius of curvature
    PetscScalar ls = ctx->ls;            // Slip length
    PetscScalar ucl = ctx->ucl;          // Contact line speed
    PetscScalar theta = ctx->thetaApp;   // Apparent contact angle

    // Current State
    Vec hn = ctx->h; // Current interface profile

    PetscFunctionBeginUser;
    PetscCall(VecGetArrayRead(h, &hp));
    PetscCall(VecGetArrayRead(hn, &hnp));
    PetscCall(VecGetArray(f, &fp));

    // Compute function

    // Boundary conditions at Contact Line
    fp[0] = hp[0] - delm;                             // microlayer height = microregion height
    fp[1] = hp[1] - hp[0] - ds * tan(theta);          // apparent contact angle condition
    fp[2] = (-hp[0] + 3 * hp[1] - 3 * hp[2] + hp[3]); // dp/ds = 0 at contact line

    // Boundary conditions at Bubble foot
    fp[n - 1] = hp[n - 1] - delM;                                         // microlayer height = matching thickness at end of microlayer
    fp[n - 2] = (hp[n - 1] - 2 * hp[n - 2] + hp[n - 3]) - (ds * ds) / rc; // curvature condition at bubble foot
    fp[n - 3] = (hp[n - 1] - 3 * hp[n - 2] + 3 * hp[n - 3] - hp[n - 4]);  // dp/ds = 0 at bubble foot

    // Interior points
    for (i = 3; i < n - 3; i++)
    {
        // fp[i] = hp[i] - hnp[i] + dt / (6 * mu * ds) * (hp[i + 1] * hp[i + 1] * (hp[i + 1] + 3 * Ls) * (sigma * (hp[i + 3] - 2 * hp[i + 2] + 2 * hp[i] - hp[i - 1]) / (2 * ds * ds * ds) + rhol * g) - hp[i - 1] * hp[i - 1] * (hp[i - 1] + 3 * Ls) * (sigma * (hp[i + 1] - 2 * hp[i] + 2 * hp[i - 2] - hp[i - 3]) / (2 * ds * ds * ds) + rhol * g)) - dt * Ucl * (hp[i + 1] - hp[i - 1]) / (2 * ds) + dt * (Tw - Tsat) / (rhol * hfg * (Ri + hp[i] / k)) - dt * sigma * Tsat / (pow(rhol * hfg, 2) * (Ri + hp[i] / k)) * (hp[i + 1] - 2 * hp[i] + hp[i - 1]) / (ds * ds);

        fp[i] = hp[i] - hnp[i] +
                dt / (2 * ds) * (sigma / (3 * mul) * hp[i + 1] * hp[i + 1] * (hp[i + 1] + 3 * ls) * (hp[i + 3] - 2 * hp[i + 2] + 2 * hp[i] - hp[i - 1]) / (2 * ds * ds * ds) - ucl * hp[i + 1]) -
                dt / (2 * ds) * (sigma / (3 * mul) * hp[i - 1] * hp[i - 1] * (hp[i - 1] + 3 * ls) * (hp[i + 1] - 2 * hp[i] + 2 * hp[i - 2] - hp[i - 3]) / (2 * ds * ds * ds) - ucl * hp[i - 1]) +
                dt * (Tw - Tsat) / (rhol * hfg * (Ri + hp[i] / kappal));
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
    const PetscScalar *hp; // hp : pointer to next interface profile (variable of nonlinear eqns)

    PetscInt i, j[7], n; // i : Row Index, j : Column Index
    PetscCall(VecGetSize(h, &n));
    PetscScalar A[7]; // Jacobian entries

    PetscCall(VecGetSize(h, &n));

    // Fluid properties
    PetscScalar sigma = ctx->sigma; // Surface tension
    PetscScalar mul = ctx->mul;     // Viscosity
    PetscScalar rhol = ctx->rhol;   // Liquid density
    PetscScalar Ri = ctx->Ri;       // Interface Heat Transfer Resistance
    PetscScalar kappal = ctx->kappal; // Thermal conductivity
    PetscScalar Tsat = ctx->Tsat;   // Saturation temperature
    PetscScalar hfg = ctx->hfg;     // Heat of vaporization

    // For now, assume constant wall temperature
    PetscScalar Tw = ctx->Tw; // Wall temperature
    // Model parameters
    PetscScalar ds = ctx->ds;            // Grid size
    PetscScalar dt = ctx->dt;            // Time step
    PetscScalar delm = ctx->delm;        // Microregion height
    PetscScalar delM = ctx->delM; // Matching thickness at end of microlayer
    PetscScalar rc = ctx->rc;            // Bubble foot radius of curvature
    PetscScalar ls = ctx->ls;            // Slip length
    PetscScalar ucl = ctx->ucl;          // Contact line speed
    PetscScalar theta = ctx->thetaApp;   // Apparent contact angle

    // Current State
    Vec hn = ctx->h; // Current interface profile

    PetscFunctionBeginUser;
    PetscCall(VecGetArrayRead(h, &hp));

    // Interior Grid Points
    for (i = 3; i < n - 3; i++)
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

        A[2] = - dt / (4 * ds * ds * ds * ds) * sigma / mul * (hp[i - 1] * (hp[i - 1] + 2 * ls) * (hp[i + 1] - 2 * hp[i] + 2 * hp[i - 2] - hp[i - 3]) + hp[i + 1] * hp[i + 1] / 3 * (hp[i + 1] + 3 * ls)) + dt * ucl / (2 * ds);

        A[3] = 1 + dt / (6 * ds * ds * ds * ds) * sigma / mul * (hp[i - 1] * hp[i - 1] * (hp[i - 1] + 3 * ls) + hp[i + 1] * hp[i + 1] * (hp[i + 1] + 3 * ls)) - dt * (Tw - Tsat) / (rhol * hfg) * (1 / (Ri + hp[i] / kappal)) * (1 / (Ri + hp[i] / kappal)) * 1 / kappal ;

        A[5] = -dt / (6 * ds * ds * ds * ds) * sigma / mul * hp[i + 1] * hp[i + 1] * (hp[i + 1] + 3 * ls);

        A[6] = dt / (12 * ds * ds * ds * ds) * sigma / mul * hp[i + 1] * hp[i + 1] * (hp[i + 1] + 3 * ls);

        PetscCall(MatSetValues(B, 1, &i, 7, j, A, INSERT_VALUES));
    }

    // Boundary Conditions at Contact Line
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

    // Boundary Conditions at Bubble Foot
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

// ************************************************************************* //