#ifndef PETSC_FUNCTIONS_H
#define PETSC_FUNCTIONS_H

#include <petscsnes.h>
#include "volFields.H"
#include "fvMesh.H"
#include <boost/math/interpolators/makima.hpp>

// Define PETSC Context struct

struct MicrolayerContext
{   
    PetscBool isMl; // Flag to indicate if microlayer is present
    // contactline position, microlayer length, desired grid size, actual grid size time step
    PetscScalar xCl, l, dsReq, ds, dt; 
    PetscInt n; // number of grid points
    Vec s, x, h, q; // grid in contactline frame, grid in global frame, microlayer thickness, heat flux
    
    PetscScalar ub, ucl, thetaApp, delm, del0, rc, delM; // bubble expansion velocity, contactline speed, apparent contact angle, microregion height, initial microlayer thickness, bubble foot radius of curvature, matching thickness at end of microlayer

    // Fluid properties
    PetscScalar rhol, mul, kappal, Tsat, sigma, hfg, Ri, ls; // density, viscosity, thermal conductivity, saturation temperature, surface tension, latent heat of vaporization, and interface thermal resistance
    // Constructor

    // Wall temperature (assumed constant for now)
    PetscScalar Tw;

    MicrolayerContext(Foam::IOdictionary& microlayerDict)
    {
        // Initialize PETSC context with microlayer properties
        xCl = microlayerDict.lookupOrDefault<Foam::scalar>("initialContactLinePosition", 0.0);
        l = microlayerDict.get<Foam::scalar>("microlayerLength");
        if (l > 0)
        {
            isMl = PETSC_TRUE; // Microlayer is present
        }

        dsReq = microlayerDict.lookupOrDefault<Foam::scalar>("desiredGridSize", 1e-6);
        n = int(l / dsReq) + 1;
        ds =  l / (n - 1);  

        dt = 0.0; // will be set during iteration
        ub = 0.0; // will be set during iteration
        
        ucl = microlayerDict.get<Foam::scalar>("contactLineVelocity");
        thetaApp = microlayerDict.get<Foam::scalar>("thetaApp");
        delm = microlayerDict.lookupOrDefault<Foam::scalar>("microregionHeight", 1e-7);
        del0 = microlayerDict.lookupOrDefault<Foam::scalar>("initialMicrolayerThickness", 1e-6);
        rc = microlayerDict.lookupOrDefault<Foam::scalar>("bubbleFootRadiusOfCurvature", 1e-6);
        delM = 0.3 * rc; 

        rhol = microlayerDict.get<Foam::scalar>("rhol");
        mul = microlayerDict.get<Foam::scalar>("mul");
        kappal = microlayerDict.get<Foam::scalar>("kappal");
        Tsat = microlayerDict.get<Foam::scalar>("Tsat");
        sigma = microlayerDict.get<Foam::scalar>("sigma");
        hfg = microlayerDict.get<Foam::scalar>("hfg");
        Ri = microlayerDict.get<Foam::scalar>("Ri");
        ls = microlayerDict.get<Foam::scalar>("ls");
        
        Tw = microlayerDict.get<Foam::scalar>("Tw"); 

        // Initialize PETSC vectors
        VecCreateSeq(PETSC_COMM_SELF, n, &s);
        VecCreateSeq(PETSC_COMM_SELF, n, &x);
        VecCreateSeq(PETSC_COMM_SELF, n, &h);

        PetscScalar *sp, *xp, *hp;
        VecGetArray(s, &sp);
        VecGetArray(x, &xp);
        VecGetArray(h, &hp);

        std::vector<double> xSp = {0, 0.1 * l, 0.2 * l, 0.22 * l, 0.3 * l, 0.4 * l, 0.5 * l, 0.6 * l, 0.7 * l, 0.8 * l, 1.001*l};
        std::vector<double> ySp = {0, 1.05 * delM, 0.52 * delM, 0.35 * delM, 0.45 * delM, 0.59 * delM, 0.73 * delM, 0.86 * delM, 0.95 * delM, 1.25 * delM, 1.5 * delM};

        auto spline = boost::math::interpolators::makima<std::vector<double>>(std::move(xSp), std::move(ySp));

        for (PetscInt i = 0; i < n; i++)
        {
            sp[i] = i * ds;
            xp[i] = xCl + i * ds;
            hp[i] = spline(i * ds); //delM * i * ds / l; // Initial microlayer thickness (wedge shape)
        }

        VecRestoreArray(s, &sp);
        VecRestoreArray(x, &xp);
        VecRestoreArray(h, &hp);
    }
};

PetscErrorCode FormFunction(SNES snes, Vec h, Vec f, void *ctx_);

PetscErrorCode FormJacobian(SNES snes, Vec h, Mat jac, Mat B, void *ctx_);

PetscErrorCode computeHeatFlux(void *ctx_);

PetscScalar computeMeanHeatFlux(PetscScalar xMin, PetscScalar xMax, void *ctx_);

#endif

// ************************************************************************* //