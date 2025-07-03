#ifndef PETSC_FUNCTIONS_H
#define PETSC_FUNCTIONS_H

#include <petscsnes.h>
#include "volFields.H"
#include "fvMesh.H"
#include <boost/math/interpolators/makima.hpp>
#include <boost/math/tools/roots.hpp>

// Define PETSC Context struct

struct MicrolayerContext
{
    PetscBool isMl; // Flag to indicate if microlayer is present
    // contactline position, microlayer length, desired grid size, actual grid size time step
    PetscScalar xCl, l, dsReq, ds, dt;
    PetscInt n, nMin; // number of grid points
    Vec s, x, h, q;   // grid in contactline frame, grid in global frame, microlayer thickness, heat flux

    PetscScalar ub, ucl, thetaApp, delm, del0, rc, delM; // bubble expansion velocity, contactline speed, apparent contact angle, microregion height, initial microlayer thickness, bubble foot radius of curvature, matching thickness at end of microlayer

    // Fluid properties
    PetscScalar rhol, mul, kappal, Tsat, sigma, hfg, Ri, ls, cpl; // density, viscosity, thermal conductivity, saturation temperature, surface tension, latent heat of vaporization, and interface thermal resistance
    // Constructor

    // Wall temperature (assumed constant for now)
    PetscScalar Tw;

    MicrolayerContext(Foam::IOdictionary &microlayerDict)
    {   
        // Read fluid properties from the dictionary
        rhol = microlayerDict.get<Foam::scalar>("rhol");
        mul = microlayerDict.get<Foam::scalar>("mul");
        kappal = microlayerDict.get<Foam::scalar>("kappal");
        Tsat = microlayerDict.get<Foam::scalar>("Tsat");
        sigma = microlayerDict.get<Foam::scalar>("sigma");
        hfg = microlayerDict.get<Foam::scalar>("hfg");
        Ri = microlayerDict.get<Foam::scalar>("Ri");
        ls = microlayerDict.get<Foam::scalar>("ls");
        cpl = microlayerDict.get<Foam::scalar>("cpl");
        Tw = microlayerDict.get<Foam::scalar>("Tw");

        // Initialize PETSC context with microlayer properties
        xCl = microlayerDict.lookupOrDefault<Foam::scalar>("initialContactLinePosition", 0.0);
        l = microlayerDict.get<Foam::scalar>("microlayerLength");
        if (l > 0)
        {
            isMl = PETSC_TRUE; // Microlayer is present
        }

        dsReq = microlayerDict.lookupOrDefault<Foam::scalar>("desiredGridSize", 1e-6);
        nMin = microlayerDict.lookupOrDefault<Foam::label>("minimumGridPoints", 101); // Minimum number of grid points
        n = Foam::max(int(l / dsReq) + 1, nMin);
        ds = l / (n - 1);

        dt = 0.0; // will be set during iteration
        ub = 0.0; // will be set during iteration

        ucl = microlayerDict.get<Foam::scalar>("contactLineVelocity");
        thetaApp = microlayerDict.get<Foam::scalar>("thetaApp");
        delm = microlayerDict.lookupOrDefault<Foam::scalar>("microregionHeight", 1e-7);
        //del0 = microlayerDict.lookupOrDefault<Foam::scalar>("initialMicrolayerThickness", 1e-6);
        //rc = microlayerDict.lookupOrDefault<Foam::scalar>("bubbleFootRadiusOfCurvature", 1e-6);

        std::cout << "Using Landau-Levich initial microlayer thickness model (Zhang et al. 2024) with Scriven solution for initial condition" << std::endl;

        const Foam::scalar alphal = kappal / (rhol * cpl);

        const Foam::scalar B = microlayerDict.get<Foam::scalar>("growthConstant");

        const Foam::scalar Rb = xCl + l;

        const Foam::scalar Rbdot = 2 * B * B * alphal / Rb;

        const Foam::scalar Rbddot = -4 * Foam::pow(B, 4) * Foam::pow(alphal, 2) / Foam::pow(Rb, 3);

        // Coefficients

        const Foam::scalar a3 = rhol / sigma * (Foam::pow(Rbdot / Rb, 2.0) - Rbddot / (3 * Rb));

        const Foam::scalar a2 = -rhol * Rbddot / (2 * sigma);

        const Foam::scalar a1 = 1 / Rb;

        const Foam::scalar a0 = -1;

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

        const Foam::scalar guess = 1e-6;
        const Foam::scalar min = -100.0;
        const Foam::scalar max = 100.0;
        const int digits = std::numeric_limits<Foam::scalar>::digits;

        Foam::scalar root = boost::math::tools::newton_raphson_iterate(
            f, guess, min, max, digits);

        const Foam::scalar xbar = root;

        const Foam::scalar meniscusCurvature = 3 * a3 * pow(xbar, 2) + 2 * a2 * xbar + a1;

        rc = 1.0 / meniscusCurvature;

        del0 = 1.34 * rc * Foam::pow(mul * Rbdot / sigma, 2.0 / 3.0);

        delM = Foam::min(0.3 * rc, 20 * del0);

        std::cout << "Initial microlayer thickness = " << del0 * 1e6 << " micrometers" << std::endl;


        // Initialize PETSC vectors
        VecCreateSeq(PETSC_COMM_SELF, n, &s);
        VecCreateSeq(PETSC_COMM_SELF, n, &x);
        VecCreateSeq(PETSC_COMM_SELF, n, &h);

        PetscScalar *sp, *xp, *hp;
        VecGetArray(s, &sp);
        VecGetArray(x, &xp);
        VecGetArray(h, &hp);

        std::vector<double> xSp = {0, 0.1 * l, 0.2 * l, 0.22 * l, 0.3 * l, 0.4 * l, 0.5 * l, 0.6 * l, 0.7 * l, 0.8 * l, 1.001 * l};
        std::vector<double> ySp = {0, 1.05 * del0, 0.52 * del0, 0.35 * del0, 0.45 * del0, 0.59 * del0, 0.73 * del0, 0.86 * del0, 0.95 * del0, 1.25 * del0, delM};

        auto spline = boost::math::interpolators::makima<std::vector<double>>(std::move(xSp), std::move(ySp));

        for (PetscInt i = 0; i < n; i++)
        {
            sp[i] = i * ds;
            xp[i] = xCl + i * ds;
            hp[i] = spline(i * ds); // delM * i * ds / l; // Initial microlayer thickness (wedge shape)
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