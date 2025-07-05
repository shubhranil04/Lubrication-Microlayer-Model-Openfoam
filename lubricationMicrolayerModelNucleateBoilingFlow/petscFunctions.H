#ifndef PETSC_FUNCTIONS_H
#define PETSC_FUNCTIONS_H

#include <petscsnes.h>
#include "volFields.H"
#include "fvMesh.H"
#include <boost/math/interpolators/makima.hpp>
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#include <boost/math/tools/roots.hpp>

// Define PETSC Context struct

struct MicrolayerContext
{
    PetscBool isMl; // Flag to indicate if microlayer is present

    PetscScalar r0, umax = 1e3; // Initial bubble radius, maximum bubble velocity
    // contactline position, microlayer length, desired grid size, actual grid size time step
    PetscScalar xCl, l, dsReq, ds, dt, dtMax;
    PetscInt n, nMin; // number of grid points
    Vec s, x, h, q;   // grid in contactline frame, grid in global frame, microlayer thickness, heat flux

    PetscScalar ub = 0, ucl, thetaApp, delm, del0, rc, curv, delM; // bubble expansion velocity, contactline speed, apparent contact angle, microregion height, initial microlayer thickness, bubble foot radius of curvature, matching thickness at end of microlayer

    // Fluid properties
    PetscScalar rhol, mul, kappal, Tsat, sigma, hfg, Ri, ls, cpl, alphal; // density, viscosity, thermal conductivity, saturation temperature, surface tension, latent heat of vaporization, and interface thermal resistance
    // Constructor

    // Wall temperature (assumed constant for now)
    PetscScalar Tw;

    Foam::word initialMicrolayerThicknessModel; // Model for initial microlayer thickness

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

        alphal = kappal / (rhol * cpl); // Thermal diffusivity of liquid

        // Read initial microlayer thickness model
        initialMicrolayerThicknessModel = microlayerDict.lookupOrDefault<Foam::word>("initialMicrolayerThicknessModel", "linear");

        // Initialize PETSC context with microlayer properties
        r0 = microlayerDict.lookupOrDefault<Foam::scalar>("initialBubbleRadius", 10e-6); // Initial bubble radius
        xCl = microlayerDict.lookupOrDefault<Foam::scalar>("initialContactLinePosition", 0.0);
        dtMax = microlayerDict.lookupOrDefault<Foam::scalar>("maximumTimeStep", 1e-8); // Maximum time step
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
        thetaApp = microlayerDict.get<Foam::scalar>("thetaApp")* Foam::constant::mathematical::pi / 180.0; // Convert to radians
        delm = 0; // microlayerDict.lookupOrDefault<Foam::scalar>("microregionHeight", 0);

        const Foam::scalar B = microlayerDict.get<Foam::scalar>("growthConstant");

        // umax = 10 * 2 * B * B * alphal / r0; // Maximum bubble velocity

        ub = 2 * B * B * alphal / r0 * Foam::sin(thetaApp); // Set bubble expansion velocity

        if (initialMicrolayerThicknessModel == "landauLevichScriven" || initialMicrolayerThicknessModel == "landauLevich")
        {
            // If initial microlayer thickness or bubble foot radius of curvature is not set, calculate them

            std::cout << "Using Landau-Levich initial microlayer thickness model (Zhang et al. 2024) with Scriven solution for initial condition" << std::endl;

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

            curv = 3 * a3 * pow(xbar, 2) + 2 * a2 * xbar + a1;

            rc = 1.0 / curv;

            del0 = 1.34 * rc * Foam::pow(mul * Rbdot / sigma, 2.0 / 3.0);

            std::cout << "Initial microlayer thickness = " << del0 * 1e6 << " micrometers" << std::endl;
        }

        else if (initialMicrolayerThicknessModel == "linear")
        {
            std::cout << "Using linear initial microlayer thickness model" << std::endl;

            const Foam::scalar cSlope = microlayerDict.lookupOrDefault<Foam::scalar>("cSlope", 4.46e-3); // Default value from Utaka et al. 2013

            del0 = cSlope * (xCl + l); // Linear model for initial microlayer thickness

            rc = 1e18;

            curv = 0;
        }

        else if (initialMicrolayerThicknessModel == "powerLaw")
        {
            std::cout << "Using power law initial microlayer thickness model" << std::endl;

            const Foam::scalar coeff = microlayerDict.get<Foam::scalar>("coeff");
            const Foam::scalar pow = microlayerDict.get<Foam::scalar>("pow");

            del0 = coeff * Foam::pow(xCl + l, pow);

            std::cout << "Initial microlayer thickness = " << del0 * 1e6 << " micrometers" << std::endl;

            rc = 1e18;

            curv = 0;
        }

        else
        {
            std::cerr << "Fatal error in function: Unknown initialMicrolayerThicknessModel: "
                      << initialMicrolayerThicknessModel << std::endl;
            std::abort();
        }

        delM = del0; // 0.3 * rc; // 10 * del0; // Foam::min(0.3 * rc, 20 * del0);

        // Initialize PETSC vectors
        VecCreateSeq(PETSC_COMM_SELF, n, &s);
        VecCreateSeq(PETSC_COMM_SELF, n, &x);
        VecCreateSeq(PETSC_COMM_SELF, n, &h);

        PetscScalar *sp, *xp, *hp;
        VecGetArray(s, &sp);
        VecGetArray(x, &xp);
        VecGetArray(h, &hp);

        // std::vector<double> xSp = {0, 0.1 * l, 0.2 * l, 0.22 * l, 0.3 * l, 0.4 * l, 0.5 * l, 0.6 * l, 0.7 * l, 0.8 * l, 1.001 * l};
        // std::vector<double> ySp = {0, 1.05 * del0, 0.52 * del0, 0.35 * del0, 0.45 * del0, 0.59 * del0, 0.73 * del0, 0.86 * del0, 0.95 * del0, 1.25 * del0, delM};
        //  auto spline = boost::math::interpolators::makima<std::vector<double>>(std::move(xSp), std::move(ySp));

        // Based on dewetting cfd
        /*
        std::vector<double> xSp = {0, 0.1*l, 0.2*l, 0.3*l, 0.4*l, 0.5*l, 0.6*l, 0.7*l, 0.8*l, 0.9*l, 1*l, 1.1*l};
        //std::vector<double> ySp = {0, 1.05*del0, 0.52*del0, 0.45*del0, 0.59*del0, 0.73*del0, 0.86*del0, 0.95*del0, 1.25*del0, 5.84*del0, 19.8*del0, 48.21*del0};
        std::vector<double> ySp = {0, 0.5*delM, 0.25*delM, 0.2*delM, 0.23*delM, 0.26*delM, 0.29*delM, 0.32*delM, 0.35*delM, 0.5*delM, delM, 1.25*delM};

        double x0Sp = 0.0;
        double dxSp = 0.1 * l; // Uniform spacing

        // Construct the spline
        auto spline = boost::math::interpolators::cardinal_cubic_b_spline<double>(
            ySp.begin(), ySp.end(), x0Sp, dxSp
        );


        for (PetscInt i = 0; i < n; i++)
        {
            sp[i] = i * ds;
            xp[i] = xCl + i * ds;
            hp[i] = spline(i * ds); // delM * i * ds / l; // Initial microlayer thickness (wedge shape)
        }
        */

        // Analytical
        /*
        double a = 2 * del0 / Foam::constant::mathematical::pi;
        double b = Foam::constant::mathematical::pi * std::tan(thetaApp) / (2 * del0);
        double c = 0.5 * curv;
        double d = l - std::sqrt(delM / c);

        int nSp = 21;
        std::vector<double> xSp(nSp);
        std::vector<double> ySp(nSp, 0.0);

        double x0Sp = 0.0;
        double dxSp = l / (nSp - 1); // Uniform spacing

        for (PetscInt i = 0; i < nSp; i++)
        {
            xSp[i] = i * dxSp;
        }

        // Compute ySp based on condition
        for (int i = 0; i < nSp; ++i)
        {
            if (xSp[i] <= d)
            {
                ySp[i] = a * std::atan(b * xSp[i]);
            }
            else
            {
                ySp[i] = c * std::pow((xSp[i] - d), 2) + a * std::atan(b * d);
            }
        }

        // Construct the spline
        auto spline = boost::math::interpolators::cardinal_cubic_b_spline<double>(
            ySp.begin(), ySp.end(), x0Sp, dxSp
        );


        for (PetscInt i = 0; i < n; i++)
        {
            sp[i] = i * ds;
            xp[i] = xCl + i * ds;
            hp[i] = spline(i * ds); // delM * i * ds / l; // Initial microlayer thickness (wedge shape)
        }
        */

        // taninv

        struct CoeffFunction
        {
            Foam::scalar theta_, delta_, lml_;

            std::pair<Foam::scalar, Foam::scalar> operator()(Foam::scalar x) const
            {
                Foam::scalar A = Foam::tan(theta_) / delta_;
                Foam::scalar u = Foam::atan(x * lml_);
                Foam::scalar du_dx = lml_ / (1 + Foam::pow(x * lml_, 2));

                Foam::scalar fx = A * u / x - 1;
                Foam::scalar dfx = A * (x * du_dx - u) / (x * x);
                return std::make_pair(fx, dfx);
            }
        };

        CoeffFunction f{thetaApp, del0, l};

        const Foam::scalar guess = 1e6;
        const Foam::scalar min = 0;
        const Foam::scalar max = 1e10;
        const int digits = std::numeric_limits<Foam::scalar>::digits;

        Foam::scalar b = boost::math::tools::newton_raphson_iterate(
            f, guess, min, max, digits);

        Foam::scalar a = Foam::tan(thetaApp) / b;

        for (PetscInt i = 0; i < n; i++)
        {
            sp[i] = i * ds;
            xp[i] = xCl + i * ds;

            hp[i] = a * Foam::atan(b * i * ds);
        }

        VecRestoreArray(s, &sp);
        VecRestoreArray(x, &xp);
        VecRestoreArray(h, &hp);

        // Write initial conditions
        PetscViewer viewer;

        Foam::fileName xml_file = "x0.csv";

        PetscViewerASCIIOpen(PETSC_COMM_WORLD, xml_file.c_str(), &viewer);
        VecView(x, viewer);
        PetscViewerDestroy(&viewer);

        Foam::fileName hml_file = "h0.csv";

        PetscViewerASCIIOpen(PETSC_COMM_WORLD, hml_file.c_str(), &viewer);
        VecView(h, viewer);
        PetscViewerDestroy(&viewer);
    }
};

PetscErrorCode FormFunction(SNES snes, Vec h, Vec f, void *ctx_);

PetscErrorCode FormJacobian(SNES snes, Vec h, Mat jac, Mat B, void *ctx_);

PetscErrorCode computeHeatFlux(void *ctx_);

PetscScalar computeMeanHeatFlux(PetscScalar xMin, PetscScalar xMax, void *ctx_);

#endif

// ************************************************************************* //