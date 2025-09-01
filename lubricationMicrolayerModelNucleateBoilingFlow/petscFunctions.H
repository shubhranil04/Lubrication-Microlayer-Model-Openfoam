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

    PetscBool ridgeDepleted = PETSC_FALSE;
    PetscScalar dewettingFlag = 1;
    PetscScalar wedgeAngle; // Wedge angle of mesh (in degrees)
    PetscScalar r0, umax = 1e3; // Initial bubble radius, maximum bubble velocity
    PetscScalar xCl, xB, l, dsReq, ds, dt = 0, dtMax; // contactline position, microlayer end position, microlayer length, desired grid size, actual grid size, time step
    PetscInt n; // number of grid points
    Vec s, x, h, q, Twall;   // grid in contactline frame, grid in global frame, microlayer thickness, heat flux, wall temperature

    PetscInt varCa; // if true, calculate contact line velocity during runtime
    PetscScalar A_theta, b_theta; // powerlaw coeffs for evaporative contact angle $ \theta_{ev} = \theta_0 + A \Delta T ^ b $ (in radians)
    PetscScalar ub = 0, ab = 0, ucl, TCl, theta0, thetaApp, delm, del0, rc, curv, delM; // bubble expansion velocity, contactline speed, apparent contact angle, microregion height, initial microlayer thickness, bubble foot radius of curvature, matching thickness at end of microlayer
    PetscScalar rbCalc;
    // Fluid properties
    PetscScalar rhol, mul, kappal, Tsat, sigma, hfg, Ri, ls, cpl, alphal; // density, viscosity, thermal conductivity, saturation temperature, surface tension, latent heat of vaporization, and interface thermal resistance
    // Constructor

    // Wall temperature (assumed constant for now)
    PetscScalar Tw;

    Foam::word initialMicrolayerThicknessModel; // Model for initial microlayer thickness

    MicrolayerContext(Foam::IOdictionary &microlayerDict, const Foam::Time &runTime)
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
	TCl = Tw;
        alphal = kappal / (rhol * cpl); // Thermal diffusivity of liquid
	dewettingFlag = microlayerDict.lookupOrDefault<Foam::scalar>("dewettingFlag", 1);
	wedgeAngle = microlayerDict.get<Foam::scalar>("wedgeAngle");
        // Read initial microlayer thickness model
        initialMicrolayerThicknessModel = microlayerDict.lookupOrDefault<Foam::word>("initialMicrolayerThicknessModel", "linear");
	
	// Read flag for variable contact line Ca
	varCa = microlayerDict.lookupOrDefault<Foam::label>("variableCapillaryNum",0);

        // Initialize PETSC context with microlayer properties
        r0 = microlayerDict.lookupOrDefault<Foam::scalar>("bubbleRadius", 10e-6); // Initial bubble radius
        rbCalc = r0;
		
	xCl = microlayerDict.get<Foam::scalar>("contactLinePosition");
        dtMax = microlayerDict.lookupOrDefault<Foam::scalar>("maximumTimeStep", 1e-8); // Maximum time step
        l = microlayerDict.get<Foam::scalar>("microlayerLength");
        if (l > 0)
        {
            isMl = PETSC_TRUE; // Microlayer is present
        }
	xB = xCl + l;

        dsReq = microlayerDict.lookupOrDefault<Foam::scalar>("desiredGridSize", 1e-6);

        ucl = microlayerDict.get<Foam::scalar>("contactLineVelocity");
        theta0 = microlayerDict.get<Foam::scalar>("theta0") * Foam::constant::mathematical::pi / 180.0; // Convert to radians
        delm = microlayerDict.get<Foam::scalar>("microregionHeight");                                                                                           
	A_theta = microlayerDict.get<Foam::scalar>("A_theta");
	b_theta = microlayerDict.get<Foam::scalar>("b_theta");
	
	thetaApp = theta0 + A_theta * Foam::pow(TCl-Tsat, b_theta);
        const Foam::scalar B = microlayerDict.get<Foam::scalar>("growthConstant");

        ub = 2 * B * B * alphal / r0; // Set bubble expansion velocity using Scriven solution

	ab = - 4 * Foam::pow(B,4) * Foam::pow(alphal,2) / Foam::pow(r0,3); // Set bubble expansion acceleration using Scriven solution

        if (runTime.value() == 0)
        {

            n = int(l / dsReq) + 1;
            ds = l / (n - 1);

            // Initialize PETSC vectors
            VecCreateSeq(PETSC_COMM_SELF, n, &s);
            VecCreateSeq(PETSC_COMM_SELF, n, &x);
            VecCreateSeq(PETSC_COMM_SELF, n, &h);

            PetscScalar *sp, *xp, *hp;
            VecGetArray(s, &sp);
            VecGetArray(x, &xp);
            VecGetArray(h, &hp);

            if (initialMicrolayerThicknessModel == "landauLevichScriven" || initialMicrolayerThicknessModel == "landauLevich")
            {
                // If initial microlayer thickness or bubble foot radius of curvature is not set, calculate them

                std::cout << "Using Landau-Levich initial microlayer thickness model (Zhang et al. 2024) with Scriven solution for initial condition" << std::endl;

                const Foam::scalar Rb = r0;

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

            delM = del0; 

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

                CoeffFunction f_{thetaApp, del0 - delm, l};

                const Foam::scalar guess_ = 1e6;
                const Foam::scalar min_ = 0;
                const Foam::scalar max_ = 1e10;
                const int digits_ = std::numeric_limits<Foam::scalar>::digits;

                Foam::scalar b = boost::math::tools::newton_raphson_iterate(
                    f_, guess_, min_, max_, digits_);

                Foam::scalar a = Foam::tan(thetaApp) / b;

                for (PetscInt i = 0; i < n; i++)
                {
                    sp[i] = i * ds;
                    xp[i] = xCl + i * ds;

                    hp[i] = delm + a * Foam::atan(b * i * ds);
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

        else // For restart of simulation from a previous time step
        {
            if (l > 0)
            {
                std::cout << "Reading microlayer information from previous time step ..." << std::endl;

                Foam::fileName timeDir = runTime.timePath();
                std::string xmlPath = (timeDir / "xml.dat").c_str();
                std::string hmlPath = (timeDir / "hml.dat").c_str();

                PetscViewer viewer;

                PetscViewerBinaryOpen(PETSC_COMM_WORLD, xmlPath.c_str(), FILE_MODE_READ, &viewer);
                VecCreateSeq(PETSC_COMM_SELF, PETSC_DECIDE, &x);
                VecLoad(x, viewer);
                PetscViewerDestroy(&viewer);

                VecGetSize(x, &n);
                ds = l / (n - 1);

                PetscViewerBinaryOpen(PETSC_COMM_WORLD, hmlPath.c_str(), FILE_MODE_READ, &viewer);
                VecCreateSeq(PETSC_COMM_SELF, PETSC_DECIDE, &h);
                VecLoad(h, viewer);
                PetscViewerDestroy(&viewer);

                VecCreateSeq(PETSC_COMM_SELF, n, &s);
                VecCopy(x, s);
                VecShift(s, -xCl);

                del0 = microlayerDict.get<Foam::scalar>("initialMicrolayerThickness");
                rc = microlayerDict.get<Foam::scalar>("bubbleFootRadiusOfCurvature");
                curv = microlayerDict.get<Foam::scalar>("meniscusCurvature");
                delM = del0;
            }
        }
    } // End of constructor
};

PetscErrorCode FormFunction(SNES snes, Vec h, Vec f, void *ctx_);

PetscErrorCode FormJacobian(SNES snes, Vec h, Mat jac, Mat B, void *ctx_);

PetscErrorCode boundThickness(void *ctx_, PetscScalar lim_);

PetscErrorCode getWallTemperature(void *ctx_, const Foam::fvPatchScalarField& Twall, const Foam::fvMesh& mesh);

PetscErrorCode computeHeatFlux(void *ctx_);

PetscScalar computeMeanHeatFlux(PetscScalar xMin, PetscScalar xMax, void *ctx_);

#endif

// ************************************************************************* //VecShift(s, -xCl);
