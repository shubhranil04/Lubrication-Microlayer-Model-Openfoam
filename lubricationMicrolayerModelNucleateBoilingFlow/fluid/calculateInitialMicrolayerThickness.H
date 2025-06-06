const word initialMicrolayerThicknessModel = microlayerDict.lookupOrDefault<word>("initialMicrolayerThicknessModel", "linear");

scalar initialMicrolayerThickness = microlayerDict.lookupOrDefault<scalar>("initialMicrolayerThickness", 0.0);

scalar meniscusSlope = microlayerDict.lookupOrDefault<scalar>("meniscusSlope", 0.0);

scalar meniscusCurvature = microlayerDict.lookupOrDefault<scalar>("meniscusCurvature", 0.0);

if (initialMicrolayerThicknessModel == "linear")
{   
    Info << "Using linear initial microlayer thickness model" << nl;

    const scalar cSlope = microlayerDict.lookupOrDefault<scalar>("cSlope",4.46e-3);  // Default value from Utaka et al. 2013

    Info << "cSlope = " << cSlope << nl;

    initialMicrolayerThickness = cSlope * bubbleRadius1;

    Info << "Initial microlayer thickness = " << initialMicrolayerThickness*1e6 << " micrometers" << nl << "" << nl;

    const scalar sigma(thermophysicalDict.get<dimensionedScalar>("sigma").value());

    const volScalarField& mu1 =
        db.objectRegistry::lookupObject<volScalarField>("thermo:mu.liquid");

    const scalar muLiquid = gSum((mu1.internalField()*mesh.V()).ref())/gSum(mesh.V());

    scalar Rbdot = microlayerDict.lookupOrDefault<scalar>("bubbleExpansionVelocity",0);

    if (Rbdot < 0)
    {
        Rbdot = 0;
    }
    // Calculate meniscus curvature based on the initial microlayer thickness
    meniscusCurvature = 1.34 * Foam::pow(muLiquid * Rbdot / sigma, 2.0/3.0) / 
                         initialMicrolayerThickness;
}

else if (initialMicrolayerThicknessModel == "landauLevich")
{
    Info << "Using Landau-Levich initial microlayer thickness model (Zhang et al. 2024)" << nl;

    const scalar sigma(thermophysicalDict.get<dimensionedScalar>("sigma").value());

    const scalar rhoLiquid = gSum((rho1.internalField()*mesh.V()).ref())/gSum(mesh.V());
    
    const volScalarField& mu1 =
        db.objectRegistry::lookupObject<volScalarField>("thermo:mu.liquid");

    const scalar muLiquid = gSum((mu1.internalField()*mesh.V()).ref())/gSum(mesh.V());

    const scalar Rb = microlayerDict.lookupOrDefault<scalar>("bubbleRadius",0) + 1e-20;

    scalar Rbdot = microlayerDict.lookupOrDefault<scalar>("bubbleExpansionVelocity",0);

    if (Rbdot < 0)
    {
        Rbdot = 0;
    }

    const scalar Rbddot = microlayerDict.lookupOrDefault<scalar>("bubbleExpansionAcceleration",0);

    // Coefficients

    const scalar a3 = rhoLiquid/sigma * (Foam::pow(Rbdot/Rb,2.0) - Rbddot/(3*Rb));

    const scalar a2 = - rhoLiquid * Rbddot / (2* sigma);

    const scalar a1 = 1 / Rb;

    const scalar a0 = -1;

    struct CubicFunction
    {
        Foam::scalar a, b, c, d;

        std::pair<Foam::scalar, Foam::scalar> operator()(Foam::scalar x) const
        {
            Foam::scalar fx = a*x*x*x + b*x*x + c*x + d;
            Foam::scalar dfx = 3*a*x*x + 2*b*x + c;
            return std::make_pair(fx, dfx);
        }
    };

    CubicFunction f{a3, a2, a1, a0}; 

    const scalar guess = 1e-6;
    const scalar min = -100.0;
    const scalar max = 100.0;
    const int digits = std::numeric_limits<Foam::scalar>::digits;

    scalar root = boost::math::tools::newton_raphson_iterate(
        f, guess, min, max, digits
    );

    const scalar xbar = root;

    const scalar temp = a3 * Foam::pow(Rb, 3) + a2 * Foam::pow(Rb, 2) + a1 * Rb + a0;

    if (1-temp*temp > 0){
        meniscusSlope = Foam::pow(temp*temp/(1-temp*temp), 0.5);
    }
    
    else
    {
        meniscusSlope = 0.0;
    }   

    meniscusCurvature = 3*a3*pow(xbar,2) + 2*a2*xbar + a1;

    const scalar meniscusRadius = 1.0 / meniscusCurvature;

    initialMicrolayerThickness = 1.34 * meniscusRadius * Foam::pow(muLiquid * Rbdot / sigma, 2.0/3.0);

    Info << "Initial microlayer thickness = " << initialMicrolayerThickness*1e6 << " micrometers" << nl;

    Info << "Meniscus slope = " << meniscusSlope << nl << "" << nl;
}

microlayerDict.set("meniscusSlope", meniscusCurvature);
microlayerDict.set("meniscusCurvature", meniscusCurvature);
microlayerDict.set("initialMicrolayerThickness", initialMicrolayerThickness);