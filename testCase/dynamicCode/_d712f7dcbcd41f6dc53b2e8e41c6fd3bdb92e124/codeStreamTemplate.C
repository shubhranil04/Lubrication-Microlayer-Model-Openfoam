/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) YEAR AUTHOR, AFFILIATION
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description
    Template for use with codeStream.

\*---------------------------------------------------------------------------*/

#include "dictionary.H"
#include "Ostream.H"
#include "Pstream.H"
#include "pointField.H"
#include "tensor.H"
#include "unitConversion.H"

//{{{ begin codeInclude
#line 25 "/home/shubhranil/OpenFOAM/AddOns-v2106/TwoPhaseFlow/solver/testCase/0/solid/T.#codeStream"
#include "fvCFD.H"
//}}} end codeInclude

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C" void codeStream_d712f7dcbcd41f6dc53b2e8e41c6fd3bdb92e124(Foam::Ostream& os, const Foam::dictionary& dict)
{
//{{{ begin code
    #line 42 "/home/shubhranil/OpenFOAM/AddOns-v2106/TwoPhaseFlow/solver/testCase/0/solid/T.#codeStream"
/* Access to internal mesh information */
    const IOdictionary& d = static_cast<const IOdictionary&>(dict);
    const fvMesh& mesh = refCast<const fvMesh>(d.db());

    scalarField T(mesh.nCells(), 0.); /* initialization */

    /* User variables initialization */
    const scalar kappa = 401;
    const scalar Tact = 379.17;
    const scalar q = 34e3;
    
    forAll(T, i) /* Loop on elements */
    {
      const scalar y = mesh.C()[i][1]; /* Access to y component of cell centers coordinates */
      T[i] = Tact-q/kappa*y; /* T assignment for all elements */
    }

    T.writeEntry("", os); /* Write output in the dictionary */
//}}} end code
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

