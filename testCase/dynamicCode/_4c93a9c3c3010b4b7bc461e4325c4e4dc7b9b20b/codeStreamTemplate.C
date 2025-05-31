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
#line 26 "/home/shubhranil/OpenFOAM/AddOns-v2106/TwoPhaseFlow/solver/testCase/0/fluid/T.liquid.#codeStream"
#include "fvCFD.H"
//}}} end codeInclude

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C" void codeStream_4c93a9c3c3010b4b7bc461e4325c4e4dc7b9b20b(Foam::Ostream& os, const Foam::dictionary& dict)
{
//{{{ begin code
    #line 43 "/home/shubhranil/OpenFOAM/AddOns-v2106/TwoPhaseFlow/solver/testCase/0/fluid/T.liquid.#codeStream"
/* Access to internal mesh information */
    const IOdictionary& d = static_cast<const IOdictionary&>(dict);
    const fvMesh& mesh = refCast<const fvMesh>(d.db());

    scalarField T(mesh.nCells(), 0.); /* initialization */

    /* User variables initialization */
 
    const scalar Tact = 379.17;
    const scalar Tbulk = 373.15;
    const scalar blThickness = 0.5e-3;	
    
    forAll(T, i) /* Loop on elements */
    {
      const scalar y = mesh.C()[i][1]; /* Access to y component of cell centers coordinates */
      
      if (y < blThickness)
      {
      	T[i] = Tact - y/blThickness*(Tact - Tbulk);
      }
      
      else
      {
      	T[i] = Tbulk;
      }	
    }

    T.writeEntry("", os); /* Write output in the dictionary */
//}}} end code
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

