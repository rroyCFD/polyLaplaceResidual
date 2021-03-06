/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is compatible with OpenFOAM.

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

Class
    Foam::polyLaplaceResidual


Description 
    Polynomial Laplace residual for LES (Return residual of the field). The 
    polynomial coefficients are derived follwing the publication by Trias et al.

    Trias, F. X., & Verstappen, R. W. C. P. (2011). On the construction of
    discrete filters for symmetry-preserving regularization models. Computers
    and Fluids, 40(1), 139–148. https://doi.org/10.1016/j.compfluid.2010.08.015


SourceFiles
    polyLaplaceResidual.C

\*---------------------------------------------------------------------------*/

#ifndef polyLaplaceResidual_H
#define polyLaplaceResidual_H

#include "LESfilter.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyLaplaceResidual Declaration
\*---------------------------------------------------------------------------*/

class polyLaplaceResidual : public LESfilter
{
    // Private data

        // polynomial coefficients
        scalar d1_;
        scalar d2_;

        surfaceScalarField deltaSquared_;

    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        polyLaplaceResidual(const polyLaplaceResidual&);
        void operator=(const polyLaplaceResidual&);


public:

    //- Runtime type information
    TypeName("polyLaplaceResidual");


    // Constructors

        //- Construct from components
        polyLaplaceResidual(const fvMesh& mesh, scalar d1, scalar d2);

        //- Construct from IOdictionary
        polyLaplaceResidual(const fvMesh& mesh, const dictionary&);


    //- Destructor
    virtual ~polyLaplaceResidual()
    {}


    // Member Functions

        //- Read the LESfilter dictionary
        virtual void read(const dictionary&);


    // Member Operators

        virtual tmp<volScalarField> operator()
        (
            const tmp<volScalarField>&
        ) const;

        virtual tmp<volVectorField> operator()
        (
            const tmp<volVectorField>&
        ) const;

        virtual tmp<volSymmTensorField> operator()
        (
            const tmp<volSymmTensorField>&
        ) const;

        virtual tmp<volTensorField> operator()
        (
            const tmp<volTensorField>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
