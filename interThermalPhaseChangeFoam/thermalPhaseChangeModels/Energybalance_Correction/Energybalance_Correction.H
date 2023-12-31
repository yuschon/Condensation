/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 Alex Rattner and Mahdi Nabil
     \\/     M anipulation  |
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

Class
    Foam::thermalPhaseChangeModel::Energybalance_Correction

Description
    This model implements the corrected energy balance model


SourceFiles
    Energybalance_Correction.C

\*--------------------------------------------------------------------*/

#ifndef Energybalance_Correction_H
#define Energybalance_Correction_H

#include "thermalPhaseChangeModel.H"
#include "MeshGraph.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalPhaseChangeModels
{

/*--------------------------------------------------------------------*\
                              Class EmpiricalRateParameter
\*--------------------------------------------------------------------*/

class Energybalance_Correction
:
    public thermalPhaseChangeModel
{

// Private Member Functions

    //- Calculate and return the phase change heat rate
    void calcQ_pc();

    //- Phase change heat
    volScalarField Q_pc_;

    //- Calculate the heat source coefficient
    dimensionedScalar me;

    //- Calculate the correction coefficient
    dimensionedScalar mcor;


public:

    //- Runtime type information
    TypeName("Energybalance_Correction");


    // Constructors

    //- construct from components
    Energybalance_Correction
    (
        const word& name,
        const dictionary& thermalPhaseChangeProperties,
        const twoPhaseThermalMixture& twoPhaseProperties,
        const volScalarField& T,
        const volScalarField& alpha1
    );


    //- Destructor
    ~Energybalance_Correction()
    {}


    // Member Functions

        //- Return phase change heat
        tmp<volScalarField> Q_pc() const
        {
            return Q_pc_;
        }

        //- Correct the thermalPhaseChange model
        void correct()
        {
            calcQ_pc();
        }

        //- Read the transportProperties dictionary and update
        bool read(const dictionary& thermalPhaseChangeProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermalPhaseChangeModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
