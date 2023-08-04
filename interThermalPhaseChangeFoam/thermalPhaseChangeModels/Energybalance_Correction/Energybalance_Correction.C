/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022 schon
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

\*---------------------------------------------------------------------------*/

#include "Energybalance_Correction.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalPhaseChangeModels
{
    defineTypeNameAndDebug(Energybalance_Correction, 0);
    addToRunTimeSelectionTable
    (
        thermalPhaseChangeModel,
        Energybalance_Correction,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalPhaseChangeModels::Energybalance_Correction::Energybalance_Correction
(
        const word& name,
        const dictionary& thermalPhaseChangeProperties,
        const twoPhaseThermalMixture& twoPhaseProperties,
        const volScalarField& T,
        const volScalarField& alpha1
)
:
    thermalPhaseChangeModel
    (
        name,
        thermalPhaseChangeProperties,
        twoPhaseProperties,
        T,
        alpha1
    ),
    Q_pc_
    (
        IOobject
        (
            "PhaseChangeHeat",
            T_.time().timeName(),
            T.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        T.mesh(),
        dimensionedScalar( "dummy", dimensionSet(1,-1,-3,0,0,0,0), 0 )
    )
{

    //reading me and mcor
    thermalPhaseChangeProperties_.lookup("me") >> me;   
    thermalPhaseChangeProperties_.lookup("mcor") >> mcor; 


    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::thermalPhaseChangeModels::Energybalance_Correction::calcQ_pc()
{
    const dimensionedScalar& kl = twoPhaseProperties_.lambda1();
    volVectorField gradT = fvc::grad(T_);
    volVectorField gradalpha1 = fvc::grad(alpha1_);
    Q_pc_ = me*pos(1-alpha1_)*(gradT&gradalpha1)*kl*(1+mcor*(T_sat_-T_)/T_sat_);

}


bool Foam::thermalPhaseChangeModels::Energybalance_Correction::
read(const dictionary& thermalPhaseChangeProperties)
{
    thermalPhaseChangeModel::read(thermalPhaseChangeProperties);
    thermalPhaseChangeProperties_.lookup("me") >> me;
    thermalPhaseChangeProperties_.lookup("mcor") >> mcor;
    return true;
}


// ************************************************************************* //
