/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::RASModels::WAone

Group
    grpRASTurbulence

Description
	Wray-Agarwal Turbulence Model with Rahman's wall distance free modification
with type-2 switch function, and Han's coefficient adjustment. 

SourceFiles
    WAone.C

\*---------------------------------------------------------------------------*/

#ifndef WAone_H
#define WAone_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                      Class WAone Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class WAone
:
    public eddyViscosity<RASModel<BasicTurbulenceModel> >
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        WAone(const WAone&);
        WAone& operator=(const WAone&);

protected:

    // Protected data

        // Model coefficients

            dimensionedScalar kappa_;
            dimensionedScalar ks_;
            dimensionedScalar Cw_;
            dimensionedScalar Cr1_;
            dimensionedScalar Cr2_;
            dimensionedScalar Cr3_;
            dimensionedScalar C1ke_;
			dimensionedScalar C1kw_;
            dimensionedScalar sigmake_;
            dimensionedScalar sigmakw_;
			dimensionedScalar C2ke_;
            dimensionedScalar C2kw_;
            dimensionedScalar Cmu_;
            dimensionedScalar Cs1_;
            dimensionedScalar Cs2_;
            dimensionedScalar Cm_;

        
        // Fields

            volScalarField Rnu_;
            volScalarField utau_;
            volScalarField wGUx_;
            volVectorField wallGradUTest_;
			volScalarField f1_;
            volScalarField S_;
            volScalarField W_;
///tran
            dimensionedScalar chi1_;
            //dimensionedScalar chi2_;
            dimensionedScalar Tu_;
            dimensionedScalar RethetaC_;
			volScalarField gamma_;
			volScalarField term1m_;
			volScalarField term2m_;
			volScalarField nuBCm_;
			volScalarField chi2m_;
			volScalarField Rethetam_;
            dimensionedScalar tm_;
            volScalarField plim_;
            dimensionedScalar Ct1_;
            dimensionedScalar CP1_;
            dimensionedScalar CP2_;
            
            
            
			const volScalarField& y_;

    // Protected Member Functions

        tmp<volScalarField> chi() const;

        tmp<volScalarField> fmu(const volScalarField& chi) const;

        tmp<volScalarField> WDF_R
        (
            const volScalarField& S,
            const volScalarField& W    
        ) const;

        tmp<volScalarField> WDF_omega
        (
            const volScalarField& S 
        ) const;

        tmp<volScalarField> WDF_k
        (
            const volScalarField& omega 
        ) const;

        tmp<volScalarField> arg1
        (
            const volScalarField& S,
            const volScalarField& W    
        ) const;

        tmp<volScalarField> calcSwitch
        (
            const volScalarField& S,
            const volScalarField& W    
        ) const;

        tmp<volScalarField> blend
        (
            const volScalarField& Switch,
            const dimensionedScalar& psi1,
            const dimensionedScalar& psi2
        ) const;

        tmp<volScalarField> sigmaR(const volScalarField& Switch) const;

        tmp<volScalarField> C1(const volScalarField& Switch) const;

        tmp<volScalarField> Rev(const volScalarField& W) const;
        
        tmp<volScalarField> Retheta(const volScalarField& Rev) const;
        
     

        void correctNut(const volScalarField& fmu);
        virtual void correctNut();

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("WAone");

    // Constructors

        //- Construct from components
        WAone
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~WAone()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for Rnu
        tmp<volScalarField> DRnuEff(volScalarField Switch) const;

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const;

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "WAone.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
