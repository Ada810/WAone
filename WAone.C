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

\*---------------------------------------------------------------------------*/

#include "WAone.H"
#include "bound.H"
#include "wallDist.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
template<class BasicTurbulenceModel>
tmp<volScalarField> WAone<BasicTurbulenceModel>::chi() const
{
    return Rnu_/this->nu();
    
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WAone<BasicTurbulenceModel>::fmu
(
    const volScalarField& chi
) const
{
    const volScalarField chi3(pow3(chi));
    return chi3/(chi3 + pow3(Cw_));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WAone<BasicTurbulenceModel>::WDF_R
(
    const volScalarField& S,
    const volScalarField& W   
) const
{
    return mag(W/S);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WAone<BasicTurbulenceModel>::WDF_omega
(
    const volScalarField& S 
) const
{
    return S/sqrt(Cmu_);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WAone<BasicTurbulenceModel>::WDF_k
(
    const volScalarField& omega
) const
{
    return this->nut_*omega;
	//return Rnu_*omega;
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WAone<BasicTurbulenceModel>::arg1
(
    const volScalarField& S,
    const volScalarField& W
) const
{
    const volScalarField R = WDF_R(S, W);
    const volScalarField omega = WDF_omega(S);
    const volScalarField k = WDF_k(omega);

    const volScalarField eta = S*max(1.0, R);

    return (this->nu()+Rnu_)/2 * sqr(eta)/max(Cmu_*k*omega,
                                              dimensionedScalar("SMALL", 
                                                                dimensionSet(0, 2, -3, 0, 0), 
                                                                SMALL)
                                             );
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WAone<BasicTurbulenceModel>::calcSwitch
(
    const volScalarField& S,
    const volScalarField& W    
) const
{
    return tanh(pow(Cs1_*arg1(S, W), Cs2_));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WAone<BasicTurbulenceModel>::blend
(
    const volScalarField& Switch,
    const dimensionedScalar& psi1,
    const dimensionedScalar& psi2
) const
{
    return Switch*(psi1 - psi2) + psi2;
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WAone<BasicTurbulenceModel>::sigmaR(const volScalarField& Switch) const
{
    return blend(Switch, sigmakw_, sigmake_);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WAone<BasicTurbulenceModel>::C1(const volScalarField& Switch) const
{
    return blend(Switch, C1kw_, C1ke_);
}

///////////////////////////////////////////////////////////////////transition
template<class BasicTurbulenceModel>
tmp<volScalarField> WAone<BasicTurbulenceModel>::Rev
(
  const volScalarField& W
) const
{
   return this->rho_*sqr(y_)*W/this->mu();
   //return this->rho_*sqr(y_)*S/this->mu();
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WAone<BasicTurbulenceModel>::Retheta
(
  const volScalarField& Rev
) const
{
   return Rev/2.193;
}



/////////////////////////////////////////////////////////////////////

template<class BasicTurbulenceModel>
void WAone<BasicTurbulenceModel>::correctNut
(
    const volScalarField& fmu
)
{
    this->nut_ = Rnu_*fmu;
    this->nut_.correctBoundaryConditions();

    BasicTurbulenceModel::correctNut();
}

template<class BasicTurbulenceModel>
void WAone<BasicTurbulenceModel>::correctNut()
{
    correctNut(fmu(this->chi()));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
WAone<BasicTurbulenceModel>::WAone
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel> >
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),
    Cw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw",
            this->coeffDict_,
            8.54
        )
    ),

    C1ke_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1ke",
            this->coeffDict_,
			0.1284
        )
    ),

    C1kw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1kw",
            this->coeffDict_,
			0.0829
        )
    ),

    sigmake_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmake",
            this->coeffDict_,
			1.0
        )
    ),

    sigmakw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmakw",
            this->coeffDict_,
			0.72
        )
    ),

    C2ke_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2ke",
            this->coeffDict_,
			C1ke_.value()/sqr(kappa_.value()) + sigmake_.value()
        )
    ),

    C2kw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2kw",
            this->coeffDict_,
			C1kw_.value()/sqr(kappa_.value()) + sigmakw_.value()
        )
    ),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
			0.09
        )
    ),

    Cs1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cs1",
            this->coeffDict_,
			1.0
        )
    ),

    Cs2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cs2",
            this->coeffDict_,
			4.0
        )
    ),

    Cm_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cm",
            this->coeffDict_,
			8.0
        )
    ),

    Rnu_
    (
        IOobject
        (
            "Rnu",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    
    utau_
    (
        IOobject
        (
            "utau",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("0.0", dimensionSet(0, 1, -1, 0, 0), 0.0)
    ),
    
    wGUx_
    (
        IOobject
        (
            "wGUx",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("0.0", dimensionSet(0, 0, -1, 0, 0), 0.0)
    ),
    
    wallGradUTest_
    (
        IOobject
        (
            "wallGradUTest",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedVector
        (
            "wallGradU",
            U.dimensions()/dimLength,
            vector::zero
        )
    ),

    f1_
    (
        IOobject
        (
            "f1",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("0.0", dimensionSet(0, 0, 0, 0, 0), 0.0)
    ),

    S_
    (
        IOobject
        (
            "S",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("0.0", dimensionSet(0, 0, -1, 0, 0), 0.0)
    ),

    W_
    (
        IOobject
        (
            "W",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("0.0", dimensionSet(0, 0, -1, 0, 0), 0.0)
    ),

///////////////////transition
    chi1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "chi1",
            this->coeffDict_,
            0.02
        )
    ),
    /*
    chi2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "chi2",
            this->coeffDict_,
            50.0
        )
    ),
    
*/
    Tu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Tu",
            this->coeffDict_,
            1.0
        )
    ),
    
    RethetaC_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "RethetaC",
            this->coeffDict_,
            //803.73*pow((Tu_.value()+0.6067),-1.027)
            803.73*pow((Tu_.value()+0.6067),-1.027)
        )
    ),
    
    gamma_
    (
        IOobject
        (
            "gamma",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    this->mesh_,
    dimensionedScalar("gamma", dimless, scalar(1.0))
    ),

    term1m_
    (
        IOobject
        (
            "term1m",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    this->mesh_,
    dimensionedScalar("1", dimless, scalar(0.0))
    ),

    term2m_
    (
        IOobject
        (
            "term2m",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    this->mesh_,
    dimensionedScalar("2", dimless, scalar(0.0))
    ),

    nuBCm_
    (
        IOobject
        (
            "nuBCm",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    this->mesh_,
    dimensionedScalar("0.0", dimensionSet(0, 0, 0, 0, 0), 1.0)
    ),
    
    chi2m_
    (
        IOobject
        (
            "chi2m",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    this->mesh_,
    dimensionedScalar("0.0", dimensionSet(0, 0, 0, 0, 0), 1.0)
    ),
    
    Rethetam_
    (
        IOobject
        (
            "Rethetam",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    this->mesh_,
    dimensionedScalar("0.0", dimensionSet(0, 0, 0, 0, 0), 1.0)
    ),
    
    tm_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "tm",
            this->coeffDict_,
            50
        )
    ),
    
    plim_
    (
        IOobject
        (
            "plim",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("0.0", dimensionSet(0, 2, -2, 0, 0), 0.0)
    ),
    
    Ct1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ct1",
            this->coeffDict_,
            1.0
        )
    ),
    
    CP1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CP1",
            this->coeffDict_,
            1.2
        )
    ), 
    
    CP2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CP2",
            this->coeffDict_,
            1.0
        )
    ), 

	y_(wallDist::New(this->mesh_).y())

{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool WAone<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel> >::read())
    {   
        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
tmp<volScalarField> WAone<BasicTurbulenceModel>::DRnuEff(volScalarField Switch) const
{
    return tmp<volScalarField>
    (
        new volScalarField("DRnuEff", Rnu_*sigmaR(Switch) + this->nu())
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> WAone<BasicTurbulenceModel>::k() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "k",
                this->runTime_.timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimensionedScalar("0", dimensionSet(0, 2, -2, 0, 0), 0)
        )
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> WAone<BasicTurbulenceModel>::epsilon() const
{
    WarningInFunction
        << "Turbulence kinetic energy dissipation rate not defined for "
        << "Spalart-Allmaras model. Returning zero field"
        << nl;

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "epsilon",
                this->runTime_.timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimensionedScalar("0", dimensionSet(0, 2, -3, 0, 0), 0)
        )
    );
}


template<class BasicTurbulenceModel>
void WAone<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
	const volVectorField& U = this->U_;


    eddyViscosity<RASModel<BasicTurbulenceModel> >::correct();


    // Calculate strain rate magnitude S
	volScalarField S2(2.0*magSqr(symm(fvc::grad(this->U_))));
	volScalarField S = sqrt(S2);
	bound(S, dimensionedScalar("0", S.dimensions(), SMALL)); // SMALL = 1e-15
	bound(S2, dimensionedScalar("0", S2.dimensions(), SMALL));
	S_ = S;

    // Calculate vorticity magnitude W
    volScalarField W2(2.0*magSqr(skew(fvc::grad(this->U_))));
    volScalarField W = sqrt(W2);
	bound(W, dimensionedScalar("0", W.dimensions(), SMALL));
	bound(W2, dimensionedScalar("0", W2.dimensions(), SMALL));
	W_ = W;

	volScalarField magU = mag(U);
	volScalarField nueff = this->nut_+this->nu();
	const volScalarField Rev(this->Rev(W));
    const volScalarField Retheta(this->Retheta(Rev));

	//volScalarField nuBC = this->nut_/(mag(this->U_)*y_);
	

	volScalarField nuBC = Rnu_/(mag(this->U_)*y_);
	//nuBCm_ = nuBC;
	volScalarField Term2 = nuBC;

	forAll(nuBC, cellI)
{
	//nuBC[cellI] = (0.3*Rnu_[cellI])/(S[cellI]*sqr(y_[cellI]));
	nuBC[cellI] = (CP2_.value()*Rnu_[cellI])/(S[cellI]*sqr(y_[cellI]));
	//Term2[cellI] = max(nuBC[cellI]-0.0005,0.0)/0.0005;
}
	nuBCm_=nuBC;
	term2m_= max(tm_*this->nut_/this->nu(),0.0);
	
	term1m_ = max(CP1_*Retheta-RethetaC_,0.0)/(chi1_*RethetaC_); //2.5
	
	gamma_ = 1.0-exp(-sqrt(Ct1_*term1m_)-sqrt(term2m_)); //0.01term1
	
	
    gamma_=min(gamma_,scalar(1.0));
    bound(gamma_,scalar(0));
    eddyViscosity<RASModel<BasicTurbulenceModel> >::correct();
  //////////transitionBC local reference
   

    
   
	
	
	chi2m_ = this->nut_/this->nu();
	Rethetam_ = Retheta;
	

    // Calculate switch function (f1)
    f1_ = calcSwitch(S, W);
    plim_ = 0.5*max(gamma_-0.2, 0.0)*(1.0-gamma_)*min(max((Rev/2420.0)-1.0, 0.0), 3.0)*max(3.0*this->nu()-this->nut_, dimensionedScalar("0", dimensionSet(0, 2, -1, 0, 0), 0))*W;

    eddyViscosity<RASModel<BasicTurbulenceModel> >::correct();

    // Define and solve R-Equation
    tmp<fvScalarMatrix> RnuEqn
    (
        fvm::ddt(alpha, rho, Rnu_)
      + fvm::div(alphaRhoPhi, Rnu_)
      - fvm::laplacian(alpha*rho*DRnuEff(f1_), Rnu_)
     ==
        alpha*rho*C1(f1_)*gamma_*fvm::Sp(S, Rnu_)
      + alpha*rho*f1_*C2kw_*fvm::Sp((fvc::grad(Rnu_)&fvc::grad(S))/S, Rnu_)
      - alpha*rho*(1.0-f1_)*min(C2ke_*Rnu_*Rnu_*magSqr(fvc::grad(S))/S2,
                                Cm_*magSqr(fvc::grad(Rnu_)))
    
                              
                             
    );


    
    
    RnuEqn().relax();
    solve(RnuEqn);
    bound(Rnu_, dimensionedScalar("0", Rnu_.dimensions(), 0.0));
    Rnu_.correctBoundaryConditions();
    correctNut();
    
   

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
