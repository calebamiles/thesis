/*---------------------------------------------------------------------------*\
=========                 |
\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
 \\    /   O peration     |
  \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

Application
	icoFoam

Description
	Transient solver for incompressible, laminar flow of Newtonian fluids.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IOEquationReader.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
	#include "setRootCase.H"
	#include "createTime.H"
	#include "createMesh.H"
	#include "initContinuityErrs.H"
	#include "createFields.H"
	#include "readTimeControls.H"
	#include "CourantNo.H"
	#include "setInitialDeltaT.H"
	#include "createEquationReader.H"

	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	fileName path(args.rootPath()/args.caseName());	
	//equationReader eqns; //commented for update

	//Create time variable "t" for use in symbolic expressions with equationReader
	eqns.scalarSources().addSource(runTime.value(), "t", dimTime);

   //Create position variables x, y, and z for use in symbolic expressions with equationReader  
	eqns.vectorSources().addSource(mesh.C());

  
  /*
  
	//Define scalar fields for validation        
	volScalarField bOneX(
		IOobject
		("bOneX",
		runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
		mesh,
		dimAcceleration
		);
		
	volScalarField bOneY(
		IOobject
		("bOneY",
		runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
		mesh,
		dimAcceleration
		);

	
	volScalarField bTwoX(
		IOobject
		("bTwoX",
		runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
		mesh,
		dimAcceleration
		);
		
	volScalarField bTwoY(
		IOobject
		("bTwoY",
		runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
		mesh,
		dimAcceleration 
		);


	volScalarField bThreeX(
		IOobject
		("bThreeX",
		runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
		mesh,
		dimAcceleration
		);

	volScalarField bThreeY(
		IOobject
		("bThreeY",
		runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
		mesh,
		dimAcceleration
		);


	volScalarField bFourX(
		IOobject
		("bFourX",
		runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
		mesh,
		dimAcceleration
		);

	volScalarField bFourY(
		IOobject
		("bFourY",
		runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
		mesh,
		dimAcceleration 
		);
*/

	Info<< "\nReading Equations\n" << endl;
	
/*	IFstream tfIF1(path/"constant/nonInertial");
	const dictionary nonInertial(tfIF1);
	
	IFstream tfIF2(path/"constant/transportProperties");
	const dictionary transportPropertiesEq(tfIF2);
	
	*/
	
	    // Add a dictionary source
        IOdictionary nonInertial
        (
            IOobject
            (
                "nonInertial",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        eqns.addSource(nonInertial);
	
	//Define a few new dimensions 	
	dimensionSet dimInvTime 	= dimensionSet(0,0,-1,0,0,0,0);
	dimensionSet dimInvTimeSqd = dimensionSet(0,0,-2,0,0,0,0);
	
	
	//Declaration of scalar components describing frame motion
	dimensionedScalar A1("A1", dimAcceleration, 0);
	dimensionedScalar A2("A2", dimAcceleration, 0);
	dimensionedScalar A3("A3", dimAcceleration, 0); 
	
	dimensionedScalar Omega1("Omega1", dimInvTime, 0);
	dimensionedScalar Omega2("Omega2", dimInvTime, 0);
	dimensionedScalar Omega3("Omega3", dimInvTime, 0); 
	
	dimensionedScalar dOmega1("dOmega1", dimInvTimeSqd, 0);
	dimensionedScalar dOmega2("dOmega2", dimInvTimeSqd, 0);
	dimensionedScalar dOmega3("dOmega3", dimInvTimeSqd, 0); 
	
	//Read scalar components describing frame motion from dictionary
	eqns.readEquation(nonInertial, "aX");
	eqns.readEquation(nonInertial, "aY");
	eqns.readEquation(nonInertial, "aZ");
	
	eqns.readEquation(nonInertial, "omegaX");
	eqns.readEquation(nonInertial, "omegaY");
	eqns.readEquation(nonInertial, "omegaZ");
	
	eqns.readEquation(nonInertial, "dOmegaX");
	eqns.readEquation(nonInertial, "dOmegaY");
	eqns.readEquation(nonInertial, "dOmegaZ");	
	
	/*
	//Read temporary validation equations
	eqns.readEquation(nonInertial, "b1X");
	eqns.readEquation(nonInertial, "b1Y");
	eqns.readEquation(nonInertial, "b2X");
	eqns.readEquation(nonInertial, "b2Y");
	eqns.readEquation(nonInertial, "b3X");
	eqns.readEquation(nonInertial, "b3Y");
	eqns.readEquation(nonInertial, "b4X");
	eqns.readEquation(nonInertial, "b4Y");
	
	//Get indices of temporary validation equations 
	int bOneXindex		=	eqns.lookup("b1X");
	int bOneYindex		=	eqns.lookup("b1Y");
	int bTwoXindex		=	eqns.lookup("b2X");
	int bTwoYindex		=	eqns.lookup("b2Y");
	int bThreeXindex	=	eqns.lookup("b3X");
	int bThreeYindex	=	eqns.lookup("b3Y");
	int bFourXindex	=	eqns.lookup("b4X");
	int bFourYindex	=	eqns.lookup("b4Y");
	*/
	
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


	Info<< "\nStarting time loop\n" << endl;

	while (runTime.run())
	{

		#include "readPISOControls.H"
		#include "readTimeControls.H"
		#include "CourantNo.H"		
		#include "setDeltaT.H"
	
		runTime++;		
		Info<< "Time = " << runTime.timeName() << nl << endl;	
		
		
		//Assemble vectors from scalar components read from dictionary
		Info<< "\nAssembling Vectors for non-Inertial Forces\n" << nl << endl;
		
		A1 = eqns.evaluateDimensionedScalar("aX");
		A2 = eqns.evaluateDimensionedScalar("aY");
		A3 = eqns.evaluateDimensionedScalar("aZ");

		Info<< "\nCreating Angular Momentum Vector Components\n" << nl << endl;		
		Omega1 = eqns.evaluateDimensionedScalar("omegaX");
		Omega2 = eqns.evaluateDimensionedScalar("omegaY");
		Omega3 = eqns.evaluateDimensionedScalar("omegaZ");

		Info<< "\nCreating Angular Acceleration Vector Components\n" << nl << endl;
		dOmega1 = eqns.evaluateDimensionedScalar("dOmegaX");
		dOmega2 = eqns.evaluateDimensionedScalar("dOmegaY");
		dOmega3 = eqns.evaluateDimensionedScalar("dOmegaZ");
		
		Info<< "\nAssembling Vectors\n" << nl << endl;
		dimensionedVector  rectAccel 	= A1*vector(1,0,0) + A2*vector(0,1,0) + A3*vector(0,0,1);
		dimensionedVector  Omega 		= Omega1*vector(1,0,0) + Omega2*vector(0,1,0) + Omega3*vector(0,0,1);
		dimensionedVector  dOmega		= dOmega1*vector(1,0,0) + dOmega2*vector(0,1,0) + dOmega3*vector(0,0,1);	
		
		Info<< "\nVector Assembly Complete\n" << nl << endl;

		/*		
		//Evaluate temporary source equations for validation
			
		eqns.evaluateField(bOneXindex, bOneX);
		eqns.evaluateField(bOneYindex, bOneY);
		eqns.evaluateField(bTwoXindex, bTwoX);
		eqns.evaluateField(bTwoYindex, bTwoY);
		eqns.evaluateField(bThreeXindex, bThreeX);
		eqns.evaluateField(bThreeYindex, bThreeY);
		eqns.evaluateField(bFourXindex, bFourX);
		eqns.evaluateField(bFourYindex, bFourY);
		
		
		//Assemble temporary vector fields for validation 		
		volVectorField bOne 		= bOneX*vector(1,0,0) 	+ bOneY*vector(0,1,0);
		volVectorField bTwo 		= bTwoX*vector(1,0,0) 	+ bTwoY*vector(0,1,0);
		volVectorField bThree 	= bThreeX*vector(1,0,0) + bThreeY*vector(0,1,0);
		volVectorField bFour 	= bFourX*vector(1,0,0) 	+ bFourY*vector(0,1,0);
		volVectorField b			= bOne + bTwo + bThree 	+ bFour; 
		*/

		/*	
			Define a placeholder variable to represent cross product in the 
			the Coriolis force because we will want to calculate 
			this implicitly as it depends on U which is unknown
		*/
		//volVectorField OmegaCrossU = (Omega ^ U);
		
		fvVectorMatrix UEqn
		(
				fvm::ddt(U)
			+ fvm::div(phi, U)
			- fvm::laplacian(nu, U)
			
			//Coriolis force/
			//-fvm::Sp(2.0,OmegaCrossU)			
			+(2.0*Omega ^ U)	
			
			//Centrifugal force
			+(Omega^(Omega ^ mesh.C()))
				
			//Euler force				
			+(dOmega ^ (mesh.C()))
				
			//Rectilinear acceleration of reference frame				
			+rectAccel	
			
			/*
			//Additional temporary terms for validation
			-b -(2.0*Omega ^ U) -(Omega^(Omega ^ mesh.C())) -(dOmega ^ (mesh.C())) -rectAccel	
			*/
		);

		solve(UEqn == -fvc::grad(p));

		// --- PISO loop

		for (int corr=0; corr<nCorr; corr++)
		{
				volScalarField rUA = 1.0/UEqn.A();

				U = rUA*UEqn.H();
				phi = (fvc::interpolate(U) & mesh.Sf())
					+ fvc::ddtPhiCorr(rUA, U, phi);

				adjustPhi(phi, U, p);

				for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
				{
					fvScalarMatrix pEqn
					(
						fvm::laplacian(rUA, p) == fvc::div(phi)
					);

					pEqn.setReference(pRefCell, pRefValue);
					pEqn.solve();

					if (nonOrth == nNonOrthCorr)
					{
						phi -= pEqn.flux();
					}
				}

				#include "continuityErrs.H"

				U -= rUA*fvc::grad(p);
				U.correctBoundaryConditions();
		}
        


		runTime.write();

		Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
				<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
				<< nl << endl;
	}

	Info<< "End\n" << endl;

	return 0;
}


// ************************************************************************* //
