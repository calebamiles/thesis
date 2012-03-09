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
	pimpleFoam

Description
	Large time-step transient solver for incompressible, flow using the PIMPLE
	(merged PISO-SIMPLE) algorithm.

	Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IOEquationReader.H"
#include "IFstream.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"
#include "IObasicSourceList.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
	#include "setRootCase.H"
	#include "createTime.H"
	#include "createMesh.H"
	#include "initContinuityErrs.H"
	#include "createFields.H"
	#include "createEquationReader.H"
	
	pimpleControl pimple(mesh);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	fileName path(args.rootPath()/args.caseName());	


	//Create time variable "t" for use in symbolic expressions with equationReader
	eqns.scalarSources().addSource(runTime.value(), "t", dimTime);
	
   //Create position variables x, y, and z for use in symbolic expressions with equationReader  
	eqns.vectorSources().addSource(mesh.C());
	
	
	Info<< "\nReading Equations\n" << endl;
	
	
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
	
	
	Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

		Info<< "Time = " << runTime.timeName() << nl << endl;

		//Assemble vectors from scalar components read from dictionary
		Info<< "\nAssembling Vectors for non-Inertial Forces\n" << endl;	
		
		A1 = eqns.evaluateDimensionedScalar("aX");
		A2 = eqns.evaluateDimensionedScalar("aY");
		A3 = eqns.evaluateDimensionedScalar("aZ");

		Info<< "\nCreating Angular Momentum Vector Components\n" << endl;		
		Omega1 = eqns.evaluateDimensionedScalar("omegaX");
		Omega2 = eqns.evaluateDimensionedScalar("omegaY");
		Omega3 = eqns.evaluateDimensionedScalar("omegaZ");

		Info<< "\nCreating Angular Acceleration Vector Components\n" << endl;
		dOmega1 = eqns.evaluateDimensionedScalar("dOmegaX");
		dOmega2 = eqns.evaluateDimensionedScalar("dOmegaY");
		dOmega3 = eqns.evaluateDimensionedScalar("dOmegaZ");
		
		Info<< "\nAssembling Vectors\n" << endl;		
		dimensionedVector  rectAccel 	= A1*vector(1,0,0) + A2*vector(0,1,0) + A3*vector(0,0,1);
		dimensionedVector  Omega 		= Omega1*vector(1,0,0) + Omega2*vector(0,1,0) + Omega3*vector(0,0,1);
		dimensionedVector  dOmega		= dOmega1*vector(1,0,0) + dOmega2*vector(0,1,0) + dOmega3*vector(0,0,1);	
		
		Info<< "\nVector Assembly Complete\n" << endl;
		
		
        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
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
