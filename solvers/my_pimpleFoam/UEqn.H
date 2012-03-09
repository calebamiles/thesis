// Solve the Momentum equation

tmp<fvVectorMatrix> UEqn
(
    fvm::ddt(U)
  + fvm::div(phi, U)
  + turbulence->divDevReff(U)
  //Coriolis force/		
  +(2.0*Omega ^ U)	
			
  //Centrifugal force
  +(Omega^(Omega ^ mesh.C()))
				
  //Euler force				
  +(dOmega ^ (mesh.C()))
				
  //Rectilinear acceleration of reference frame				
  +rectAccel	
  ==
    sources(U)
);

UEqn().relax();

sources.constrain(UEqn());

volScalarField rAU(1.0/UEqn().A());

if (pimple.momentumPredictor())
{
    solve(UEqn() == -fvc::grad(p));
}
