/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    SRFInitUrel

Description
    Initialize Urel field at time 0 using the SRF model

\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "fvc.H"

#include "SRFModel.H"

#include "wallFvPatch.H" // to ID wall patches

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{
    if (runTime.timeOutputValue() > 0) 
    {
        Info<< "Skipping time " << runTime.timeOutputValue() << " > 0" << endl;
        return;
    }

    IOobject Uheader
    (
        "Urel",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    if (Uheader.headerOk())
    {
        Info<< "    Reading Urel" << endl;
        volVectorField Urel(Uheader, mesh); // GeometricField<vector, fvPatchField, volMesh>
        vectorField& Uif = Urel.internalField();

        Info<< "    Creating SRF model\n" << endl;
        autoPtr<SRF::SRFModel> SRF
        (
             SRF::SRFModel::New(Urel)
        );

        // Get reference to the SRF model
        const SRF::SRFModel& srf =
            Urel.db().lookupObject<SRF::SRFModel>("SRFProperties");

        // Check internal field uniformity to get freestream velocity
        Info<< "    Checking internal field with size " << Uif.size() << endl;
        vector UInf(Uif[0]);
        forAll(Uif, i)
        {
            if (Uif[i] != UInf)
            {
                Info<<"NOTE: non-uniform internal field?" << endl;
                break;
            }
        }
        Info<< "    Uniform velocity : " << UInf << endl;

        Info<< "    Calculating internal field" << endl;
        Uif = UInf - srf.velocity(Urel.mesh().C()); // evaluate SRF at cell centers

        Info<< "    Calculating boundary fields" << endl;
        const fvPatchList& patches = mesh.boundary();
        forAll(patches, patchI)
        {
            const fvPatch& currPatch = patches[patchI];

            if (!isA<wallFvPatch>(currPatch))
            {
                Info<< "     - " << currPatch.name() << endl;
                Urel.boundaryField()[patchI] == // need '==' to force assignment, e.g. with fixed patch
                    UInf - srf.velocity(mesh.Cf().boundaryField()[patchI]); // at face centers
            }
        }

        Info<< "    Writing Urel" << endl;
        Urel.write();
    }
    else
    {
        Info<< "    No Urel" << endl;
    }

    Info<< "\nEnd\n" << endl;
}


// ************************************************************************* //
