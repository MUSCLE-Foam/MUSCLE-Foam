/*---------------------------------------------------------------------------*\

    MUSCLE-Foam: an OpenFOAM-based multiphysics solver for reactor physics
    Copyright (C) 2026 NRG PALLAS

\-----------------------------------------------------------------------------/

License
    This file is part of MUSCLE-Foam. It is published under the terms of the
    GNU General Public License. See <http://www.gnu.org/licenses/>.

Application
    generateCrossSectionFiles

Description
    This utility can be used to quickly generate cross-section files when
    setting up a MUSCLE-Foam case. Specify the number of energy groups and the
    neutronics solver and the script will generate all files. Each cross-
    section file will have the name of the file plus "_VALUE" as an initial
    internalField value. These values should be replaced by the appropriate
    values before running a simulation.

Usage
    \b generateCrossSectionFiles [OPTION]

      - \par -energyGroups <n of energy groups>
        Number of energy groups

      - \par -neutronicsSolver <name>
        Neutronics solver name

      - \par -noAlpha
        Do not generate thermal feedback coefficient files

    Example usage:
      - To generate cross-section files for 6 energy groups in the SP3 solver
        run the following from the command-line in your case directory:
        \verbatim
            generateCrossSectionFiles -energyGroups 6 -neutronicsSolver SP3
        \endverbatim

      - To run the utility based on the entries in controlDict:
        \verbatim
            neutronicsSolver diffusion;

            energyGroups 7;
        \endverbatim
        then execute \c generateCrossSectionFiles

\*---------------------------------------------------------------------------*/

#include "fv.H"
#include "argList.H"
#include "IOstreams.H"
#include "dimensionSets.H"
#include "OFstream.H"
#include "Time.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "energyGroups",
        "number",
        "Number of energy groups"
    );

    argList::addOption
    (
        "neutronicsSolver",
        "name",
        "Neutronics solver name"
    );

    argList::addBoolOption
    (
        "noAlpha",
        "Do not generate thermal feedback coefficient files"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    label energyGroups
    (
        runTime.controlDict().lookupOrDefault("energyGroups", 0)
    );

    // Read the neutronicsSolverName from the optional solver entry in
    // controlDict
    word neutronicsSolverName
    (
        runTime.controlDict().lookupOrDefault("neutronicsSolver", word::null)
    );

    bool noAlpha(args.optionFound("noAlpha"));

    // Optionally reset the energyGroups from the -energyGroups command-line
    // argument
    args.optionReadIfPresent("energyGroups", energyGroups);

    // Optionally reset the neutronicsSolver name from the -solver
    // command-line argument
    args.optionReadIfPresent("neutronicsSolver", neutronicsSolverName);

    // Directory setup
    fileName zeroDir = "0";

    if (!isDir(zeroDir))
    {
        Info << "Creating 0/ directory...\n";
        mkDir(zeroDir);
    }

    // Helper function to write a file
    auto writeFile = [](const fileName fullPath, const word name, const dimensionSet dim)
    {
        // Open file for writing
        OFstream file(fullPath);

        // Write OpenFOAM header
        file << "// ************************************************************************** //\n";
        file << "FoamFile\n";
        file << "{\n";
        file << "    format      ascii;\n";
        file << "    class       volScalarField;\n";
        file << "    location    \"0\";\n";
        file << "    object      " << name << ";\n";
        file << "}\n";
        file << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";

        // Write dimensions
        file << "dimensions      " << dim << ";\n\n";

        // Write internalField
        file << "internalField   uniform 0;\n\n";

        // Write empty boundaryField
        file << "boundaryField\n";
        file << "{\n";
        file << "    #include \"XS\"\n";
        file << "}\n\n";
        file << "// ************************************************************************** //";

        Info << "Created file: " << fullPath << "\n";
    };

    // Scattering cross-sections
    for (int i = 0; i < energyGroups; i++)
    {
        for (int j = 0; j < energyGroups; j++)
        {
            word S = "S"+name(i)+name(j)+"_0";
            writeFile(zeroDir / S, S, inv(dimLength));

            if (!noAlpha)
            {
                word alphaS = "alphaS"+name(i)+name(j);
                writeFile(zeroDir / alphaS, alphaS, inv(dimLength));
            }

            if (neutronicsSolverName == "SP3")
            {
                word S2 = "S2"+name(i)+name(j)+"_0";
                writeFile(zeroDir / S2, S2, inv(dimLength));

                word S3 = "S3"+name(i)+name(j)+"_0";
                writeFile(zeroDir / S3, S3, inv(dimLength));

                if (!noAlpha)
                {
                    word alphaS2 = "alphaS2"+name(i)+name(j);
                    writeFile(zeroDir / alphaS2, alphaS2, inv(dimLength));

                    word alphaS3 = "alphaS3"+name(i)+name(j);
                    writeFile(zeroDir / alphaS3, alphaS3, inv(dimLength));
                }
            }
        }
    }

    // Other cross-sections
    wordList crossSectionNames({"FN", "M", "SP", "TR", "IV"});

    forAll(crossSectionNames, c)
    {
        for (int i = 0; i < energyGroups; i++)
        {
            word XS = crossSectionNames[c]+name(i)
                + (crossSectionNames[c] != "IV" ? "_0" : "");
            if (crossSectionNames[c] == "SP")
            {
                writeFile(zeroDir / XS, XS, dimForce);
            }
            else if (crossSectionNames[c] == "IV")
            {
                writeFile(zeroDir / XS, XS, inv(dimVelocity));
            }
            else
            {
                writeFile(zeroDir / XS, XS, inv(dimLength));
            }

            if (!noAlpha)
            {
                word alphaXS = "alpha"+crossSectionNames[c]+name(i);
                if (crossSectionNames[c] == "SP")
                {
                    writeFile(zeroDir / alphaXS, alphaXS, dimForce);
                }
                else if (crossSectionNames[c] != "IV")
                {
                    writeFile(zeroDir / alphaXS, alphaXS, inv(dimLength));
                }
            }
        }
    }

    // Prompt/delayed yields
    for (int i = 0; i < energyGroups; i++)
    {
        word chiPrompt = "chiPrompt"+name(i);
        writeFile(zeroDir / chiPrompt, chiPrompt, dimless);

        word chiDel = "chiDel"+name(i);
        writeFile(zeroDir / chiDel, chiDel, dimless);
    }

    // Reference density and temperatures
    writeFile(zeroDir / "rhoXS", "rhoXS", dimDensity);
    writeFile(zeroDir / "TXS", "TXS", dimTemperature);

    // Misc files
    writeFile(zeroDir / "DAlbedo", "DAlbedo", dimLength);
    writeFile(zeroDir / "fluxAlbedo", "fluxAlbedo", inv(dimArea*dimTime));

    // Write XS file
    OFstream XSfile(zeroDir / "XS");

    XSfile << "\".*\"\n";
    XSfile << "{\n";
    XSfile << "    type    zeroGradient;\n";
    XSfile << "}";

    Info << "Initialization complete. Files created in 0/ directory.\n";

    return 0;
}
