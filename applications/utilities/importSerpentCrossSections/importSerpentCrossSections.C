/*---------------------------------------------------------------------------*\

    MUSCLE-Foam: an OpenFOAM-based multiphysics solver for reactor physics
    Copyright (C) 2026 NRG PALLAS

\-----------------------------------------------------------------------------/

License
    This file is part of MUSCLE-Foam. It is published under the terms of the
    GNU General Public License. See <http://www.gnu.org/licenses/>.

Application
    importSerpentCrossSections

Description
    This utility can be used to copy Serpent cross-sections to a setFieldsDict
    which can be used to set cross-sections in different predefined universes.

Usage
    \b importSerpentCrossSections [OPTION]

      - \par -energyGroups <n of energy groups>
        Number of energy groups

      - \par -outputFile <Serpent output file>
        Serpent output file name

      - \par -neutronicsSolver <name>
        Neutronics solver name

      - \par -noAlpha
        Do not generate thermal feedback coefficient files

    Example usage:
      - To copy cross-section data for the diffusion solver with 6 energy
        groups, run the following from the command-line in your case directory:
        \verbatim
            importSerpentCrossSections -energyGroups 6 -outputFile case_res.m \
                -neutronicsSolver SP3
        \endverbatim

      - To run the utility based on the entries in controlDict:
        \verbatim
            neutronicsSolver diffusion;

            energyGroups 7;
        \endverbatim
        then execute \c generateCrossSectionFiles -outputFile case_res.m

\*---------------------------------------------------------------------------*/

#include "fv.H"
#include "argList.H"
#include "IOstreams.H"
#include "dimensionSets.H"
#include "OFstream.H"
#include "IFstream.H"
#include "Time.H"

using namespace Foam;

// Conversion factors for length (cm to m) and energy (MeV to J)
const static scalar lengthConversion = 100.0;
const static scalar energyConversion = 1.60218e-13;

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
        "outputFile",
        "name",
        "Serpent output file name"
    );

    argList::addOption
    (
        "neutronicsSolver",
        "name",
        "Neutronics solver name"
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

    // Optionally reset the energyGroups from the -energyGroups command-line
    // argument
    args.optionReadIfPresent("energyGroups", energyGroups);

    // Optionally reset the neutronicsSolver name from the -solver
    // command-line argument
    args.optionReadIfPresent("neutronicsSolver", neutronicsSolverName);

    // Serpent output file
    fileName outputFile = args.optionRead<fileName>("outputFile");
    IFstream serpentOutput(outputFile);

    // Dict file
    fileName setFieldsDict = "system/setFieldsDictXS";

    // Output stream
    OFstream file(setFieldsDict);

    // Helper function to write a cross-section entry
    auto writeXS = []
    (
        OFstream& file,
        const string line,
        const word XSname,
        const label energyGroups,
        bool suffix,
        bool lenConv,
        bool eneConv
    )
    {
        size_t start = line.rfind('[');
        size_t end   = line.rfind(']');

        if (start == std::string::npos || end == std::string::npos)
        {
            FatalErrorInFunction
                << "Could not find square brackets in line: "
                << nl << line
                << nl << XSname
                << nl << start
                << nl << end
                << exit(FatalError);
        }

        string data = line.substr(start + 1, end - start - 1);
        std::stringstream ss(data);

        List<scalar> xs(energyGroups, Zero);

        scalar value;
        label index = 0;
        label i = 0;

        while (ss >> value)
        {
            // Only store every other value (cross-sections, not uncertainties)
            if (index % 2 == 0)
            {
                if (i >= energyGroups)
                {
                    FatalErrorInFunction
                        << "More cross-section entries found than expected ("
                        << energyGroups << ")"
                        << exit(FatalError);
                }
                xs[i++] = value * (lenConv ? lengthConversion : 1.0)
                    * (eneConv ? energyConversion : 1.0);
            }
            index++;
        }

        if (i < energyGroups)
        {
            FatalErrorInFunction
                << "Fewer cross-section entries found (" << i
                << ") than expected (" << energyGroups << ")"
                << exit(FatalError);
        }

        forAll(xs, e)
        {
            file << "            volScalarFieldValue "<< XSname << e;
            if (suffix)
            {
                file << "_0";
            }
            file << "    " << xs[e] << "\n";
        }

        file << "\n";
    };

    // Write OpenFOAM header
    file << "// ************************************************************************** //\n";
    file << "FoamFile\n";
    file << "{\n";
    file << "    format      ascii;\n";
    file << "    class       dictionary;\n";
    file << "    location    \"system\";\n";
    file << "    object      setFieldsDict;\n";
    file << "}\n";
    file << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";

    file << "regions\n";
    file << "(\n";

    // Additional variables
    string line;
    word universe;
    bool readyForNextUniv = true;
    List<scalar> fiss(energyGroups*2, Zero);
    bool dnpRead = false;
    label DNPGroups = -1;
    List<scalar> betaEff;
    List<scalar> lambda;

    // Loop through Serpent output file lines
    while (std::getline(serpentOutput.stdStream(), line))
    {
        // Find universe name
        if (line.find("GC_UNIVERSE_NAME") != std::string::npos)
        {
            if (!readyForNextUniv)
            {
                FatalErrorInFunction
                << "New universe found in line: " << nl
                << line << nl
                << "while universe " << universe
                << " not fully initialized yet."
                << exit(FatalError);
            }

            // Not ready for next universe untill all XS are set
            readyForNextUniv = false;

            size_t firstQuote = line.find('\'');
            size_t secondQuote = line.find('\'', firstQuote + 1);

            if
            (
                   firstQuote != std::string::npos
                && secondQuote != std::string::npos
            )
            {
                string univ = line.substr
                (
                    firstQuote + 1,
                    secondQuote - firstQuote - 1
                );

                universe = word(univ);

                Info << "Universe: " << universe << nl;
            }
            else
            {
                FatalErrorInFunction
                << "Failed to find universe name in line: " << nl
                << line << nl
                << exit(FatalError);
            }

            // New universe entry
            file << "    zoneToCell\n";
            file << "    {\n";
            file << "        zone univ" << universe << ";\n";
            file << "        fieldValues \n";
            file << "        (\n";
        }

        // Total cross-sections
        if (line.find("INF_TOT") != std::string::npos)
        {
            writeXS(file, line, "TOT", energyGroups, true, true, false);
        }

        // Fission XS times neutrons per fission
        if (line.find("INF_NSF") != std::string::npos)
        {
            writeXS(file, line, "FN", energyGroups, true, true, false);
        }

        // Inverse neutron velocities
        if (line.find("INF_INVV") != std::string::npos)
        {
            writeXS(file, line, "IV", energyGroups, false, true, false);
        }

        // Transport cross-sections
        if (line.find("INF_TRANSPXS") != std::string::npos)
        {
            writeXS(file, line, "TR", energyGroups, true, true, false);
        }

        // Removal cross-sections
        if (line.find("INF_REMXS") != std::string::npos)
        {
            writeXS(file, line, "M", energyGroups, true, true, false);
        }

        // Prompt yield
        if (line.find("INF_CHIP") != std::string::npos)
        {
            writeXS(file, line, "chiPrompt", energyGroups, false, false, false);
        }

        // Delayed yield
        if (line.find("INF_CHID") != std::string::npos)
        {
            writeXS(file, line, "chiDel", energyGroups, false, false, false);
        }

        // Scattering XS
        for (int ord = 0; ord < 8; ord++)
        {
            word scatteringName("INF_S"+Foam::name(ord));

            if (line.find(scatteringName) != std::string::npos)
            {
                size_t start = line.rfind('[');
                size_t end   = line.rfind(']');

                if (start == std::string::npos || end == std::string::npos)
                {
                    FatalErrorInFunction
                        << "Could not find square brackets in line: "
                        << nl << line
                        << exit(FatalError);
                }

                string data = line.substr(start + 1, end - start - 1);
                std::stringstream ss(data);

                for (int i = 0; i < energyGroups; i++)
                {
                    string energyGroupData = "[";

                    for (int j = 0; j < energyGroups; j++)
                    {
                        scalar s, unc;

                        if (!(ss >> s >> unc))
                        {
                            FatalErrorInFunction
                                << "Unexpected end of data when reading "
                                << "scattering matrix, group "
                                << i << " target " << j << exit(FatalError);
                        }

                        energyGroupData += Foam::name(s) + " "
                            + Foam::name(unc) + " ";
                    }

                    energyGroupData += "]";

                    writeXS
                    (
                        file,
                        energyGroupData,
                        ord > 0 ?
                          word("S"+Foam::name(ord)+Foam::name(i))
                        : word("S"+Foam::name(i)),
                        energyGroups,
                        true,
                        true,
                        false
                    );
                }

                if (ord == 7)
                {
                    // End universe entry
                    file << "        );\n";
                    file << "    }\n\n";

                    readyForNextUniv = true;
                }
            }
        }

        // Fission cross-section times energy per fission
        if (line.find("INF_FISS") != std::string::npos)
        {
            size_t start = line.rfind('[');
            size_t end   = line.rfind(']');

            if (start == std::string::npos || end == std::string::npos)
            {
                FatalErrorInFunction
                    << "Could not find square brackets in line: "
                    << nl << line
                    << exit(FatalError);
            }

            string data = line.substr(start + 1, end - start - 1);
            std::stringstream ss(data);

            forAll(fiss, i)
            {
                if (!(ss >> fiss[i]))
                {
                    FatalErrorInFunction
                        << "Unexpected end of data when reading "
                        << "fission cross-sections."
                        << exit(FatalError);
                }
            }
        }
        if (line.find("INF_KAPPA") != std::string::npos)
        {
            size_t start = line.rfind('[');
            size_t end   = line.rfind(']');

            if (start == std::string::npos || end == std::string::npos)
            {
                FatalErrorInFunction
                    << "Could not find square brackets in line: "
                    << nl << line
                    << exit(FatalError);
            }

            string data = line.substr(start + 1, end - start - 1);
            std::stringstream ss(data);

            string SPData = "[";

            scalar kappa, unc;

            for (int i = 0; i < energyGroups; i++)
            {
                if (!(ss >> kappa >> unc))
                {
                    FatalErrorInFunction
                        << "Unexpected end of data when reading "
                        << "kappa values."
                        << nl << line
                        << nl << data
                        << nl << i
                        << nl << kappa
                        << nl << unc
                        << exit(FatalError);
                }

                fiss[i*2] *= kappa;

                SPData += Foam::name(fiss[i*2]) + " 0.0 ";
            }

            SPData += "]";

            writeXS
            (
                file,
                SPData,
                "SP",
                energyGroups,
                true,
                true,
                true
            );
        }

        // Delayed neutron precursors
        if (!dnpRead)
        {
            if (line.find("PRECURSOR_GROUPS") != std::string::npos)
            {
                std::size_t eqPos = line.find('=');
                std::size_t semicolonPos = line.find(';', eqPos);

                if
                (
                       eqPos != std::string::npos
                    && semicolonPos != std::string::npos
                )
                {
                    string valueStr = line.substr
                    (
                        eqPos + 1, semicolonPos - eqPos - 1
                    );
                    std::stringstream ss(valueStr);
                    ss >> DNPGroups;

                    if (!ss)
                    {
                        FatalErrorInFunction
                            << "Failed to parse PRECURSOR_GROUPS from line: " << line
                            << exit(FatalError);
                    }
                }
                else
                {
                    FatalErrorInFunction
                        << "Malformed PRECURSOR_GROUPS line: " << line
                        << exit(FatalError);
                }
            }

            if (line.find("ADJ_NAUCHI_BETA_EFF") != std::string::npos)
            {
                std::size_t lb = line.rfind('[');
                std::size_t rb = line.rfind(']');

                if (lb == std::string::npos || rb == std::string::npos)
                {
                    FatalErrorInFunction
                        << "Malformed ADJ_NAUCHI_BETA_EFF line:\n" << line
                        << exit(FatalError);
                }

                string contents = line.substr(lb + 1, rb - lb - 1);
                std::stringstream ss(contents);
                scalar val, unc;

                betaEff.setSize(DNPGroups);

                for (label i = 0; i < DNPGroups; i++)
                {
                    if (!(ss >> val >> unc))
                    {
                        FatalErrorInFunction
                            << "Failed reading BETA_EFF group " << i
                            << " from line:\n" << line
                            << exit(FatalError);
                    }
                    betaEff[i] = val;
                }
            }

            if (line.find("ADJ_NAUCHI_LAMBDA") != std::string::npos)
            {
                std::size_t lb = line.rfind('[');
                std::size_t rb = line.rfind(']');

                if (lb == std::string::npos || rb == std::string::npos)
                {
                    FatalErrorInFunction
                        << "Malformed ADJ_NAUCHI_LAMBDA line:\n" << line
                        << exit(FatalError);
                }

                string contents = line.substr(lb + 1, rb - lb - 1);
                std::stringstream ss(contents);
                scalar val, unc;

                lambda.setSize(DNPGroups);

                for (label i = 0; i < DNPGroups; i++)
                {
                    if (!(ss >> val >> unc))
                    {
                        FatalErrorInFunction
                            << "Failed reading LAMBDA group " << i
                            << " from line:\n" << line
                            << exit(FatalError);
                    }
                    lambda[i] = val;
                }

                dnpRead = true;
            }
        }
    }

    file << ");\n\n";
    file << "// ************************************************************************* //";

    Info << "setFieldsDict created from " << outputFile << endl;

    // Delayed neutron data
    fileName delayedNeutronConstants = "constant/delayedNeutronConstants";

    // Output stream
    OFstream dncFile(delayedNeutronConstants);

    // Write OpenFOAM header
    dncFile << "// ************************************************************************** //\n";
    dncFile << "FoamFile\n";
    dncFile << "{\n";
    dncFile << "    format      ascii;\n";
    dncFile << "    class       dictionary;\n";
    dncFile << "    location    \"constant\";\n";
    dncFile << "    object      delayedNeutronConstants;\n";
    dncFile << "}\n";
    dncFile << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";

    dncFile << "DNPGroups   " << DNPGroups << ";\n\n";

    forAll(lambda, i)
    {
        dncFile << "lambda" << i << "     " << lambda[i] << ";\n";
    }

    dncFile << "\n";

    forAll(betaEff, i)
    {
        dncFile << "beta" << i << "       " << betaEff[i] << ";\n";
    }

    dncFile << "\n";
    dncFile << "// ************************************************************************* //";

    Info << "delayedNeutronConstants created from " << outputFile << endl;

    return 0;
}
