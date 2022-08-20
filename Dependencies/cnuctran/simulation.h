/*

      This file is part of the CNUCTRAN library

      @license  MIT
      @link     https://github.com/rabieomar92/cnuctran

      Copyright (c) 2022 M. R. Omar

      This header file contains the definitions of all methods for preparing
      the simulation associated to the proposed probabilistic method.

 */

#ifndef SIMULATION_H
#define SIMULATION_H


#include <iostream>
#include <iomanip>
#include <mpreal.h>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <pugixml/pugixml.hpp>
#include <cnuctran/smatrix.h>
#include <cnuctran/solver.h>
#include <cnuctran/cnuctran.h>
#include <ctime>

using namespace pugi;
using namespace mpfr;

namespace cnuctran {


    class simulation
    {

    public:

        static void build_chains(solver& s, std::map<std::string, std::map<std::string, mpreal>>& rxn_rates,
            std::string xml_data_location)
        {
            std::vector<std::string> species_names = s.species_names;

//..........Loads the nuclides data from the XML source file. 
            xml_document file;
            xml_parse_result open_success = file.load_file(xml_data_location.c_str());
            if (!open_success)
            {
                std::cout << "INFO\t<cnuctran::depletion_scheme::build_chains(...)> Nuclides data file is not provided." << std::endl;
                return;
            }

//..........Unpack the XML root.
            xml_node root = file.child("depletion");

            for (xml_node species : root.children())
            {
                std::string species_name = species.attribute("name").value();
                std::vector<std::string>::iterator it = std::find(species_names.begin(), species_names.end(), species_name);
                if (it == species_names.end())
                    continue;

                mpreal decay_rate;
                if (species.attribute("half_life"))
                    decay_rate = mpfr::log(mpreal("2")) / mpreal(species.attribute("half_life").value());
                else
                    decay_rate = mpreal("0");

                for (xml_node removal : species.children())
                {
                    if (std::string(removal.name()) == "decay_type")
                    {
                        mpreal decay_rate_adjusted = mpreal(removal.attribute("branching_ratio").value()) * decay_rate;
                        std::string parent = species_name;
                        std::string daughter = removal.attribute("target").value();
                        std::vector<std::string>::iterator it_parent = std::find(species_names.begin(), species_names.end(), parent);
                        int parent_id = distance(species_names.begin(), it_parent);
                        std::vector<std::string>::iterator it_daughter = std::find(species_names.begin(), species_names.end(), daughter);
                        if (it_daughter != species_names.end())
                        {
                            int daughter_id = distance(species_names.begin(), it_daughter);
                            s.add_removal(parent_id, decay_rate_adjusted, std::vector<int>({ daughter_id }));
                        }
                        else
                        {
                            s.add_removal(parent_id, decay_rate_adjusted, std::vector<int>({ __nop__ }));
                        }

                    }

                    if (rxn_rates.size() != 0)
                    {
                        if (rxn_rates.count(species_name))
                        {
                            if (std::string(removal.name()) == "reaction_type" && removal.attribute("target"))
                            {
                                std::string parent = species_name;
                                std::vector<std::string>::iterator it_parent = std::find(species_names.begin(), species_names.end(), parent);
                                int parent_id = distance(species_names.begin(), it_parent);
                                if (rxn_rates[parent].count(removal.attribute("type").value()) &&
                                    removal.attribute("type").value() != "fission")
                                {
                                    std::string daughter = removal.attribute("target").value();
                                    mpreal removal_rate = mpreal(rxn_rates[parent][removal.attribute("type").value()]);
                                    std::vector<std::string>::iterator it_daughter = std::find(species_names.begin(), species_names.end(), daughter);
                                    if (it_daughter != species_names.end())
                                    {
                                        int daughter_id = distance(species_names.begin(), it_daughter);
                                        s.add_removal(parent_id, removal_rate, std::vector<int>({ daughter_id }));
                                    }
                                    else
                                    {
                                        s.add_removal(parent_id, removal_rate, std::vector<int>({ __nop__ }));
                                    }
                                }
                            }

                            if (std::string(removal.name()) == "neutron_fission_yields")
                            {
                                std::string parent = species_name;
                                std::vector<std::string>::iterator it_parent = std::find(species_names.begin(), species_names.end(), parent);
                                int parent_id = distance(species_names.begin(), it_parent);
                                mpreal energy = mpreal("0");
                                std::vector<std::string> products;
                                std::vector<mpreal> yields;
                                if (rxn_rates[parent].count("fission") != 0)
                                {
                                    for (xml_node data : removal.children())
                                    {
                                        if (std::string(data.name()) == "energies")
                                        {
                                            std::stringstream ss(data.child_value()); std::string token;
                                            std::vector<std::string> energies;
                                            while (std::getline(ss, token, ' '))
                                                if (token != "")
                                                    energies.push_back(token);
                                            energy = energies[0];
                                        }

                                        if (std::string(data.name()) == "fission_yields")
                                        {
                                            if (mpreal(data.attribute("energy").value()) == energy)
                                            {
                                                for (xml_node param : data.children())
                                                {
                                                    if (std::string(param.name()) == "products")
                                                    {
                                                        std::stringstream ss(param.child_value()); std::string token;
                                                        while (std::getline(ss, token, ' '))
                                                            if (token != "")
                                                                products.push_back(token);
                                                    }

                                                    if (std::string(param.name()) == "data")
                                                    {
                                                        std::stringstream ss(param.child_value()); std::string token;
                                                        while (std::getline(ss, token, ' '))
                                                            if (token != "")
                                                                yields.push_back(mpreal(token));
                                                    }

                                                }

                                                mpreal total_fission_rate = rxn_rates[parent]["fission"];
                                                std::vector<mpreal> yields_to_add;
                                                std::vector<int> daughters_id_to_add;
                                                std::vector<std::string>::iterator it_product;
                                                for (std::string product : products)
                                                {

                                                    it_product = find(species_names.begin(), species_names.end(), product);
                                                    if (it_product != species_names.end())
                                                    {
                                                        int product_id = distance(species_names.begin(), it_product);
                                                        daughters_id_to_add.push_back(product_id);
                                                        it_product = find(products.begin(), products.end(), product);
                                                        product_id = distance(products.begin(), it_product);
                                                        yields_to_add.push_back(yields[product_id]);
                                                    }
                                                }
                                                it_parent = find(species_names.begin(), species_names.end(), parent);
                                                parent_id = distance(species_names.begin(), it_parent);

                                                s.add_removal(parent_id, total_fission_rate, daughters_id_to_add, yields_to_add);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            return;
        }

        //Gets the nuclide symbols i.e. Am241, Pu240 etc. based on the presecribed parameters.
        //Amin = Minimum atomic number.
        //Amax = Maximum atomic number.
        //xml_data_location = the location of the nuclides data file.
        static std::vector<std::string> get_nuclide_names(std::string xml_data_location, int AMin = -1, int AMax = -1)
        {
            xml_document file;
            xml_parse_result open_success = file.load_file(xml_data_location.c_str());
            if (!open_success)
            {
                std::cout << "WARNING\t<cnuctran::depletion_scheme::get_nuclide_names(...)> Fail retrieving data from " << xml_data_location << "." << std::endl;
                return std::vector<std::string>();
            }


            xml_node root = file.child("depletion");

            std::vector<std::string> species_names;

            for (xml_node species : root.children())
            {
                std::string name = species.attribute("name").value();
                std::stringstream ss(name); std::string token;
                std::getline(ss, token, '_');
                std::string x = "";
                for (char c : token)
                    if (isdigit(c)) x += c;
                if (AMin == AMax == -1)
                    species_names.push_back(name);
                else
                {
                    int A = std::stoi(x);
                    if (A >= AMin && A <= AMax)
                        species_names.push_back(name);
                }
            }
            return species_names;
        }


        /*
            Reads the input XML file (input.xml) and obtains all simulation parameters. Finally, this
            routine runs the simulation.
        */
        static void from_input(std::string xml_input)
        {

            try
            {
//..............Loads the XML input file. Handle loading errors.
                xml_document input_file;
                xml_parse_result open_success = input_file.load_file(xml_input.c_str());
                if (!open_success) {
                    std::cout << "FATAL-ERROR\t<cnuctran::simulation::from_input()> " << open_success.description() << std::endl;
                    exit(1);
                }

                xml_node root = input_file.child("problem");

//..............Reads simulation parameters.

                mpreal n;
                mpreal t;
                int precision_digits = 400;
                int output_digits = 15;
                double min_rate;
                const char_t* tmp;
                int AMin = -1;
                int AMax = -1;
                std::string output_location = "";

                //Obtains the verbosity level from the input file.
                tmp = root.child("simulation_params").child("verbosity").child_value();
                tmp != "" ? __vbs__ = std::stoi(tmp) : __vbs__ = 0;

                //Obtains the precision digits from the input file.
                tmp = root.child("simulation_params").child("precision_digits").child_value();
                tmp != "" ? precision_digits = std::stoi(tmp) : precision_digits = __dps__;
                if (precision_digits < 30)
                {
                    precision_digits = 30;
                    std::cout << "WARNING\t<cnuctran::simulation::from_input()> A precision < 30 digits is vulnerable to errorneous arithmetics that lead to fatal error." << std::endl;
                    std::cout << "The precision was set to 30 digits." << std::endl;
                }

                //Obtains epsilon from the user.
                tmp = root.child("simulation_params").child("epsilon").child_value();
                if (tmp != "") __eps__ = mpreal(tmp, digits2bits(__dps__));
                if (__eps__ < mpreal("0.0", digits2bits(__dps__))) std::cout << "WARNING\t<cnuctran::simulation::from_input()> Epsilon must be greater than zero.\n";

                //Obtains the minimum rate from the input file.
                tmp = root.child("simulation_params").child("min_rate").child_value();
                if (tmp != "") __mnr__ = mpreal(tmp).toDouble();

                //Obtains the minimum rate from the input file.
                tmp = root.child("simulation_params").child("max_rate").child_value();
                if (tmp != "") __mxr__ = mpreal(tmp).toDouble();

                //Obtains the output precision digits from the input file.
                tmp = root.child("simulation_params").child("output_digits").child_value();
                tmp != "" ? output_digits = std::stoi(tmp) : output_digits = __dop__;


                //IMPORTANT! Reads the precision before declaring any high-precision float vars.
                mpreal::set_default_prec(digits2bits(precision_digits));
                std::cout.precision(output_digits);

                //Reads the substep interval in seconds.
                tmp = root.child("simulation_params").child("n").child_value();
                tmp != "" ? n = mpreal(tmp) : throw (int)errex::MISSING_SUBSTEP_SIZE;

                //Reads the time step in seconds.
                tmp = root.child("simulation_params").child("time_step").child_value();
                tmp != "" ? t = mpreal(tmp) : throw (int)errex::MISSING_STEP_SIZE;

                //Reads the final councentration output file location.
                tmp = root.child("simulation_params").child("output").child_value();
                tmp != "" ? output_location = tmp : output_location = ".//output.xml";
                if (__vbs__) std::cout << "Final nuclide concentrations will be written in " << output_location << std::endl;



//..............LOOP OVER ALL ZONES and read the species, initial concentrations and reaction rates (for the zone).
                int nzones = 0;
                //Initialize file stream of the output file (output.xml).
                std::ofstream file;
                file.open(output_location, std::ios::out);
                file << "<output>" << std::endl;

                // Initialize file stream of the user-friendly output file (output.out).
                std::ofstream uffile;
                uffile.open("./output.out", std::ios::out);
                uffile << "Final nuclides concentration computed using CNUCTRAN" << std::endl;
                auto current_t = std::time(nullptr);
                auto tm = *std::localtime(&current_t);
                uffile << "Date: " << std::put_time(&tm, "%d-%m-%Y %H-%M-%S") << std::endl << std::endl;

                //Loop over all relevant XML child nodes for each zone.
                for (xml_node zone : root.children())
                {

                    // Filters XML nodes that are not a zone node.
                    if (std::string(zone.name()) != "zone") continue;

                    // Printing concentrations table header for each calculation zone.
                    // This table will be printed into ./output.out.
                    uffile << "Zone: " << std::string(zone.name()) << std::endl;
                    uffile << "+--------------+---------------------------------------------------+" << std::endl;
                    uffile << "| Nuclide Name | Final Concentration                               |" << std::endl;
                    uffile << "+--------------+---------------------------------------------------+" << std::endl;

                    if (__vbs__) std::cout << "Processing zone '" << zone.attribute("name").value() << "'" << std::endl;

                    // Reads the species names.
                    if (!zone.child("species")) throw (int)errex::MISSING_SPECIES_NAMES;
                    std::vector<std::string> species_names;
                    std::string species = zone.child("species").child_value();
                    if (strlen(zone.child("species").attribute("amin").value()) > 0) {
                        xml_document source;
                        xml_parse_result open_success = source.load_file(zone.child("species").attribute("source").value());
                        if (open_success)
                        {

                            if (strlen(zone.child("species").attribute("amin").value()) > 0)
                                AMin = std::stoi(zone.child("species").attribute("amin").value());
                            else
                            {
                                std::cout << "WARNING\t<cnuctran::simulation::from_input()> AMin attribute is missing. Considering all nuclides with A > 0." << std::endl;
                                AMin = 0;
                            }
                            if (strlen(zone.child("species").attribute("amax").value()) > 0)
                                AMax = std::stoi(zone.child("species").attribute("amax").value());
                            else
                            {
                                std::cout << "WARNING\t<cnuctran::simulation::from_input()> AMin attribute is missing. Considering all nuclides with A > 0." << std::endl;
                                AMax = 400;
                            }
                            species_names = get_nuclide_names(zone.child("species").attribute("source").value(), AMin, AMax);
                            if (__vbs__) std::cout << "Building chains... A = [" << AMin << "," << AMax << "]. Total no. of nuclides = " << species_names.size() << std::endl;
                        }
                        else
                            throw (int)errex::MISSING_SPECIES_NAMES;
                    }
                    else
                    {
                        std::stringstream ss = std::stringstream(species); std::string token;
                        while (std::getline(ss, token, ' '))
                        {
                            if (token != "" || token != NULL || strlen(token.c_str()) > 2)
                                species_names.push_back(beautify(token));
                        }
                    }

//..................Reads the initial concentrations for each zone.
                    std::map<std::string, mpreal> w0;
                    const char* w0_source = zone.child("initial_concentrations").attribute("source").value();
                    const char* override_species_names = zone.child("initial_concentrations").attribute("override_species_names").value();
                    if (w0_source != "")
                    {
                        xml_document w0_doc;
                        xml_parse_result load_success = w0_doc.load_file(w0_source);
                        if (override_species_names == "true") species_names.clear();
                        if (load_success)
                        {
                            for (xml_node concs : w0_doc.child("output").children())
                            {
                         
                                if (std::string(concs.name()) != "species_concentrations") continue;
                                if (std::string(concs.attribute("zone").value()) != std::string(zone.attribute("name").value())) continue;
                                if (__vbs__) std::cout << "Reading the initial nuclide concentrations from " << w0_source << " for zone '" << concs.attribute("zone").value() << "'." << std::endl;
                                for (xml_node nuclide : concs)
                                {
                                    if (std::string(nuclide.name()) != "concentration") continue;
                                    std::string name = nuclide.attribute("species").value();
                                    mpreal concentration = mpreal(nuclide.attribute("value").value());
                                    if (override_species_names == "true")
                                    {
                                        species_names.push_back(name);
                                        if (concentration != mpreal("0"))
                                            w0[name] = concentration;
                                    }
                                    else
                                    {
                                        if (find(species_names.begin(), species_names.end(), name) != species_names.end())
                                            w0[name] = concentration;
                                    }
                                }
                            }
                        }
                        else
                        {
                            throw (int)errex::MISSING_W0_SOURCE;
                        }
                    }
                    else
                    {

                        if (!zone.child("initial_concentrations").child_value()) throw (int)errex::MISSING_W0;
                        for (xml_node item : zone.child("initial_concentrations").children())
                        {
                            if (std::string(item.name()) != "concentration") continue;
                            mpreal concentration = mpreal(item.attribute("value").value());
                            w0[item.attribute("species").value()] = concentration;
                            if (__vbs__ == 2) std::cout << "INPUT\t<initial_concentrations> species = " << item.attribute("species").value() << " w0 = " << concentration << std::endl;
                        }
                    }

//..................INITIALIZES the solver.
                    solver sol = solver(species_names);

//..................Reads the user defined removals.
                    std::string daughter_species, parent_species; 
                    for (xml_node removal : zone.child("removals").children())
                    {
                        std::vector<int> daughters;
                        std::vector<mpreal> yields;
                        mpreal rate = mpreal(removal.attribute("rate").value());
                        parent_species = removal.attribute("parent").value();
                        // Reads the product(s)/daughters of the removal.
                        std::stringstream ss(removal.attribute("daughters").value()); std::string token;
                        std::stringstream sy(removal.attribute("yields").value()); std::string ytoken;

                        while (std::getline(ss, token, ' '))
                        {
                            std::getline(sy, ytoken, ' ');
                            if (token != "") 
                            {
                                std::vector<std::string>::iterator it_daughter = find(sol.species_names.begin(), sol.species_names.end(), token);
                                if (it_daughter != species_names.end())
                                {
                                    int daughter_id = distance(sol.species_names.begin(), it_daughter);
                                    daughters.push_back(daughter_id);
                                    if (beautify(ytoken) != NULL || beautify(ytoken) != "")
                                        yields.push_back(mpreal(beautify(ytoken)));
                                    std::vector<std::string>::iterator it_parent = find(sol.species_names.begin(), sol.species_names.end(), parent_species);
                                    if (it_parent != species_names.end())
                                    {
                                       int parent_id = distance(sol.species_names.begin(), it_parent);
                                        sol.add_removal(parent_id, rate, daughters, yields);
                                    }
                                    else
                                    {
                                        std::cout << "WARNING\t<cnuctran::simulation::from_input()> The parent species \"" << parent_species << "\" for the specified removal is not defined in the current calculation.\n";
                                    }
                                    
                                }
                            }
                        }
                        if (__vbs__ == 2) std::cout << "INPUT\t<removal> parent = " << parent_species << " rate = " << rate << std::endl;
                    }
                   
                    


//..................Reads the rxn rates.
                    std::map<std::string, std::map<std::string, mpreal>> rxn_rates;

                    for (xml_node reaction : zone.child("reaction_rates").children())
                    {
                        if (std::string(reaction.name()) != "reaction") continue;
                        mpreal rate = mpreal(reaction.attribute("rate").value());
                        rxn_rates[reaction.attribute("species").value()][reaction.attribute("type").value()] = rate;
                        if (__vbs__ == 2) std::cout << "INPUT\t<rxn_rates> species = " << reaction.attribute("species").value() << " type = " << reaction.attribute("type").value() << " rate = " << rate << std::endl;
                    }

//..................Runs the simulation.

                    build_chains(sol, rxn_rates, zone.child("species").attribute("source").value());
                    std::map<std::string, mpreal> w;
                    w = sol.solve(w0, n, t);
//..................Prints to output file.
                    std::stringstream ss("");
                    ss << "\t<!-- NOTE: All species concentrations less than the epsilon (" << __eps__ << ") \n\t     are not included in this file. -->\n";
                    ss << "\t<species_concentrations zone=\"" << zone.attribute("name").value() << "\" amin = \"" << AMin << "\" amax=\"" << AMax << "\" n=\"" << sol.species_names.size() << "\" time_step=\"" << t << "\">" << std::endl;
                    for (std::string species : sol.species_names)
                    {
                        mpreal c = w[species];
                        if (c > __eps__)
                        {
                            ss << "\t\t<concentration species=\"" << species << "\" value=\"" << std::scientific << std::setprecision(output_digits) << w[species] << "\"/>" << std::endl;
                            uffile << "| " << std::setw(13) << std::left << species << "| " << std::scientific << std::setw(50) << std::left << std::setprecision(output_digits) << w[species] << "|" << std::endl;
                        }
                    }
                    ss << "\t</species_concentrations>" << std::endl;
                    file << ss.str();
                    uffile << "+------------------------------------------------------------------+" << std::endl << std::endl;
                }

                //Closes the output file stream.
                file << "</output>" << std::endl;
                file.close();
                uffile.close();


            }
            catch (int e)
            {
                switch (e)
                {
                case errex::MISSING_SUBSTEP_SIZE:
                    std::cout << "FATAL-ERROR\t<cnuctran::simulation::from_input()> Precision order, n, is not supplied." << std::endl;
                    exit(1);
                case errex::MISSING_STEP_SIZE:
                    std::cout << "FATAL-ERROR\t<cnuctran::simulation::from_input()> Time step, t, is not supplied." << std::endl;
                    exit(1);
                case errex::MISSING_SPECIES_NAMES:
                    std::cout << "FATAL-ERROR\t<cnuctran::simulation::from_input()> Species names are not supplied." << std::endl;
                    exit(1);
                case errex::MISSING_W0:
                    std::cout << "FATAL-ERROR\t<cnuctran::simulation::from_input()> Initial nuclide concentrations are not supplied." << std::endl;
                    exit(1);
                case errex::NUCLIDES_DATA_LOAD_FAILED:
                    std::cout << "FATAL-ERROR\t<cnuctran::simulation::from_input()> Could open source file." << std::endl;
                    exit(1);
                case errex::XML_READING_ERROR:
                    std::cout << "FATAL-ERROR\t<cnuctran::simulation::from_input()> An error has occurred while reading the XML input file: " << xml_input << std::endl;
                    exit(1);
                case errex::MISSING_W0_SOURCE:
                    std::cout << "FATAL-ERROR\t<cnuctran::simulation::from_input()> Could not open the initial nuclide concentrations XML file." << std::endl;
                    exit(1);

                default:
                    std::cout << "FATAL-ERROR\t<cnuctran::simulation::from_input()> Unexpected error has occurred." << std::endl;
                    exit(1);
                }
            }
        }
    };
}


#endif
