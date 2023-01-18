/*

      This file is part of the CNUCTRAN library

      @author   M. R. Omar (rabieomar@usm.my)
      @license  MIT
      @link     https://github.com/rabieomar92/cnuctran

      Copyright (c) 2023, Universiti Sains Malaysia

      This header file contains the definitions of all methods for preparing
      the simulation associated to the proposed probabilistic method.

 */

#ifndef SIMULATION_H
#define SIMULATION_H

#include <iomanip>
#include <fstream>
#include <pugixml.hpp>
#include <solver.h>

using namespace pugi;
using namespace mpfr;
using namespace std;

namespace cnuctran {


    class simulation
    {

    public:

        /*
            This sub-routine reads nuclide transmutation information from the XML 
            nuclides data library.

        */
        static void build_chains(solver& s, map<string, map<string, mpreal>>& rxn_rates,
            string xml_data_location)
        {
            vector<string> species_names = s.species_names;

            //..........Loads the nuclides data from the XML source file. 
            xml_document file;
            xml_parse_result open_success = file.load_file(xml_data_location.c_str());
            if (!open_success)
            {
                cout << "INFO\t<cnuctran::depletion_scheme::build_chains(...)> Nuclides data file is not provided." << endl;
                return;
            }

            //..........Unpack the XML root.
            xml_node root = file.child("depletion");

            for (xml_node species : root.children())
            {
                string species_name = species.attribute("name").value();
                vector<string>::iterator it = std::find(species_names.begin(), species_names.end(), species_name);
                if (it == species_names.end())
                    continue;


                mpreal decay_rate;
                if (species.attribute("half_life"))
                    decay_rate = mpfr::log(mpreal("2")) / mpreal(species.attribute("half_life").value());
                else
                    decay_rate = mpreal("0");

                for (xml_node removal : species.children())
                {
                    if (string(removal.name()) == "decay_type")
                    {
                        mpreal decay_rate_adjusted = mpreal(removal.attribute("branching_ratio").value()) * decay_rate;
                        string parent = species_name;
                        string daughter = removal.attribute("target").value();
                        vector<string>::iterator it_parent = std::find(species_names.begin(), species_names.end(), parent);
                        int parent_id = distance(species_names.begin(), it_parent);
                        vector<string>::iterator it_daughter = std::find(species_names.begin(), species_names.end(), daughter);
                        if (it_daughter != species_names.end())
                        {
                            int daughter_id = distance(species_names.begin(), it_daughter);
                            s.add_removal(parent_id, decay_rate_adjusted, vector<int>({ daughter_id }));
                        }
                        else
                        {
                            s.add_removal(parent_id, decay_rate_adjusted, vector<int>({ __nop__ }));
                        }

                    }

                    if (rxn_rates.size() != 0)
                    {
                        if (rxn_rates.count(species_name))
                        {
                            if (string(removal.name()) == "reaction_type" && removal.attribute("target"))
                            {
                                string parent = species_name;
                                vector<string>::iterator it_parent = std::find(species_names.begin(), species_names.end(), parent);
                                int parent_id = distance(species_names.begin(), it_parent);
                                if (rxn_rates[parent].count(removal.attribute("type").value()) &&
                                    removal.attribute("type").value() != "fission")
                                {
                                    string daughter = removal.attribute("target").value();
                                    mpreal removal_rate = mpreal(rxn_rates[parent][removal.attribute("type").value()]);
                                    vector<string>::iterator it_daughter = std::find(species_names.begin(), species_names.end(), daughter);
                                    if (it_daughter != species_names.end())
                                    {
                                        int daughter_id = distance(species_names.begin(), it_daughter);
                                        s.add_removal(parent_id, removal_rate, vector<int>({ daughter_id }));
                                    }
                                    else
                                    {
                                        s.add_removal(parent_id, removal_rate, vector<int>({ __nop__ }));
                                    }
                                }
                            }

                            if (string(removal.name()) == "neutron_fission_yields")
                            {
                                string parent = species_name;
                                vector<string>::iterator it_parent = std::find(species_names.begin(), species_names.end(), parent);
                                int parent_id = distance(species_names.begin(), it_parent);
                                mpreal energy = mpreal("0");
                                vector<string> products;
                                vector<mpreal> yields;
                                if (rxn_rates[parent].count("fission") != 0)
                                {
                                    for (xml_node data : removal.children())
                                    {
                                        if (string(data.name()) == "energies")
                                        {
                                            stringstream ss(data.child_value()); string token;
                                            vector<string> energies;
                                            while (getline(ss, token, ' '))
                                                if (token != "")
                                                    energies.push_back(token);
                                            energy = energies[0];
                                        }

                                        if (string(data.name()) == "fission_yields")
                                        {
                                            if (mpreal(data.attribute("energy").value()) == energy)
                                            {
                                                for (xml_node param : data.children())
                                                {
                                                    if (string(param.name()) == "products")
                                                    {
                                                        stringstream ss(param.child_value()); string token;
                                                        while (getline(ss, token, ' '))
                                                            if (token != "")
                                                                products.push_back(token);
                                                    }

                                                    if (string(param.name()) == "data")
                                                    {
                                                        stringstream ss(param.child_value()); string token;
                                                        while (getline(ss, token, ' '))
                                                            if (token != "")
                                                                yields.push_back(mpreal(token));
                                                    }

                                                }

                                                mpreal total_fission_rate = rxn_rates[parent]["fission"];
                                                vector<mpreal> yields_to_add;
                                                vector<int> daughters_id_to_add;
                                                vector<string>::iterator it_product;
                                                for (string product : products)
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
        static vector<string> get_nuclide_names(string xml_data_location, int AMin = -1, int AMax = -1)
        {
            xml_document file;
            xml_parse_result open_success = file.load_file(xml_data_location.c_str());
            if (!open_success)
            {
                cout << "ERROR <cnuctran.depletion_scheme.get_nuclide_names(...)>\nFail retrieving data from " << xml_data_location << "." << endl;
                return vector<string>();
            }


            xml_node root = file.child("depletion");

            vector<string> species_names;

            for (xml_node species : root.children())
            {
                string name = species.attribute("name").value();
                stringstream ss(name); string token;
                getline(ss, token, '_');
                string x = "";
                for (char c : token)
                    if (isdigit(c)) x += c;
                if (AMin == AMax == -1)
                    species_names.push_back(name);
                else
                {
                    int A = stoi(x);
                    if (A >= AMin && A <= AMax)
                        species_names.push_back(name);
                }
            }
            return species_names;
        }


        /*
            Reads the input XML file (input.xml) and obtains all simulation parameters. Finally, this
            routine runs the simulation.

            P.S. Many apologies for the complicated code.
        */
        static void from_input(string xml_input)
        {

            try
            {
                //..............Loads the XML input file. Handle loading errors.
                xml_document input_file;
                xml_parse_result open_success = input_file.load_file(xml_input.c_str());
                if (!open_success) {
                    cout << "fatal-error <cnuctran.simulation.from_input()>\n" << open_success.description() << endl;
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
                string output_location = "output.xml";

                //Obtains the verbosity level from the input file.
                tmp = root.child("simulation_params").child("verbosity").child_value();
                tmp != "" ? __vbs__ = stoi(tmp) : __vbs__ = 0;

                //Obtains the precision digits from the input file.
                tmp = root.child("simulation_params").child("precision_digits").child_value();
                if (tmp != "")
                {
                    precision_digits = stoi(tmp);
                    __dps__ = precision_digits;
                } 
                else
                    precision_digits = __dps__;
              
                if (precision_digits < 30)
                {
                    precision_digits = 30;
                    cout << "warning <cnuctran.simulation.from_input()>\nA precision < 30 digits is vulnerable to errorneous arithmetics that lead to fatal error." << endl;
                    cout << "The precision was set to 30 digits." << endl;
                }

                //Obtains the minimum rate from the input file.
                tmp = root.child("simulation_params").child("min_rate").child_value();
                if (tmp != "") __mnr__ = mpreal(tmp).toDouble();

                //Obtains the minimum rate from the input file.
                tmp = root.child("simulation_params").child("max_rate").child_value();
                if (tmp != "") __mxr__ = mpreal(tmp).toDouble();

                //Obtains the epsilon from the input file.
                tmp = root.child("simulation_params").child("epsilon").child_value();
                if (tmp != "") __eps__ = mpreal(tmp);

                //Obtains the output precision digits from the input file.
                tmp = root.child("simulation_params").child("output_digits").child_value();
                tmp != "" ? output_digits = stoi(tmp) : output_digits = __dop__;


                //IMPORTANT! Reads the precision before declaring any high-precision float vars.
                mpreal::set_default_prec(digits2bits(precision_digits));
                cout.precision(output_digits);

                //Reads the substep interval in seconds.
                tmp = root.child("simulation_params").child("n").child_value();
                tmp != "" ? n = mpreal(tmp) : throw (int)errex::MISSING_SUBSTEP_SIZE;

                //Reads the time step in seconds.
                tmp = root.child("simulation_params").child("time_step").child_value();
                tmp != "" ? t = mpreal(tmp) : throw (int)errex::MISSING_STEP_SIZE;

                //Reads the final councentration output file location.
                tmp = root.child("simulation_params").child("output").child_value();
                tmp != "" ? output_location = tmp : output_location = ".//output.xml";
                if (__vbs__) cout << "Final nuclide concentrations will be written in " << output_location << endl;

//..............LOOP OVER ALL ZONES and read the species, initial concentrations and reaction rates (for the zone).
                int nzones = 0;
                //Initialize file stream of the output file (output.xml).
                ofstream file_xml;
                ofstream file_out;
                file_xml.open(output_location, ios::out);
                file_out.open(output_location + ".out", ios::out);
                file_xml << "<output>" << endl;
                
                //Loop over all relevant XML child nodes for each zone.
                for (xml_node zone : root.children())
                {
                    // Filters XML nodes that are not a zone node.
                    if (string(zone.name()) != "zone") continue;
                    if (__vbs__) cout << "Processing zone '" << zone.attribute("name").value() << "'" << endl;

                    // Reads the species names.
                    if (!zone.child("species")) throw (int)errex::MISSING_SPECIES_NAMES;
                    vector<string> species_names;
                    auto species = zone.child("species").child_value();
                    if (strlen(zone.child("species").attribute("amin").value()) > 0) {
                        xml_document source;
                        auto open_success = source.load_file(zone.child("species").attribute("source").value());
                        if (open_success)
                        {   
                            AMax = strlen(zone.child("species").attribute("amax").value()) > 0 ? stoi(zone.child("species").attribute("amax").value()) : 400;
                            AMin = strlen(zone.child("species").attribute("amin").value()) > 0 ? stoi(zone.child("species").attribute("amin").value()) : 0;
                            species_names = get_nuclide_names(zone.child("species").attribute("source").value(), AMin, AMax);
                            if (__vbs__) cout << "Building chains... A = [" << AMin << "," << AMax << "]. Total no. of nuclides = " << species_names.size() << endl;
                        }
                        else
                            throw (int)errex::MISSING_SPECIES_NAMES;
                    }
                    else
                    {
                        
                        stringstream ss = stringstream(species); string token;
                        while (getline(ss, token, ' '))
                        {
                            token.erase(remove_if(token.begin(), token.end(), ispunct), token.end());
                            if (token != "") species_names.push_back(trim(token));
                        }
                        if (__vbs__) cout << "Building chains... Total no. of nuclides = " << species_names.size() << endl;
                    }

//..................Reads the initial concentrations for each zone.
                    map<string, mpreal> w0;
                    const char* w0_source = zone.child("initial_concentrations").attribute("source").value();
                    const char* override_species_names = zone.child("initial_concentrations").attribute("override").value();
                    if (w0_source != "")
                    {
                        xml_document w0_doc;
                        auto load_success = w0_doc.load_file(w0_source);
                        if (override_species_names == "true") species_names.clear();
                        if (load_success)
                        {
                            for (xml_node concs : w0_doc.child("output").children())
                            {

                                if (string(concs.name()) != "nuclide_concentrations") continue;
                                if (string(concs.attribute("zone").value()) != string(zone.attribute("name").value())) continue;
                                if (__vbs__) cout << "Reading the initial nuclide concentrations from " << w0_source << " for zone '" << concs.attribute("zone").value() << "'." << endl;
                                for (xml_node nuclide : concs)
                                {
                                    if (string(nuclide.name()) != "concentration") continue;
                                    auto name = nuclide.attribute("species").value();
                                    auto concentration = mpreal(nuclide.attribute("value").value());
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
                            if (string(item.name()) != "concentration") continue;
                            mpreal concentration = mpreal(item.attribute("value").value());
                            w0[item.attribute("species").value()] = concentration;
                        }
                    }

//..................Reads the rxn rates.
                    map<string, map<string, mpreal>> rxn_rates;

                    for (xml_node reaction : zone.child("reaction_rates").children())
                    {
                        mpreal rate = mpreal(reaction.attribute("rate").value());
                        rxn_rates[reaction.attribute("species").value()][reaction.attribute("type").value()] = rate; 
                    }

//ATTENTION!........This is where the code solves for nuclides concentrations.
                    solver sol = solver(species_names);
                    build_chains(sol, rxn_rates, zone.child("species").attribute("source").value());
                    auto w = sol.solve(w0, n, t);

//..................Prints to output file.
                    stringstream ss_xml("");
                    stringstream ss_out("");
                    ss_xml << "\t<nuclide_concentrations zone=\"" 
                           << zone.attribute("name").value() 
                           << "\" amin = \"" << AMin 
                           << "\" amax=\"" << AMax 
                           << "\" total_nuclides=\"" << sol.species_names.size() 
                           << "\" time_step=\"" << t << "\">" << endl;

                    file_out << "CNUCTRAN v1.1 OUTPUT." << endl;
                    file_out << setw(20) << left << "time step" << "= " << scientific << t << "s" << endl;
                    file_out << setw(20) << left << "order (n)" << "= " << (int)n.toFloat() << endl;
                    file_out << setw(20) << left << "total nuclides" << "= " << w.size() << endl;
                    int k = int(floor(log(t / pow(mpreal("10"), -n)) / log(mpreal("2.0"))));
                    mpreal dt = t / pow(mpreal("2.0"), k);
                    file_out << setw(20) << left << "substep" << "= " << scientific << dt << "s (" << k << " sparse mults.)" << endl;
                    file_out << setw(20) << left << "precision" << "= " << precision_digits << " digits." << endl;
                    file_out << setw(8) << left << "Species" << setw(10) << left << "Non-zero" << setw(output_digits + 10) << left << "Concentration" << endl;

                    for (string species : sol.species_names)
                    {
                        mpreal c = w[species];
                        ss_xml << "\t\t<concentration species=\"" << species 
                            << "\" value=\"" << scientific 
                            << setprecision(output_digits) << c 
                            << "\" />" << endl;
                        ss_out << setw(8) << left << species << " " << setw(10);
                        w[species] > mpreal("0.0") ? ss_out << left << "yes" : ss_out << "";
                        ss_out << setw(output_digits + 10) << scientific
                               << setprecision(output_digits) << left << c << endl;

                    }
                    ss_xml << "\t</nuclide_concentrations>" << endl;
                    file_out << ss_out.str();
                    file_xml << ss_xml.str();

                }

                //Closes the output file stream.
                file_xml << "</output>" << endl;
                file_xml.close();
                file_out.close();


            }
            catch (int e)
            {
                switch (e)
                {
                case errex::MISSING_SUBSTEP_SIZE:
                    cout << "fatal-error <cnuctran.simulation.from_input()>\nPrecision order, n, is not supplied." << endl;
                    exit(1);
                case errex::MISSING_STEP_SIZE:
                    cout << "fatal-error <cnuctran.simulation.from_input()>\nTime step, t, is not supplied." << endl;
                    exit(1);
                case errex::MISSING_SPECIES_NAMES:
                    cout << "fatal-error <cnuctran.simulation.from_input()>\nSpecies names are not supplied." << endl;
                    exit(1);
                case errex::MISSING_W0:
                    cout << "fatal-error <cnuctran.simulation.from_input()>\nInitial nuclide concentrations are not supplied." << endl;
                    exit(1);
                case errex::NUCLIDES_DATA_LOAD_FAILED:
                    cout << "fatal-error <cnuctran.simulation.from_input()>\nCould open source file." << endl;
                    exit(1);
                case errex::XML_READING_ERROR:
                    cout << "fatal-error <cnuctran.simulation.from_input()>\nAn error has occurred while reading the XML input file: " << xml_input << endl;
                    exit(1);
                case errex::MISSING_W0_SOURCE:
                    cout << "fatal-error <cnuctran.simulation.from_input()>\nCould not open the initial nuclide concentrations XML file." << endl;
                    exit(1);
                default:
                    cout << "fatal-error <cnuctran.simulation.from_input()>\nUnexpected error has occurred." << endl;
                    exit(1);
                }
            }
        }
    };
}


#endif
