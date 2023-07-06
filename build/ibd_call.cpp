#include "SyllableQuery.cpp"
#include <math.h>
#include <iostream>
#include <string>
#include <ctime>

void call_ibds(char *input_ref_file, std::string query_file, const double L, const int B) {
    char *wdc = getcwd(NULL, 0);
    string swd = string(wdc) + "/";
    free(wdc);

    if (B == 64) {
        SyllableQuery<unsigned long long> sq;

        // compress and make the ref panel ready
        int code = sq.precompute(input_ref_file);
        if (code == 1) {
            cerr << "Input file does not exist or cannot be read." << endl;
            return;
        }

        if (code == 2) {
            cerr << "Input file not in VCF format." << endl;
            return;
        }
        cout << "Precomputation complete" << endl;

        // make the queries
        std::string output_file = swd + "sq-results.txt";
        int q_code;
        const clock_t START = clock();
        q_code = sq.query(query_file.c_str(), output_file.c_str(), (int) ceil(L));

        if (q_code == 1) {
            std::cout << "Query file does not exist or cannot be read.\n";
            exit(1);
        } else if (code == 4) {
            std::cout << "The query length, " + to_string(L) + " sites, is less than the minimum allowable query length, " + 
            to_string(sq.minSiteL) + " sites.\n";
        } else {
            std::cout << "Query results outputted to " + output_file + "\nElapsed CPU time (s): " + 
            to_string(double(clock() - START) / CLOCKS_PER_SEC);
        }

    } else if (B == 128) {
        SyllableQuery<unsigned __int128> sq;

        // compress and make the ref panel ready
        int code = sq.precompute(input_ref_file);
        if (code == 1) {
            cerr << "Input file does not exist or cannot be read." << endl;
            return;
        }

        if (code == 2) {
            cerr << "Input file not in VCF format." << endl;
            return;
        }
        cout << "Precomputation complete" << endl;

        // make the queries
        std::string output_file = swd + "sq-results.txt";
        int q_code;
        const clock_t START = clock();
        q_code = sq.query(query_file.c_str(), output_file.c_str(), (int) ceil(L));

        if (q_code == 1) {
            std::cout << "Query file does not exist or cannot be read.\n";
            exit(1);
        } else if (code == 4) {
            std::cout << "The query length, " + to_string(L) + " sites, is less than the minimum allowable query length, " + 
            to_string(sq.minSiteL) + " sites.\n";
        } else {
            std::cout << "Query results outputted to " + output_file + "\nElapsed CPU time (s): " + 
            to_string(double(clock() - START) / CLOCKS_PER_SEC);
        }
    }
}