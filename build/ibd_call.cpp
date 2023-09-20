#include "SyllableQuery.cpp"
#include <math.h>
#include <iostream>
#include <string>
#include <ctime>
#include <filesystem>

void call_ibds(char *input_ref_file, std::string query_file, std::string output, const double L, const int B, const std::string chrom) {
    // char *wdc = getcwd(NULL, 0);
    // string swd = string(wdc) + "/";
    // free(wdc);

    std::string out_dir = output + "/sq_ibds" ;
    std::error_code err;
    if (!std::filesystem::exists(out_dir)) {
        if (!std::filesystem::create_directories(out_dir, err)) {
            std::cout << "CreateDirectoryRecuresive: FAILED to create " + out_dir + ", err: " + err.message() + "\n";
            exit(1);
        }
    }


    if (B == 64) {
        SyllableQuery<unsigned long long> sq;

        // compress and make the ref panel ready
        int code = sq.precompute(input_ref_file);
        if (code == 1) {
            std::cerr << "Input file does not exist or cannot be read.\n";
            return;
        }

        if (code == 2) {
            std::cerr << "Input file not in VCF format.\n";
            return;
        }
        cout << "Precomputation complete" << endl;

        // make the queries
        std::string output_file = out_dir + "/" + chrom + "-sq-ibd-results.txt";
        int q_code;
        const clock_t START = clock();
        q_code = sq.query(query_file.c_str(), output_file.c_str(), (int) ceil(L));

        if (q_code == 1) {
            std::cout << "Query file does not exist or cannot be read.\n";
            exit(1);
        } else if (q_code == 2) {
            std::cout <<  "Output file " + output_file + " cannot be written to.\n";
            exit(1);
        } else if (q_code == 4) {
            std::cout << "The query length, " + to_string(L) + " sites, is less than the minimum allowable query length, " + 
            to_string(sq.minSiteL) + " sites.\n";
            exit(1);
        } else if (q_code == 5) {
            std::cout << "Query file not in VCF format or has a different number of sites from the initial panel.\n";
            exit(1);
        } else {
            std::cout << "Query results outputted to " + output_file + ": Elapsed CPU time (s): " + 
            to_string(double(clock() - START) / CLOCKS_PER_SEC) << "\n";
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
        std::string output_file = out_dir + "/" + chrom + "-sq-ibd-results.txt";
        int q_code;
        const clock_t START = clock();
        q_code = sq.query(query_file.c_str(), output_file.c_str(), (int) ceil(L));

        if (q_code == 1) {
            std::cout << "Query file does not exist or cannot be read.\n";
            exit(1);
        } else if (q_code == 2) {
            std::cout <<  "Output file " + output_file + " cannot be written to.\n";
            exit(1);
        } else if (q_code == 4) {
            std::cout << "The query length, " + to_string(L) + " sites, is less than the minimum allowable query length, " + 
            to_string(sq.minSiteL) + " sites.\n";
        } else if (q_code == 5) {
            std::cout << "Query file not in VCF format or has a different number of sites from the initial panel.\n";
            exit(1);
        } else {
            std::cout << "Query results outputted to " + output_file + ": Elapsed CPU time (s): " + 
            to_string(double(clock() - START) / CLOCKS_PER_SEC) << "\n";
        }
    }
}
