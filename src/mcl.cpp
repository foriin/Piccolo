#include <iostream>
#include <fstream>
#include <string>

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input_bed_file> <output_bed_file>\n";
        return 1;
    }

    std::ifstream infile(argv[1]);
    if (!infile.is_open()) {
        std::cerr << "Error opening input file\n";
        return 1;
    }

    std::ofstream outfile(argv[2]);
    if (!outfile.is_open()) {
        std::cerr << "Error opening output file\n";
        return 1;
    }

    std::string line;
    while (std::getline(infile, line)) {
        outfile << line << "\thaha\n";
    }

    infile.close();
    outfile.close();

    return 0;
}

