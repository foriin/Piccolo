#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>

struct Cluster {
    std::string chromosome;
    int start;
    int end;
    double score;
    char strand;
    int count;
};

bool overlaps(const Cluster& a, const Cluster& b) {
    return std::max(a.start, b.start) < std::min(a.end, b.end);
}

bool sameChromosome(const Cluster& a, const Cluster& b) {
    return a.chromosome == b.chromosome;
}

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

    std::vector<Cluster> clusters;
    Cluster temp;
    while (infile >> temp.chromosome >> temp.start >> temp.end >> temp.score >> temp.strand) {
        temp.count = 1;
        clusters.push_back(temp);
    }

    std::sort(clusters.begin(), clusters.end(), [](const Cluster& a, const Cluster& b) {
        return a.chromosome < b.chromosome || (a.chromosome == b.chromosome && a.start < b.start);
    });

    std::vector<Cluster> merged;
    for (const auto& cluster : clusters) {
        bool foundOverlap = false;
        for (auto& m : merged) {
            if (sameChromosome(cluster, m) && overlaps(cluster, m) && cluster.strand != m.strand) {
                m.start = std::min(m.start, cluster.start);
                m.end = std::max(m.end, cluster.end);
                m.score = (m.score * m.count + cluster.score) / (m.count + 1);
                m.count += 1;
                m.strand = '.';
                foundOverlap = true;
                break;
            }
        }
        if (!foundOverlap) {
            merged.push_back(cluster);
        }
    }

    int rowNumber = 0;
    for (const auto& m : merged) {
        std::string clusterName = "cluster_" + std::to_string(++rowNumber);
        outfile << m.chromosome << "\t" << m.start << "\t" << m.end << "\t"
                << clusterName << "\t" << m.score << "\t"
                << (m.strand != '.' ? std::string(1, m.strand) : ".") << "\n";
    }

    infile.close();
    outfile.close();

    return 0;
}

