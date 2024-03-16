#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <iostream>

void PrintVectorSizes(const char* fileName = "itsdigits.root") {
    TFile* file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Could not open the file or the file is not valid." << std::endl;
        return;
    }

    TTree* tree = nullptr;
    file->GetObject("o2sim", tree);  // Adjust the tree name if needed
    if (!tree) {
        std::cerr << "Could not retrieve the TTree." << std::endl;
        file->Close();
        delete file;
        return;
    }

    TObjArray* branches = tree->GetListOfBranches();
    for (int i = 0; i < branches->GetEntries(); ++i) {
        TBranch* branch = (TBranch*)branches->At(i);
        TLeaf* leaf = branch->GetLeaf(branch->GetName());

        if (leaf) {
            // Read the first entry of the branch to populate the leaf values
            branch->GetEntry(0);

            // For vector-like structures, the leaf count should reflect the number of elements
            TLeaf* leafCount = leaf->GetLeafCount();
            Int_t count = leafCount ? leafCount->GetMaximum() : 1;

            if (count > 1 || leafCount) {  // Check if it's likely a vector
                std::cout << "Branch: " << branch->GetName() << ", Likely Vector Size: " << count << std::endl;
            } else {
                // If count is 1 and there's no leaf count, it's likely not a vector or a single-element vector
                std::cout << "Branch: " << branch->GetName() << ", Not a vector or Single-element vector" << std::endl;
            }
        }
    }

    file->Close();
    delete file;
}

