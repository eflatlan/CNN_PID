#include <TFile.h>
#include <TKey.h>
#include <TList.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <iostream>

void PrintFileContent(const char* fileName = "itsdigits.root") {
    TFile* file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Could not open the file or the file is not valid." << std::endl;
        return;
    }

    TList* list = file->GetListOfKeys();
    TIter iter(list);
    TKey* key;

    while ((key = (TKey*)iter())) {
        if (strcmp(key->GetClassName(), "TTree") == 0) {
            TTree* tree = (TTree*)file->Get(key->GetName());
            if (!tree) {
                std::cerr << "Could not retrieve the TTree." << std::endl;
                continue;
            }

            std::cout << "Tree found: " << tree->GetName() << std::endl;
            std::cout << "Number of entries: " << tree->GetEntries() << std::endl;

            TObjArray* branchList = tree->GetListOfBranches();
            if (branchList) {
                TIter branchIter(branchList);
                TBranch* branch;

                while ((branch = (TBranch*)branchIter())) {
                    std::string branchName = branch->GetName();
                    std::string branchType = branch->GetClassName() ? branch->GetClassName() : "";

                    std::cout << "  Branch: " << branchName << ", Type: " << branchType;

                    if (branchType.find("vector") != std::string::npos) {
                        // Attempt to get the size of the vector
                        TLeaf* leaf = branch->GetLeaf(branchName.c_str());
                        if (leaf) {
                            tree->GetEntry(0); // Load the first entry
                            std::cout << ", Vector Size: " << leaf->GetLen();
                        } else {
                            std::cout << " (No leaf available to determine size)";
                        }
                    }
                    std::cout << std::endl;
                }
            }
            std::cout << std::endl; // Add a blank line for readability
        }
    }

    file->Close();
    delete file;
}

