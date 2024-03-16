#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <iostream>

void PrintFileContent2(const char* fileName = "hmpidclusters.root") {
    TFile* file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Could not open the file or the file is not valid." << std::endl;
        return;
    }

    TTree* tree = nullptr;
    file->GetObject("o2hmp", tree); // Adjust the tree name if needed
    if (!tree) {
        std::cerr << "Could not retrieve the TTree." << std::endl;
        file->Close();
        delete file;
        return;
    }

    TTreeReader reader(tree);

    TObjArray* branchList = tree->GetListOfBranches();
    for (int i = 0; i < branchList->GetEntries(); ++i) {
        TBranch* branch = (TBranch*)branchList->At(i);
        TString branchName = branch->GetName();
        TString branchType = branch->GetClassName();
        std::cout << "Branch: " << branchName << ", Type: " << branchType << std::endl;
        

        if (branchType.Contains("vector")) {
            // Use TTreeReaderArray for vectors
            TTreeReaderArray<double> array(reader, branchName);  // Using double as a placeholder type

            reader.SetEntry(0); // Access the first entry to infer the size
            std::cout << "Branch: " << branchName << ", Type: " << branchType << ", Approx. Size: " << array.GetSize() << std::endl;
        } 
    }

    file->Close();
    delete file;
}

