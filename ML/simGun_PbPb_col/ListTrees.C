#include <TFile.h>
#include <TKey.h>
#include <TList.h>
#include <iostream>

void ListTrees(const char* fileName) {
    TFile* file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Could not open the file or the file is not valid." << std::endl;
        return;
    }

    TList* list = file->GetListOfKeys();
    TIter iter(list);
    TKey* key;

    while ((key = (TKey*)iter())) {
        // Check if the object is a TTree
        if (strcmp(key->GetClassName(), "TTree") == 0) {
            std::cout << "Tree found: " << key->GetName() << std::endl;
        }
    }

    file->Close();
    delete file;
}

