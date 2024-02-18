

#ifndef TEST_POPULATE

#define TEST_POPULATE

// #define _USE_MATH_DEFINES

#include "Sigma.cpp"

#include "Sigma2.cpp"

#include "populate.cpp"

#include "populate2.cpp" // TODO: change name of class and file here

#include <cmath>

#include <iostream>

#include <random>

#include <chrono>

#include <iostream>

#include <thread>

#include <array>

#include <cassert>

#include <cmath>

#include <iostream>

#include <math.h>

#include <vector>

// #include "ReconE.cpp"

// #include "ReconG.cpp"

// namespace ParticleUtils

using namespace o2;

struct Bin
{

    float x;

    float y;
};

class ArrAndMap
{

public:
    using vecArray2 = std::vector<std::array<double, 2>>;

    vecArray2 arrMaxPionPos, arrMaxKaonPos, arrMaxProtonPos;

    vecArray2 arrMinPionPos, arrMinKaonPos, arrMinProtonPos;

    TVector2 *errPos = new TVector2(0, 0);

    void setMaxArrays(const vecArray2 &_arrMaxPionPos,

                      const vecArray2 &_arrMaxKaonPos,

                      const vecArray2 &_arrMaxProtonPos)
    {

        arrMaxPionPos = _arrMaxPionPos;

        arrMaxKaonPos = _arrMaxKaonPos;

        arrMaxProtonPos = _arrMaxProtonPos;
    }

    std::unique_ptr<TH2F> limMin;

    std::unique_ptr<TH2F> limMax;

    void setDrawLimits(const TH2F *limMinIn, const TH2F *limMaxIn)
    {

        //// limMin.reset(new TH2F(*// limMinIn));

        //// limMax.reset(new TH2F(*// limMaxIn));
    }

    void setErrorPos(const TVector2 &posPhoton) { errPos->Set(posPhoton); }

    void setMinArrays(const vecArray2 &_arrMinPionPos,

                      const vecArray2 &_arrMinKaonPos,

                      const vecArray2 &_arrMinProtonPos)
    {

        arrMinPionPos = _arrMinPionPos;

        arrMinKaonPos = _arrMinKaonPos;

        arrMinProtonPos = _arrMinProtonPos;
    }

    void fillMapFromVec(TH2F *map, const vecArray2 &vecArr)
    {

        for (const auto vE : vecArr)
        {

            // Printf("fillMapFromVec El %.2f %.2f", vE[0], vE[1]);

            map->Fill(vE[0], vE[1]);
        }
    }

    int scale = 10;



    void fillmipSizeFilter(double x, double y)
    {

        mipSizeFilter.emplace_back(std::array<double, 2>{x, y});
    }

    void fillckovCandMapOutRange(double x, double y)
    {

        ckovCandMapOutRange.emplace_back(std::array<double, 2>{x, y});
    }

    void fillMipCharge(double x, double y)
    {

        mipChargeFilter.emplace_back(std::array<double, 2>{x, y});
    }

    vecArray2 mipChargeFilter, ckovCandMapRange, mipSizeFilter,

        ckovCandMapOutRange;


    int eventCnt;

    std::unique_ptr<Populate2> populatePtr;

    void setPopulatePtr(std::unique_ptr<Populate2> &&_populatePtr)
    {

        populatePtr = std::move(_populatePtr);
    }

    int eventCount = 0;

    void setEventCount(int _eventCount) { eventCount = _eventCount; }

    void drawTotalMap(std::vector<o2::hmpid::ClusterCandidate> &clusterTrack,

                      int &plotNumber, float xMip, float yMip,

                      vecArray2 pionCandidates, vecArray2 kaonCandidates,

                      vecArray2 protonCandidates, vecArray2 canCombined,

                      int trackPdg, vecArray2 allCand, float sigSep)
    {

        const auto trkRad = populatePtr->getTrackPos();

        const auto trkPC = populatePtr->getPcImp();

        auto xr = trkRad.X();

        auto yr = trkRad.Y();

        TVector2 mip(xMip, yMip);

        auto mipPhi = (mip - trkRad).Phi();

        auto pcPhi = (trkPC - trkRad).Phi();

        // Printf(" mip { %.2f,  %.2f}  trkPC { %.2f,  %.2f} trkRad { %.2f,  %.2f}",

        // mip.X(), mip.Y(), trkPC.X(), trkPC.Y(), trkRad.X(), trkRad.Y());

        // Printf(" (mip-trkRad) %.2f(trkPC-trkRad)  %.2f", (mip-trkRad).Phi(),

        // (trkPC-trkRad).Phi()); Printf(" (mip-trkRad) { %.2f,  %.2f}

        // (trkPC-trkRad) { %.2f,  %.2f}", (mip-trkRad).Phi(), (mip-trkRad).Theta(),

        // (trkPC-trkRad).Phi(), (trkPC-trkRad).Theta());

        auto len1 = 15.;

        auto len2 = 10.;

        std::unique_ptr<TLine> tLineMIP(new TLine(

            xr - len2 * TMath::Cos(mipPhi), yr - len2 * TMath::Sin(mipPhi),

            xr + len1 * TMath::Cos(mipPhi), yr + len1 * TMath::Sin(mipPhi)));

        std::unique_ptr<TLine> tLineTRK(new TLine(

            xr - len2 * TMath::Cos(pcPhi), yr - len2 * TMath::Sin(pcPhi),

            xr + len1 * TMath::Cos(pcPhi), yr + len1 * TMath::Sin(pcPhi)));

        auto distPC2MIP = (mip - trkPC).Mod();

        int numPion = pionCandidates.size();

        int numKaon = kaonCandidates.size();

        int numProton = protonCandidates.size();

        int numTotal = canCombined.size();

        const auto trkdir = populatePtr->getTrkDir();

        auto theta = trkdir.Theta();

        auto st = Form(

            "\#Pi %d K %d Pr %d T %d | theta %.3f pPhi %.3f | MIP %.2f %.2f",

            numPion, numKaon, numProton, numTotal, theta, pcPhi, mip.X(), mip.Y());

        std::unique_ptr<TH2F> hCkovCandMapRange(

            new TH2F(st, st, 1600, 0, 300, 1440, 0, 143));

        std::unique_ptr<TH2F> hmipSizeFilter(

            new TH2F("mipSizeFilter", "mipSizeFilter", 1600, 0, 300, 1440, 0, 143));

        std::unique_ptr<TH2F> hCkovCandMapOutRange(

            new TH2F("ckovCandMapOutRange", "ckovCandMapOutRange", 1600, 0, 300,

                     1440, 0, 143));

        std::unique_ptr<TH2F> hmipChargeFilter(new TH2F(

            "mipChargeFilter", "mipChargeFilter", 1600, 0, 300, 1440, 0, 143));

        std::unique_ptr<TH2F> hSignalAndNoiseMap(

            new TH2F("Signal and Noise ", "Signal and Noise ; x [cm]; y [cm]", 1600,

                     0, 300, 1440, 0, 1433));

        std::unique_ptr<TH2F> hSignalMIP(new TH2F("hmip ", "hmip ; x [cm]; y [cm]",

                                                  1600, 0., 300., 1440, 0, 143));

        std::unique_ptr<TH2F> hSignalMIPpc(new TH2F(

            "hmip pc", "hmip pc; x [cm]; y [cm]", 1600, 0., 300., 1440, 0, 143));

        std::unique_ptr<TH2F> hMaxProton(

            new TH2F("maxPoss Ckov Proton", "maxPoss Ckov Proton; x [cm]; y [cm]",

                     1600, 0, 300, 1440, 0, 143));

        std::unique_ptr<TH2F> hMaxPion(new TH2F("maxPoss Ckov Pion",

                                                "maxPoss Ckov Pion; x [cm]; y [cm]",

                                                1600, 0, 300, 1440, 0, 143));

        std::unique_ptr<TH2F> hMaxPionMinL(new TH2F(

            "maxPoss Ckov Pion min L", "maxPoss Ckov Pion min L; x [cm]; y [cm]",

            1600, 0, 300, 1440, 0, 143));

        std::unique_ptr<TH2F> hMinPionMaxL(new TH2F(

            "minPoss Ckov Pion max L", "minPoss Ckov Pion max L; x [cm]; y [cm]",

            1600, 0, 300, 1440, 0, 143));

        std::unique_ptr<TH2F> hMaxKaon(new TH2F("maxPoss Ckov Kaon",

                                                "maxPoss Ckov Kaon; x [cm]; y [cm]",

                                                1600, 0, 300, 1440, 0, 143));

        std::unique_ptr<TH2F> hMinProton(

            new TH2F("minPoss Ckov Proton", "minPoss Ckov Proton; x [cm]; y [cm]",

                     1600, 0, 300, 1440, 0, 143));

        std::unique_ptr<TH2F> hMinPion(new TH2F("minPoss Ckov Pion",

                                                "minPoss Ckov Pion; x [cm]; y [cm]",

                                                1600, 0, 300, 1440, 0, 143));

        std::unique_ptr<TH2F> hMinKaon(new TH2F("minPoss Ckov Kaon",

                                                "minPoss Ckov Kaon; x [cm]; y [cm]",

                                                1600, 0, 300, 1440, 0, 143));

        // ... similarly for other TH2F objects you might have

        std::unique_ptr<TH2F> herrPos(new TH2F("errPos", "errPos; x [cm]; y [cm]",

                                               1600, 0, 300, 1440, 0, 143));


        // array");

        fillMapFromVec(hMaxPion.get(), arrMaxPionPos); // map, array

        fillMapFromVec(hMaxKaon.get(), arrMaxKaonPos); // map, array

        fillMapFromVec(hMaxProton.get(), arrMaxProtonPos); // map, array

        fillMapFromVec(hMinPion.get(), arrMinPionPos); // map, array

        fillMapFromVec(hMinKaon.get(), arrMinKaonPos); // map, array

        fillMapFromVec(hMinProton.get(), arrMinProtonPos);

        //// limMin->SetMarkerStyle(2);

        //// limMax->SetMarkerStyle(2);

        // limMin->SetMarkerColor(kOrange-3);

        // limMax->SetMarkerColor(kOrange+2);

        // limMin->SetMarkerColor(kCyan);

        hCkovCandMapRange->SetMarkerColor(kGreen);

        hCkovCandMapOutRange->SetMarkerColor(kGreen + 2);

        // herrPos->SetMarkerColor(kRed+3);

        // herrPos->SetMarkerStyle(3);

        // herrPos->Fill(errPos->X(), errPos->Y());

        hmipSizeFilter->SetMarkerColor(kOrange - 1);

        hmipChargeFilter->SetMarkerColor(kOrange + 1);

        hCkovCandMapRange->SetMarkerStyle(3);

        hCkovCandMapOutRange->SetMarkerStyle(2);

        hmipSizeFilter->SetMarkerStyle(2);

        hmipChargeFilter->SetMarkerStyle(2);

        int xrange = 100;

        int yrange = 100;

        int xMin = xMip - xrange;

        int xMax = xMip + xrange;

        int yMin = yMip - yrange;

        int yMax = yMip + yrange;

        // xMin = 0;

        // xMax = 160*0.8;

        // yMin = 0;

        // yMax = 144*0.84;

        auto st2 = Form("All clusters pdg %d ; x [cm]; y [cm]", trackPdg);

        std::unique_ptr<TH2F> totalCluMap(

            new TH2F("All clusters", st2, 320, 0, 159, 288, 0, 143));

        totalCluMap->SetAxisRange(xMin, xMax, "X");

        totalCluMap->SetAxisRange(yMin, yMax, "Y");

        totalCluMap->SetMarkerStyle(2);

        for (const auto &c : clusterTrack)
        {

            totalCluMap->Fill(c.mX, c.mY, c.mQ);
        }

        hCkovCandMapRange->SetAxisRange(xMin, xMax, "X");

        hCkovCandMapRange->SetAxisRange(yMin, yMax, "Y");

        hmipSizeFilter->SetAxisRange(xMin, xMax, "X");

        hmipSizeFilter->SetAxisRange(yMin, yMax, "Y");

        hCkovCandMapOutRange->SetAxisRange(xMin, xMax, "X");

        hCkovCandMapOutRange->SetAxisRange(yMin, yMax, "Y");

        hmipChargeFilter->SetAxisRange(xMin, xMax, "X");

        hmipChargeFilter->SetAxisRange(yMin, yMax, "Y");

        hSignalAndNoiseMap->SetAxisRange(xMin, xMax, "X");

        hSignalAndNoiseMap->SetAxisRange(yMin, yMax, "Y");

        hSignalMIP->SetAxisRange(xMin, xMax, "X");

        hSignalMIP->SetAxisRange(yMin, yMax, "Y");

        hSignalMIPpc->SetAxisRange(xMin, xMax, "X");

        hSignalMIPpc->SetAxisRange(yMin, yMax, "Y");

        hMaxProton->SetAxisRange(xMin, xMax, "X");

        hMaxProton->SetAxisRange(yMin, yMax, "Y");

        hMaxPion->SetAxisRange(xMin, xMax, "X");

        hMaxPion->SetAxisRange(yMin, yMax, "Y");

        hMaxPionMinL->SetAxisRange(xMin, xMax, "X");

        hMaxPionMinL->SetAxisRange(yMin, yMax, "Y");

        hMinPionMaxL->SetAxisRange(xMin, xMax, "X");

        hMinPionMaxL->SetAxisRange(yMin, yMax, "Y");

        hMaxKaon->SetAxisRange(xMin, xMax, "X");

        hMaxKaon->SetAxisRange(yMin, yMax, "Y");

        hMinProton->SetAxisRange(xMin, xMax, "X");

        hMinProton->SetAxisRange(yMin, yMax, "Y");

        hMinPion->SetAxisRange(xMin, xMax, "X");

        hMinPion->SetAxisRange(yMin, yMax, "Y");

        hMinKaon->SetAxisRange(xMin, xMax, "X");

        hMinKaon->SetAxisRange(yMin, yMax, "Y");

        fillMapFromVec(hCkovCandMapRange.get(), ckovCandMapRange); // map, array

        fillMapFromVec(hCkovCandMapOutRange.get(),

                       ckovCandMapOutRange); // map, array

        fillMapFromVec(hmipSizeFilter.get(), mipSizeFilter); // map, array

        fillMapFromVec(hmipChargeFilter.get(), mipChargeFilter);

        hMaxProton->SetMarkerColor(kGreen + 4);

        hMaxKaon->SetMarkerColor(kRed + 1);

        hMinProton->SetMarkerColor(kGreen + 3);

        hMinKaon->SetMarkerColor(kRed);

        hMaxPion->SetMarkerColor(kBlue + 4); // max ckov max L

        hMinPion->SetMarkerColor(kBlue + 3); // min ckov min L/**/

        hMinPionMaxL->SetMarkerColor(kBlue); // min ckov max L

        hMaxPionMinL->SetMarkerColor(kBlue + 4); // max ckov min L

        /*

        hMaxPion->SetMarkerStyle(3);

        hMinPionMaxL->SetMarkerStyle(3);*/

        hSignalMIPpc->SetMarkerStyle(3);

        hSignalMIP->SetMarkerStyle(3);

        hSignalMIPpc->SetMarkerStyle(3);

        hSignalMIPpc->SetMarkerColor(kRed);

        hSignalAndNoiseMap->SetMarkerStyle(2);

        auto trkPCMap =

            std::make_unique<TH2F>("trkPCMap ", "trkPCMap; x [cm]; y [cm]",

                                   160 * 10, 0., 300., 144 * 10, 0, 143);

        auto trkRadMap =

            std::make_unique<TH2F>("trkRadMap ", "trkRadMap; x [cm]; y [cm]",

                                   160 * 10, 0., 300., 144 * 10, 0, 143);

        auto MIP = std::make_unique<TH2F>("MIP ", "MIP; x [cm]; y [cm]", 160 * 10,

                                          0., 300., 144 * 10, 0, 143);

        trkRadMap->Fill(trkRad.X(), trkRad.Y());

        trkPCMap->Fill(trkPC.X(), trkPC.Y());

        MIP->Fill(xMip, yMip);

        MIP->SetMarkerColor(kCyan); // max ckov max L

        MIP->SetMarkerStyle(3); // max ckov max L

        trkRadMap->SetMarkerStyle(3);

        trkPCMap->SetMarkerStyle(2);

        tLineMIP->SetLineColor(kCyan);

        tLineTRK->SetLineColor(kRed);

        auto tcnvRane =

            std::make_unique<TCanvas>(Form("tcnvRane%d", plotNumber),

                                      Form("tcnvRane%d", plotNumber), 1600, 800);

        tcnvRane->Divide(2);

        tcnvRane->cd(1);

        hMaxProton->SetMarkerColor(kGreen + 4);

        hMaxKaon->SetMarkerColor(kRed + 1);

        hMinProton->SetMarkerColor(kGreen + 3);

        hMinKaon->SetMarkerColor(kRed);

        hMaxPion->SetMarkerColor(kBlue + 4); // max ckov max L

        hMinPion->SetMarkerColor(kBlue + 3); // min ckov min L/**/

        hMinPionMaxL->SetMarkerColor(kBlue); // min ckov max L

        hMaxPionMinL->SetMarkerColor(kBlue + 4); // max ckov min L

        TLegend *legend = new TLegend(0.758, 0.758, 0.975, 0.925); //

        // Add entries to the legend

        legend->AddEntry(trkRadMap.get(), "rad", "p"); //

        legend->AddEntry(trkPCMap.get(), "MIP", "p");

        legend->AddEntry(hCkovCandMapOutRange.get(), "Out of region", "p");

        legend->AddEntry(hmipSizeFilter.get(), "Size > 2 && Charge > 200", "p");

        legend->AddEntry(hCkovCandMapRange.get(), "Cluster Candidates", "p");

        legend->AddEntry(hMaxProton.get(), "hMaxProton", "p"); //

        legend->AddEntry(hMaxKaon.get(), "hMaxKaon", "p");

        legend->AddEntry(hMinProton.get(), "OhMinProton", "p");

        hCkovCandMapRange->Draw();

        hCkovCandMapRange->SetStats(kFALSE);

        hCkovCandMapOutRange->Draw("same");

        hmipSizeFilter->Draw("same");

        hmipChargeFilter->Draw("same");

        // MIP->Draw("same");

        legend->Draw("same");

        hMaxPion->Draw("same"); //  Printf("makeEvent()  hMaxPion->Draw");

        hMinPion->Draw("same");

        hMaxProton->Draw("same"); // Printf("makeEvent()  hMaxProton->Draw");

        hMinProton->Draw("same");

        hMinKaon->Draw("same");

        hMaxKaon->Draw("same");

        tLineTRK->Draw("same");

        // tLineMIP->Draw("same");

        trkRadMap->SetMarkerStyle(20); // E.g., a filled circle

        trkRadMap->SetMarkerColor(kRed); // E.g., a red color

        trkRadMap->SetContour(99); // Use a high number for detailed contours

        trkRadMap->SetLineWidth(4); ///

        trkRadMap->SetLineColor(

            kRed); // Setting a distinct color for the contour lines

        trkRadMap->Draw("CONT3 same");

        trkPCMap->Draw("same");

        // herrPos->Draw("same");

        // printf("opulatePtr->// limMin->GetEntries() num BinEntries = %f ",

        // populatePtr->// limMin->GetEntries());

        // printf("end// limMax num BinEntries = %f ", // limMax->GetEntries());

        // printf("emd// limMin num BinEntries = %f ", // limMin->GetEntries());

        auto textNumPion = new TLatex(

            xMin, yMin + 5, ("Pions: " + std::to_string(numPion)).c_str());

        textNumPion->SetTextSize(0.04);

        //  textNumPion->Draw("same");

        auto textNumKaon = new TLatex(xMin + 10, yMin + 5,

                                      (" K: " + std::to_string(numKaon)).c_str());

        textNumKaon->SetTextSize(0.04);

        //  textNumKaon->Draw("same");

        auto textNumProton = new TLatex(

            xMin + 20, yMin + 5, (" Pr: " + std::to_string(numProton)).c_str());

        textNumProton->SetTextSize(0.04);

        //  textNumProton->Draw("same");

        auto textNumTotal = new TLatex(

            xMin + 30, yMin + 5, ("To : " + std::to_string(numTotal)).c_str());

        textNumTotal->SetTextSize(0.04);

        // textNumTotal->Draw("same");

        TPad *pad2 = static_cast<TPad *>(tcnvRane->cd(2));

        pad2->SetRightMargin(.055 + pad2->GetRightMargin());

        pad2->SetLeftMargin(-.055 + pad2->GetLeftMargin());

        pad2->SetLogz(1);

        totalCluMap->SetStats(kFALSE);

        totalCluMap->Draw("Colz");

        //	pad2->Clear();

        trkRadMap->SetLineWidth(4); ///

        trkRadMap->SetMarkerStyle(20); // E.g., a filled circle

        trkRadMap->SetMarkerColor(kRed); // E.g., a red color

        trkRadMap->SetContour(99); // Use a high number for detailed contours

        trkRadMap->SetLineColor(

            kRed); // Setting a distinct color for the contour lines

        trkRadMap->Draw("CONT3 same"); // Drawing the contours on the same pad

        tcnvRane->SaveAs(Form("Segmented_%dCkov%.2f.png", plotNumber, sigSep));

        // limMin->Draw("same");

        // limMax->Draw("same");

        // tcnvRane->SaveAs(Form("Segmented%d.png", plotNumber));

        plotNumber++;
    }

    /*

          void drawMaxRegions()

          {





                  const auto trkPC = populatePtr->getPcImp();

                  const auto trkRad = populatePtr->getTrackPos();

                  TH2F* trkPCMap = new TH2F("trkPCMap ", "trkPCMap; x [cm]; y

       [cm]",160*20,0.,300.,144*20,0,143); TH2F* trkRadMap = new TH2F("trkRadMap

       ", "trkRadMap; x [cm]; y [cm]",160*20,0.,300.,144*20,0,143);

            trkRadMap->Fill(trkRad.X(), trkRad.Y());

            trkPCMap->Fill(trkPC.X(), trkPC.Y());





                  trkRadMap->SetMarkerStyle(3);

       trkRadMap->SetMarkerColor(kGreen+4); trkPCMap->SetMarkerStyle(3);

       trkPCMap->SetMarkerColor(kGreen+2);



                  TCanvas *thSignalAndNoiseMap = new

       TCanvas(Form("hSignalAndNoiseMap%d", eventCnt),Form("hSignalAndNoiseMap%d",

       eventCnt),800,800); thSignalAndNoiseMap->cd(); hSignalAndNoiseMap->Draw();





                  hMaxPion->Draw("same");     //  Printf("makeEvent()

       hMaxPion->Draw"); hMinPion->Draw("same"); hMaxProton->Draw("same");

       //Printf("makeEvent()  hMaxProton->Draw"); hMinProton->Draw("same");

                  hMinKaon->Draw("same");

                  hMaxKaon->Draw("same");



                  hSignalMIP->Draw("same");

                  hSignalMIPpc->Draw("same");



                  trkRadMap->Draw("same");

                  trkPCMap->Draw("same");

          }







          void drawTotalMapAndMaxRegions()

          {

                  TCanvas* tcnvRane = new TCanvas(Form("TotalMap%d", eventCnt),

       Form("TotalMap%d", eventCnt), 1600, 800);



                                  const auto trkPC = populatePtr->getPcImp();

                  const auto trkRad = populatePtr->getTrackPos();

                  TH2F* trkPCMap = new TH2F("trkPCMap ", "trkPCMap; x [cm]; y

       [cm]",160*20,0.,159.,144*20,0,143); TH2F* trkRadMap = new TH2F("trkRadMap

       ", "trkRadMap; x [cm]; y [cm]",160*20,0.,159.,144*20,0,143);

            trkRadMap->Fill(trkRad.X(), trkRad.Y());

            trkPCMap->Fill(trkPC.X(), trkPC.Y());





                  trkRadMap->SetMarkerStyle(3);

       trkRadMap->SetMarkerColor(kGreen+4); trkPCMap->SetMarkerStyle(3);

       trkPCMap->SetMarkerColor(kGreen+2);



                  tcnvRane->cd();

                  ckovCandMapRange->Draw();

                  ckovCandMapOutRange->Draw("same");

                  mipSizeFilter->Draw("same");

                  mipChargeFilter->Draw("same");



                  hSignalAndNoiseMap->Draw();





                  hMaxPion->Draw("same");     //  Printf("makeEvent()

       hMaxPion->Draw"); hMinPion->Draw("same"); hMaxProton->Draw("same");

       //Printf("makeEvent()  hMaxProton->Draw"); hMinProton->Draw("same");

                  hMinKaon->Draw("same");

                  hMaxKaon->Draw("same");



                  hSignalMIP->Draw("same");

                  hSignalMIPpc->Draw("same");

                                  trkRadMap->Draw("same");

                  trkPCMap->Draw("same");

          } */
};

class CkovTools
{

private:
    // for storing candidates...

    struct Candidate
    {

        double x, y = 0.;

        double R = 0.;

        double phiL = 0., phi = 0.;

        bool isCandidate = false;
    };

    bool print = false; // ef: TODO : later pass this in ctor

    std::unique_ptr<Populate2> populatePtr, populatePtrCp, populate2Ptr,

        populatePtrOuter, populatePtrInner;

    // using array = std::array;

    using vecArray4 = std::vector<std::array<double, 4>>;

    using vecArray3 = std::vector<std::array<double, 3>>;

    using vecArray2 = std::vector<std::array<double, 2>>;

    using vecArrayPair = std::vector<std::pair<double, double>>;

    // using arrArray3 = std::vector<std::array<double,3>>;

    using segType = std::vector<std::pair<double, double>>;

    segType segPionLocal = {{}};

    bool kaonStatus = true, pionStatus = true, protonStatus = true;

    // TLine* tlinePion;

    float mSigmaSep = 1.5;

    // used in SigCrom

    //  double f = 0.00928*(7.75-5.635)/TMath::Sqrt(12.);

    // static constexpr double f = 0.0172*(7.75-5.635)/TMath::Sqrt(24.);

    static constexpr double sq6 = 2.44948974278;

    static constexpr double f = 0.0172 * (7.75 - 5.635) / (2 * sq6);

    static constexpr double PI = M_PI;

    static constexpr double halfPI = M_PI / 2;

    static constexpr double twoPI = M_PI * 2;

    static constexpr double stdDevPion = 0.001;

    static constexpr double stdDevKaon = 0.001;

    static constexpr double stdDevProton = 0.001;

    static constexpr float tGap = 8;

    static constexpr float rW = 1.5; // was 1?,

    static constexpr float qW = 0.5;

    static constexpr float lMax = 1.5;

    static constexpr float CH4GapWidth = 8;

    static constexpr float RadiatorWidth = 1.5; // was 1?

    static constexpr float QuartzWindowWidth = 0.5;

    static constexpr float L_CONST = rW / 2;

    static constexpr double dnDE =

        0.0172; // σE = (dn/dE)σdet dn/dE is 0.0172 eV−1.

    // σdet represents the standard deviation of the detected

    // Cherenkov photon spectrum

    // for simulation :  mean photonergny assumed to be cosntant

    static constexpr double photEnergyMean =

        6.82 / 1000.; // jsut based on some runs

    // 1: chromacity error

    // std-dev of single-photon angular resolution, contribution from fluctuation

    // in photon-energy

    static constexpr double sigmaE = dnDE * photEnergyMean; // given in Rad

    // 4 :

    // The track incidence angle error, related to the particle angle θp and to

    // the precision of the tracking devices. In the following, the θp error,

    // assumed to be of the order of 2 mrad at the considered incidence angles,

    // will not be quoted in tables and plots but simply included in the

    // calculation of the total angular resolution.

    static constexpr double sigmaThetaP = 0.002; // given in Rad

    // combined angular resolution

    static constexpr float sigmaAnglesSq =

        sigmaE * sigmaE + sigmaThetaP * sigmaThetaP;

    const float sigmaAngles =

        TMath::Sqrt(sigmaE * sigmaE + sigmaThetaP * sigmaThetaP);

    // (3) The localization error, related to the precision with which the photon

    // and particle impact coordinates can be measured. It is determined by the

    // detector characteristics (pad size, sense wire pitch) and by the photon

    // feedback.

    // error of position of clusters

    static constexpr float sigmaR = .2; // TDR : sizeY / sqrt(12)

    static constexpr float sigmaRsq = .2 * .2; // TDR : sizeY / sqrt(12)

    // (2) The geometric error, related to the spread of the emission point along

    // the particle path in the Cherenkov radiator. It depends on the ratio RW/GAP

    // between the radiator thickness, RW, and he proximity gap width, GAP; it can

    // be minimized by increasing GAP and reducing RW, provided the number of

    // photoelectrons per ring is sufficient for pattern recognition (Chapter 4).

    const float sigmaLfactor =

        rW / TMath::Sqrt(12); // sigmaL = sigmaLfactor/ cosThetaP

    float sigmaL = 0.0;

    float sigmaLsq = 0.0;

    // ef : set this constexpr ins

    float L = rW / 2;

    // L value for reconstruction

    static constexpr float EmissionLenght = RadiatorWidth / 2;

    float thetaP, phiP, xPC, yPC, xRad, yRad;

    float nF, nQ, nG;

    std::array<float, 3> ckovHypsMin, ckovHypsMax;

    std::vector<std::pair<double, double>> photons;

    double ckovProton = 999, ckovKaon = 999, ckovPion = 999;

    // instantiate ckovhyp boundaries to "invalid values"

    // --> later set them if they exceed momentum threshold

    double ckovPionMin = 999, ckovPionMax = 0, ckovKaonMin = 999, ckovKaonMax = 0,

           ckovProtonMin = 999, ckovProtonMax = 0;

    // double mRMax, mL1Max, mL2Max, mRMin, mL1Min, mL2Min;

    double momentum, mass;

    float cosThetaP, sinThetaP, tanThetaP;

    float cosPhiP, sinPhiP, tanPhiP;

    float trackCkov;

    double xMipLocal, yMipLocal;

    float phiLCurrent, etaCCurrent;

    float xMip, yMip, qMip;

    TVector2 trkPos; // trk at RAD

    TVector3 trkDir; // trk mag theta phi

    TVector2 mipPos; // MIP PC

    int trackPdg;

    string trackPdgString;

    int eventCnt;

    TVector2 trkPC;

public:
    typedef std::vector<std::pair<double, double>> MapType;



    // set new values of phiL and etaC for new photon



    //   double pc[5] = {xRad,yRad,L,thetaP, phiP};

    double refIndexes[3] = {nF, nQ, nG};

    ~CkovTools() {}

    // ef: let trackCkov be empty aon; TODO: add this based on pdg!

    CkovTools(double radParams[7], double refIndexes[3], double MIP[3],

              std::array<float, 3> ckovHypsMin, std::array<float, 3> ckovHypsMax,

              float trackCkov, int eventCnt, int _trackPdg, float sigmaSep)

        : ckovHypsMin(ckovHypsMin), ckovHypsMax(ckovHypsMax),

          trackCkov(trackCkov), eventCnt(eventCnt)
    {

        mSigmaSep = sigmaSep;

        trackPdg = _trackPdg;

        trackPdgString = getPDG(trackPdg);

        // double radParams[6] = {xRad,yRad,L,thetaP, phiP, randomValue.momentum};

        xMip = MIP[0], yMip = MIP[1], qMip = MIP[2];

        xRad = radParams[0];

        yRad = radParams[1];

        L = radParams[2];

        thetaP = radParams[3];

        phiP = radParams[4];

        momentum = radParams[5];

        mass = radParams[6]; // ef; let this be empty aon! TODO: get this from pdg

        nF = refIndexes[0];

        trkPos.Set(xRad, yRad); // track positon in LORS at RAD   // XY mag

        const double winThick = 0.5, radThick = 1.5;

        const int gapThick = 8;

        const double gapIdx = 1.0005, winIdx = 1.583; // inIdx = 1.5787;

        double zRad =

            -0.5 * radThick - 0.5 * winThick; // z position of middle of RAD

        TVector3 rad(trkPos.X(), trkPos.Y(), zRad); // impact point at middle of RAD

        TVector3 pc(xMip, yMip, 0.5 * winThick + gapThick); // mip at PC

        // Printf("Phi  %.2f Thta %.2f of Track", phiP, thetaP);

        phiP = (pc - rad).Phi();

        thetaP = (pc - rad).Theta();

        Printf("Phi  %.2f Thta %.2f of rad--MIP", phiP, thetaP);

        trkPC.Set(xMip, yMip); // MIP pos at PC

        mipPos.Set(xMip, yMip); // MIP pos at PC

        trkDir;

        trkDir.SetMagThetaPhi(1, thetaP, phiP); // track direction in LORS at RAD

        nF = 1.2928 - 0.0025; // ef got this from 1 run .. assuming T = 20 for sim

        auto nFstd2 = 0.01; // 2 times std-dev for tejh run


        nQ = 1.583; // inIdx = 1.5787;

        nG = 1.005;

        trackPdgString = getPDG(trackPdg);

        if (TMath::IsNaN(ckovHypsMax[0]))
        {

            // printf("Pion CkovHyps is Nan!");

            setPionStatus(false);

            // setIsNan()?
        }
        else
        {

            setPionStatus(true);

            // printf("Pion CkovHyps %.2f", ckovHyps[0]);
        }

        if (TMath::IsNaN(ckovHypsMax[1]))
        {

            // printf("Kaon CkovHyps is Nan!");

            setKaonStatus(false);

            // setIsNan()?
        }
        else
        {

            // printf("Kaon CkovHyps %.2f", ckovHyps[1]);
        }

        if (TMath::IsNaN(ckovHypsMax[2]))
        {

            // printf("Proton CkovHyps is Nan!");

            setProtonStatus(false);
        }



        cosThetaP = TMath::Cos(thetaP);

        sinThetaP = TMath::Sin(thetaP);

        tanThetaP = TMath::Tan(thetaP);

        cosPhiP = TMath::Cos(phiP);

        sinPhiP = TMath::Sin(phiP);

        // set the geometric error std-dev

        sigmaL = sigmaLfactor / cosThetaP;

        sigmaLsq = sigmaL * sigmaL;

        TRotation rotZ;

        rotZ.RotateZ(phiP);

        TRotation rotY;

        rotY.RotateY(thetaP);

        TVector3 ip(0, 0, rW - L + tGap + qW);

        TVector3 op;

        op = rotZ * rotY * ip;

        xMipLocal = tanThetaP * cosPhiP * (rW - L + tGap + qW);

        yMipLocal = tanThetaP * sinPhiP * (rW - L + tGap + qW);

        auto dX = tanThetaP * cosPhiP * (rW - L + tGap + qW);

        auto dY = tanThetaP * sinPhiP * (rW - L + tGap + qW);

        xPC = dX + xRad; // bruke PC eller MIP?

        yPC = dY + yRad;
    }

    void setPionStatus(bool status) { pionStatus = status; }

    bool getPionStatus() const { return pionStatus; }

    void setKaonStatus(bool status) { kaonStatus = status; }

    bool getKaonStatus() const { return kaonStatus; }

    void setProtonStatus(bool status) { protonStatus = status; }

    bool getProtonStatus() const { return protonStatus; }

    double getCkovPion() { return ckovPion; }

    double getCkovKaon() { return ckovKaon; }

    double getCkovProton() { return ckovProton; }

    double getMinCkovKaon() { return ckovKaonMin; }

    double getMaxCkovKaon() { return ckovKaonMax; }

    double getMinCkovPion() { return ckovPionMin; }

    double getMaxCkovPion() { return ckovPionMax; }

    double getMinCkovProton() { return ckovProtonMin; }

    double getMaxCkovProton() { return ckovProtonMax; }

    double getPhiP() { return phiP; }

    double getThetaP() { return thetaP; }


    std::vector<std::pair<double, double>>

    segment(std::vector<o2::hmpid::ClusterCandidate> &clusterTrack,

            std::array<int, 4> &arrayInfo, int trackIndex,float mipX, float mipY,

            float mipCharge, const int mcTrackPdg,

            const o2::dataformats::MatchInfoHMP &track, int trackNumber,

            int &plotNumber)
    {


        // Loop over all clusters in chamber     
        for (auto &photons : clusterTrack)
        {


            auto pdgString = getPDG(photons.mPDG);

            trackPdgString = getPDG(trackPdg);


            double thetaCer, phiCer;

            if (findPhotCkov(

                    photons.mX, photons.mY, thetaCer,

                    phiCer))
            { // find ckov angle for this  photon candidate

                // increment counter of photon candidates

                auto sigmaRing = aliSigma2(thetaP, phiP, thetaCer,

                                            phiCer); // rms of all contributing errors


                photons.mThetaCer = thetaCer;

                photons.mPhiCer = phiCer;

                photons.mSigmaRing = sigmaRing;


            }


        } 

    } // end segment



}; // end CkovTools class 

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/*
    string getPDG(int pdg)
    {
        std::string pdgString;

        switch (TMath::Abs(pdg))
        {
        case 11:
            pdgString = "Electron";
            break;
        case 211:
            pdgString = "Pion";
            break;
        case 321:
            pdgString = "Kaon";
            break;
        case 2212:
            pdgString = "Proton";
            break;
        case 50000050:
            pdgString = "Photon";
        case 50000051:
            pdgString = "Photon";
        case 22:
            pdgString = "Photon";
            break;
        }

        return pdgString;
    }
    
    Double_t RadThick() const { return 1.5; } // Radiator thickness

    Double_t GapThick() const { return 8.0; } // Proximity gap thickness

    Double_t GetRefIdx() const { return 1.2905; } // running refractive index

    bool findPhotCkov(double cluX, double cluY, double &thetaCer, double &phiCer)
    {
        TVector3 dirCkov;
        double zRad = -0.5 * RadThick() - 0.5 * WinThick();
        TVector3 rad(trkPos.X(), trkPos.Y(), zRad);
        TVector3 pc(cluX, cluY, 0.5 * WinThick() + GapThick());
        double cluR = TMath::Sqrt((cluX - mipPos.X()) * (cluX - mipPos.X()) + (cluY - mipPos.Y()) * (cluY - mipPos.Y()));
        double phi = (pc - rad).Phi();
        double ckov1 = 0;
        double ckov2 = 0.75 + thetaP;
        const double kTol = 0.01;
        Int_t iIterCnt = 0;

        while (1)
        {
            if (iIterCnt >= 50)
            {
                return kFALSE;
            }

            double ckov = 0.5 * (ckov1 + ckov2);
            dirCkov.SetMagThetaPhi(1, ckov, phi);
            TVector2 posC = populatePtrInner->traceForward(dirCkov, 0.75);
            double dist = cluR - (posC - mipPos).Mod();

            if (posC.X() == -999)
            {
                dist = -999;
            }

            iIterCnt++;

            if (dist > kTol)
            {
                ckov1 = ckov;
            }
            else if (dist < -kTol)
            {
                ckov2 = ckov;
            }
            else
            {
                dirCkov.SetMagThetaPhi(1, ckov, phi);
                populatePtrInner->lors2Trs(dirCkov, thetaCer, phiCer);
                return kTRUE;
            }
        }
    }
    */