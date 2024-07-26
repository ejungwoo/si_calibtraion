#include "si_calibration.h"

void run_energy_calibration(int run=-1, int chooseGate = 1, bool drawExample = true)
{
    MakeRun(run);

    auto fileHist = new TFile(fHistFileName,"read");
    for (auto det=0; det<40; ++det)
        for (auto side=0; side<2; ++side)
            for (auto strip=0; strip<8; ++strip)
                fHistEnergy[det][side][strip] = (TH1D*) fileHist -> Get(MakeHistName("Energy",det,side,strip));

    for (auto det=0; det<40; ++det)
    {
        cout << det << endl;

        TCanvas* cvs = nullptr;

        int iCvs = 1;
        bool drawCurrentDet = false;
        if (drawExample && (det==fExampleDet3))
            drawCurrentDet = true;

        if (drawCurrentDet) {
            cvs = MakeCanvas("LeftRightGateFit",det);
        }

        for (auto strip=0; strip<8; ++strip)
        {
            for (auto side=0; side<2; ++side)
            {
                auto hist = fHistEnergy[det][side][strip];

                if (drawCurrentDet) {
                    cvs -> cd(iCvs++);
                    hist -> Draw();
                }
            }
        }
    }
}
