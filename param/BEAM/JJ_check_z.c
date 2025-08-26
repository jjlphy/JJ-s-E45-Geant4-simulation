// check_z.C
void JJ_check_z() {
    // 1. 원본 빔 프로파일 파일을 엽니다.
    TFile *f = new TFile("beam_profile_run332.root");
    if (!f || f->IsZombie()) {
        cout << "Error: Cannot open beam_profile_run332.root" << endl;
        return;
    }
    TTree *tree = (TTree*)f->Get("tr");
    if (!tree) {
        cout << "Error: Could not find TTree named 'tr'." << endl;
        f->Close();
        return;
    }

    // 2. pointInz 값을 담을 히스토그램을 생성합니다.
    //    범위는 넉넉하게 -5000mm ~ 0mm로 잡겠습니다.
    TH1D *h_z = new TH1D("h_z", "Beam Start Z Distribution; Z [mm]; Counts", 100, -5000, 0);

    // 3. TTree에서 pointInz 데이터만 뽑아서 히스토그램을 채웁니다.
    tree->Draw("pointInz >> h_z");

    // 4. 결과를 화면에 그립니다.
    TCanvas *c1 = new TCanvas("c1", "Z Position Check", 600, 400);
    h_z->Draw();

    // 5. 히스토그램의 평균값을 출력합니다 (모든 값이 같으므로 이 값이 바로 Z좌표입니다).
    double z_position = h_z->GetMean();
    cout << "\n=======================================" << endl;
    cout << "  Beam Start Z-Position is: " << z_position << " mm" << endl;
    cout << "=======================================" << endl;
}