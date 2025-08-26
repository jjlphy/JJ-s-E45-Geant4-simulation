// JJ_analyze_beam.C (최종 수정본)
void JJ_analyze_beam() {
    // 1. 파일 및 TTree 열기
    TFile *f = new TFile("beam_profile_run332.root");
    if (!f || f->IsZombie()) {
        cout << "오류: beam_profile_run332.root 파일을 열 수 없습니다." << endl;
        return;
    }
    TTree *tree = (TTree*)f->Get("tr");
    if (!tree) {
        cout << "오류: 파일 안에서 'tr'이라는 이름의 TTree를 찾을 수 없습니다." << endl;
        f->Close();
        return;
    }

    // 2. TTree의 Branch와 연결할 C++ 변수 선언
    double pointInx, pointIny, pInx, pIny, pInz;

    // 3. Branch 주소 설정: TTree의 각 Branch를 C++ 변수와 연결합니다.
    tree->SetBranchAddress("pointInx", &pointInx);
    tree->SetBranchAddress("pointIny", &pointIny);
    tree->SetBranchAddress("pInx", &pInx);
    tree->SetBranchAddress("pIny", &pIny);
    tree->SetBranchAddress("pInz", &pInz);

    // 4. 히스토그램 생성
    TH1D *h_x = new TH1D("h_x", "Beam X Distribution; X [mm]; Counts", 100, -50, 70);
    TH1D *h_y = new TH1D("h_y", "Beam Y Distribution; Y [mm]; Counts", 100, -50, 50);
    TH1D *h_u = new TH1D("h_u", "Beam U Distribution; U (p_x/p_z); Counts", 100, -0.01, 0.05);
    TH1D *h_v = new TH1D("h_v", "Beam V Distribution; V (p_y/p_z); Counts", 100, -0.01, 0.01);

    // 5. 이벤트 루프: TTree의 모든 이벤트를 하나씩 순회하며 데이터를 처리합니다.
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i); // i번째 이벤트 데이터를 불러와 C++ 변수에 채웁니다.

        // 히스토그램 채우기
        h_x->Fill(pointInx);
        h_y->Fill(pointIny);

        // pInz가 0이 아닐 때만 u, v를 계산 (0으로 나누는 오류 방지)
        if (pInz != 0) {
            h_u->Fill(pInx / pInz);
            h_v->Fill(pIny / pInz);
        }
    }

    // 6. 가우시안 피팅
    h_x->Fit("gaus");
    h_y->Fit("gaus");
    h_u->Fit("gaus");
    h_v->Fit("gaus");

    // 7. 결과 출력 (이전과 동일)
    cout << "\n=======================================" << endl;
    cout << "  Beam Profile Fitting Results" << endl;
    cout << "=======================================" << endl;
    cout << "Variable |    Mean (mu)    |  Sigma (sigma)" << endl;
    cout << "---------------------------------------" << endl;
    printf("   x     |  %12.4f   |  %12.4f \n", h_x->GetFunction("gaus")->GetParameter(1), h_x->GetFunction("gaus")->GetParameter(2));
    printf("   y     |  %12.4f   |  %12.4f \n", h_y->GetFunction("gaus")->GetParameter(1), h_y->GetFunction("gaus")->GetParameter(2));
    printf("   u     |  %12.4f   |  %12.4f \n", h_u->GetFunction("gaus")->GetParameter(1), h_u->GetFunction("gaus")->GetParameter(2));
    printf("   v     |  %12.4f   |  %12.4f \n", h_v->GetFunction("gaus")->GetParameter(1), h_v->GetFunction("gaus")->GetParameter(2));
    cout << "=======================================" << endl;

    // 8. 캔버스에 그리기 (이전과 동일)
    TCanvas *c1 = new TCanvas("c1", "Beam Profile Analysis", 800, 800);
    c1->Divide(2, 2);
    c1->cd(1); h_x->Draw();
    c1->cd(2); h_y->Draw();
    c1->cd(3); h_u->Draw();
    c1->cd(4); h_v->Draw();
}