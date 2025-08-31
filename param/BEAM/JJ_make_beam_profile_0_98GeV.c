// make_beam_profile_0.98GeV.C
void JJ_make_beam_profile_0_98GeV() {
    // 1. 새로운 ROOT 파일을 생성합니다. (파일 이름에 운동량 명시)
    TFile *new_file = new TFile("beam_profile_0_98GeV.root", "RECREATE");
    TTree *new_tree = new TTree("tr", "Ideal Gaussian Beam Profile for E45 at 0.98GeV/c");

    // 2. TTree에 저장할 C++ 변수들을 선언합니다.
    double pointInx, pointIny, pointInz;
    double pInx, pIny, pInz;

    // 3. Branch 생성
    new_tree->Branch("pointInx", &pointInx, "pointInx/D");
    new_tree->Branch("pointIny", &pointIny, "pointIny/D");
    new_tree->Branch("pointInz", &pointInz, "pointInz/D");
    new_tree->Branch("pInx", &pInx, "pInx/D");
    new_tree->Branch("pIny", &pIny, "pIny/D");
    new_tree->Branch("pInz", &pInz, "pInz/D");

    // 4. 빔 파라미터 (Mean은 0으로, Sigma는 0.735GeV/c 데이터 유지)
    int n_events = 500000;
    double x_mean = 0.0;        // 타겟 중심을 향하도록 0으로 설정
    double x_sigma = 23.6732;
    double y_mean = 0.0;        // 타겟 중심을 향하도록 0으로 설정
    double y_sigma = 16.9384;
    double u_mean = 0.0;        // 빔이 z축에 평행하도록 0으로 설정
    double u_sigma = 0.00585364;
    double v_mean = 0.0;        // 빔이 z축에 평행하도록 0으로 설정
    double v_sigma = 0.00365566;

    // [가장 중요한 수정]
    // 빔의 전체 운동량 크기를 0.98 GeV/c로 변경합니다.
    double p_total = 0.98;

    // 5. 이벤트 루프
    for (int i = 0; i < n_events; ++i) {
        pointInx = gRandom->Gaus(x_mean, x_sigma);
        pointIny = gRandom->Gaus(y_mean, y_sigma);
        pointInz = -1100.0; // Z 시작점은 그대로 유지

        double u = gRandom->Gaus(u_mean, u_sigma);
        double v = gRandom->Gaus(v_mean, v_sigma);

        // 새로운 p_total 값으로 운동량 성분 재계산
        pInz = p_total / sqrt(u*u + v*v + 1.0);
        pInx = u * pInz;
        pIny = v * pInz;

        new_tree->Fill();
    }

    // 6. 파일 저장
    new_tree->Write();
    new_file->Close();

    cout << "New beam profile 'beam_profile_0.98GeV.root' created." << endl;
}
