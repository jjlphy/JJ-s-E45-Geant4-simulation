// make_beam_profile.C (최종 수정본)
void JJ_make_beam_profile() {
    // 1. 새로운 ROOT 파일을 생성합니다.
    TFile *new_file = new TFile("beam_profile_e45_gaus.root", "RECREATE");
    TTree *new_tree = new TTree("tr", "Ideal Gaussian Beam Profile for E45"); // Tree 이름을 'tr'로 변경!

    // 2. TTree에 저장할 C++ 변수들을 선언합니다. Branch 이름과 맞춰줍니다.
    double pointInx, pointIny, pointInz;
    double pInx, pIny, pInz;

    // 3. Branch 생성: Geant4가 기대하는 이름으로 Branch를 만듭니다.
    new_tree->Branch("pointInx", &pointInx, "pointInx/D");
    new_tree->Branch("pointIny", &pointIny, "pointIny/D");
    new_tree->Branch("pointInz", &pointInz, "pointInz/D");
    new_tree->Branch("pInx", &pInx, "pInx/D");
    new_tree->Branch("pIny", &pIny, "pIny/D");
    new_tree->Branch("pInz", &pInz, "pInz/D");

    // 4. 피팅 결과로 얻은 빔 파라미터
    int n_events = 500000;
    double x_mean = 12.6788;
    double x_sigma = 23.6732;
    double y_mean = -0.290646;
    double y_sigma = 16.9384;
    double u_mean = 0.02115;
    double u_sigma = 0.00585364;
    double v_mean = 0.000409;
    double v_sigma = 0.00365566;

    // 빔의 전체 운동량 크기 (단위: GeV/c)
    double p_total = 0.735;

    // 5. 이벤트 루프: 랜덤 값을 생성하고 계산하여 Tree를 채웁니다.
    for (int i = 0; i < n_events; ++i) {
        // 위치(x, y)와 방향(u, v)을 가우시안 분포로 랜덤하게 생성
        pointInx = gRandom->Gaus(x_mean, x_sigma);
        pointIny = gRandom->Gaus(y_mean, y_sigma);
        pointInz = -1100.0; // 👈 0.0이 아니라, 빔의 시작점 Z좌표(단위: mm)로 수정!

        double u = gRandom->Gaus(u_mean, u_sigma);
        double v = gRandom->Gaus(v_mean, v_sigma);

        // u, v와 전체 운동량(p_total)으로부터 운동량 성분(pInx, pIny, pInz) 계산
        // p_total^2 = pInx^2 + pIny^2 + pInz^2
        // u = pInx/pInz, v = pIny/pInz
        pInz = p_total / sqrt(u*u + v*v + 1.0);
        pInx = u * pInz;
        pIny = v * pInz;

        new_tree->Fill(); // 계산된 값들로 Tree의 한 줄을 채웁니다.
    }

    // 6. Tree를 파일에 쓰고 파일을 닫습니다.
    new_tree->Write();
    new_file->Close();

    cout << "New beam profile 'beam_profile_e45_gaus.root' created with correct branches." << endl;
}