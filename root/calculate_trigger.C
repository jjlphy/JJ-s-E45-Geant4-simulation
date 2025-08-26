#include <TFile.h>
#include <TTree.h>
#include <TParticle.h>
#include <vector>
#include <iostream>
#include <string>
#include <iomanip> // 출력 포맷팅을 위해 추가

#pragma link C++ class vector<TParticle>+;

void calculate_trigger() {
    TFile *f = gFile;
    if (!f || f->IsZombie()) { /* ... */ return; }
    TTree *tree = (TTree*)f->Get("g4hyptpc");
    if (!tree) { /* ... */ return; }

    // 1. 필요한 모든 브랜치에 변수 연결
    std::vector<TParticle> *t0_hits = nullptr;
    std::vector<TParticle> *bh2_hits = nullptr;
    std::vector<TParticle> *htof_hits = nullptr;
    std::vector<TParticle> *bvhd_hits = nullptr;
    std::vector<TParticle> *bvhu_hits = nullptr;

    tree->SetBranchAddress("T0", &t0_hits);
    tree->SetBranchAddress("BH2", &bh2_hits);
    tree->SetBranchAddress("HTOF", &htof_hits);
    tree->SetBranchAddress("BVH_D", &bvhd_hits);
    tree->SetBranchAddress("BVH_U", &bvhu_hits);

    // 2. 각 단계별 카운터 변수 초기화
    long long total_events = tree->GetEntries();
    long long count_step1 = 0; // T0 x BH2 (분석 기준)
    long long count_step2 = 0; // ... x HTOF multi>=2
    long long count_step3 = 0; // ... x BVH_D
    long long count_step4 = 0; // ... x BVH_U

    std::cout << "총 " << total_events << "개의 이벤트를 분석합니다..." << std::endl;

    // 3. 모든 이벤트를 순회하며 단계별 조건 확인 (계산 로직은 이전과 동일)
    for (long long i = 0; i < total_events; ++i) {
        tree->GetEntry(i);

        if (t0_hits && t0_hits->size() > 0 && bh2_hits && bh2_hits->size() > 0) {
            count_step1++;
            if (htof_hits && htof_hits->size() >= 2) {
                count_step2++;
                if (bvhd_hits && bvhd_hits->size() > 0) {
                    count_step3++;
                    if (bvhu_hits && bvhu_hits->size() > 0) {
                        count_step4++;
                    }
                }
            }
        }
    }
    std::cout << "분석 완료!" << std::endl;

    // 4. ## 수정된 최종 결과 테이블 ##
    // 분모(Baseline)를 Beam Trigger를 통과한 이벤트 수로 설정
    long long baseline_events = count_step1;

    std::cout << "\n--- 최종 트리거 효율 분석 (기준: T0 && BH2) ---" << std::endl;
    std::cout << std::fixed << std::setprecision(3); // 소수점 3자리까지

    std::cout << "분석 기준(Baseline) 이벤트 수: " << baseline_events << " 개 (100.000%)" << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    
    double eff_step2 = (baseline_events > 0) ? (double)count_step2 / baseline_events * 100.0 : 0;
    std::cout << " + [HTOF Multiplicity >= 2] 효율: " << std::setw(7) << eff_step2 << "%  (" << count_step2 << " / " << baseline_events << ")" << std::endl;

    double eff_step3 = (baseline_events > 0) ? (double)count_step3 / baseline_events * 100.0 : 0;
    std::cout << " + [... && BVH_D] 효율          : " << std::setw(7) << eff_step3 << "%  (" << count_step3 << " / " << baseline_events << ")" << std::endl;

    double eff_step4 = (baseline_events > 0) ? (double)count_step4 / baseline_events * 100.0 : 0;
    std::cout << " + [... && BVH_U] 효율          : " << std::setw(7) << eff_step4 << "%  (" << count_step4 << " / " << baseline_events << ")" << std::endl;
    
    std::cout << "--------------------------------------------------------" << std::endl;
}