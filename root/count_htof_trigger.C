#include <TFile.h>
#include <TTree.h>
#include <TParticle.h>
#include <vector>
#include <iostream>
#include <string>
#include <iomanip> // 출력 포맷팅을 위해 추가

#pragma link C++ class vector<TParticle>+;

void count_htof_trigger() {
    TFile *f = gFile;
    if (!f || f->IsZombie()) { /* ... */ return; }
    TTree *tree = (TTree*)f->Get("g4hyptpc");
    if (!tree) { /* ... */ return; }

    // 1. 필요한 브랜치에 변수 연결
    std::vector<TParticle> *t0_hits = nullptr;
    std::vector<TParticle> *bh2_hits = nullptr;
    std::vector<TParticle> *htof_hits = nullptr;

    tree->SetBranchAddress("T0", &t0_hits);
    tree->SetBranchAddress("BH2", &bh2_hits);
    tree->SetBranchAddress("HTOF", &htof_hits);

    // 2. 카운터 변수 초기화
    long long total_events = tree->GetEntries();
    long long beam_trigger_count = 0;   // 분모가 될 변수
    long long htof_trigger_count = 0;     // 분자가 될 변수

    std::cout << "총 " << total_events << "개의 이벤트를 분석합니다..." << std::endl;

    // 3. 모든 이벤트를 순회하며 조건 확인
    for (long long i = 0; i < total_events; ++i) {
        tree->GetEntry(i);

        // 단계 1: Beam Trigger 조건 확인 (T0 && BH2)
        if (t0_hits && t0_hits->size() > 0 && bh2_hits && bh2_hits->size() > 0) {
            beam_trigger_count++; // 분모 카운트

            // 단계 2: HTOF Multiplicity 조건 확인
            // Beam Trigger를 통과한 이벤트 중에서만 확인
            if (htof_hits && htof_hits->size() >= 2) {
                htof_trigger_count++; // 분자 카운트
            }
        }
    }
    std::cout << "분석 완료!" << std::endl;

    // 4. 최종 결과 출력
    // 사용자가 명시한 분모(467065)와 실제 계산된 분모가 일치하는지 확인
    if (beam_trigger_count != 467065) {
        std::cout << "\n경고: 계산된 Beam Trigger 이벤트 수(" << beam_trigger_count 
                  << ")가 예상값(467065)과 다릅니다. 결과를 주의해서 해석하세요." << std::endl;
    }

    double ratio = 0.0;
    if (beam_trigger_count > 0) {
        ratio = (double)htof_trigger_count / beam_trigger_count * 100.0;
    }

    std::cout << "\n--- HTOF 트리거 효율 분석 결과 ---" << std::endl;
    std::cout << "분석 기준 이벤트 수 (T0 && BH2): " << beam_trigger_count << "개" << std::endl;
    std::cout << "HTOF Multiplicity >= 2 조건 만족 이벤트 수: " << htof_trigger_count << "개" << std::endl;
    std::cout << std::fixed << std::setprecision(3); // 소수점 3자리까지
    std::cout << "효율 (비율): " << ratio << " %" << std::endl;
    std::cout << "------------------------------------" << std::endl;
}