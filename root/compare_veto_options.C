#include <TFile.h>
#include <TTree.h>
#include <TParticle.h>
#include <vector>
#include <iostream>
#include <string>
#include <iomanip> // 출력 포맷팅을 위해 추가

#pragma link C++ class vector<TParticle>+;

void compare_veto_options() {
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

    // 2. 카운터 변수 초기화
    long long total_events = tree->GetEntries();
    long long base_trigger_count = 0; // 분석 기준 (T0 x BH2 x HTOF-2)
    long long veto1_count = 0;        // 1안 Veto에 걸린 이벤트 수
    long long veto2_count = 0;        // 2안 Veto에 걸린 이벤트 수

    std::cout << "총 " << total_events << "개의 이벤트를 분석합니다..." << std::endl;

    // 3. 모든 이벤트를 순회하며 조건 확인
    for (long long i = 0; i < total_events; ++i) {
        tree->GetEntry(i);

        // 기본 트리거(T0 x BH2 x HTOF-2) 조건 확인
        if (t0_hits && t0_hits->size() > 0 &&
            bh2_hits && bh2_hits->size() > 0 &&
            htof_hits && htof_hits->size() >= 2) {
            
            base_trigger_count++;

            // 1안 Veto 조건 확인: BVH_U && BVH_D && BH2
            // (BH2는 이미 기본 트리거에서 확인되었으므로 BVH_U && BVH_D만 확인하면 됨)
            if (bvhu_hits && bvhu_hits->size() > 0 &&
                bvhd_hits && bvhd_hits->size() > 0) {
                veto1_count++;
            }

            // 2안 Veto 조건 확인: T0 && BVH_D && BH2
            // (T0, BH2는 이미 기본 트리거에서 확인되었으므로 BVH_D만 확인하면 됨)
            if (bvhd_hits && bvhd_hits->size() > 0) {
                veto2_count++;
            }
        }
    }
    std::cout << "분석 완료!" << std::endl;

    // 4. 최종 결과 테이블 출력
    std::cout << "\n--- Veto 방안 비교 분석 결과 ---" << std::endl;
    std::cout << "분석 기준 이벤트 수 (T0 x BH2 x HTOF-2): " << base_trigger_count << " 개" << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    
    // 1안 결과 계산 및 출력
    long long survived1 = base_trigger_count - veto1_count;
    double suppression1 = (base_trigger_count > 0) ? (double)veto1_count / base_trigger_count * 100.0 : 0;
    std::cout << "[ 1안: Veto on (BVH_U && BVH_D && BH2) ]" << std::endl;
    std::cout << "  - Veto된 이벤트 수: " << veto1_count << " 개" << std::endl;
    std::cout << "  - 최종 생존 이벤트 수: " << survived1 << " 개" << std::endl;
    std::cout << "  - 억제율: " << std::fixed << std::setprecision(3) << suppression1 << " %" << std::endl;
    
    std::cout << "--------------------------------------------------------" << std::endl;

    // 2안 결과 계산 및 출력
    long long survived2 = base_trigger_count - veto2_count;
    double suppression2 = (base_trigger_count > 0) ? (double)veto2_count / base_trigger_count * 100.0 : 0;
    std::cout << "[ 2안: Veto on (T0 && BVH_D && BH2) ]" << std::endl;
    std::cout << "  - Veto된 이벤트 수: " << veto2_count << " 개" << std::endl;
    std::cout << "  - 최종 생존 이벤트 수: " << survived2 << " 개" << std::endl;
    std::cout << "  - 억제율: " << std::fixed << std::setprecision(3) << suppression2 << " %" << std::endl;

    std::cout << "--------------------------------------------------------" << std::endl;
}