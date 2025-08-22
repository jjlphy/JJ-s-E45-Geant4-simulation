#!/usr/bin/env bash
set -euo pipefail

# === go to repo root ===
ROOT="$(git rev-parse --show-toplevel 2>/dev/null || true)"
if [[ -z "${ROOT}" || ! -d "${ROOT}/.git" ]]; then
  echo "❌ 여기 Git 저장소가 아닙니다. 프로젝트 루트에서 실행하세요."
  exit 1
fi
cd "${ROOT}"

echo "📍 Repo root: ${ROOT}"
echo

# === 1) .gitignore 업데이트 (중복 없이 안전 추가) ===
echo "📝 .gitignore 업데이트 중..."
touch .gitignore

add_ignore() {
  local pattern="$1"
  grep -qxF "$pattern" .gitignore || echo "$pattern" >> .gitignore
}

# 대규모/빈번 충돌 유발 경로들
add_ignore "fieldmap/"
add_ignore ".build/"
add_ignore "bin/"
add_ignore "**/__pycache__/"
add_ignore "runmanager/module/__pycache__/"
add_ignore "*.pyc"
add_ignore "*.o"
add_ignore "*.obj"
add_ignore "*.pcm"
add_ignore "*.rootmap"
add_ignore "*.so"
add_ignore "*.dylib"
add_ignore "*.a"
add_ignore "*.d"
add_ignore "*.log"
add_ignore "*.tmp"
add_ignore ".DS_Store"
# CMake 흔적(루트에 빌드했던 경우 대비; 이미 .build/로도 대부분 커버됨)
add_ignore "CMakeFiles/"
add_ignore "CMakeCache.txt"
add_ignore "cmake_install.cmake"
add_ignore "install_manifest.txt"
add_ignore "module.modulemap"
add_ignore "Makefile"
add_ignore "cmake-build-*/"

git add .gitignore
echo "✅ .gitignore 갱신 완료"
echo

# === 2) 이미 추적 중인 파일을 '인덱스에서만' 제거 (로컬 파일은 유지) ===
echo "🧹 Git 인덱스에서 빌드/대용량 경로 untrack 처리(로컬 파일은 유지)..."

untrack_path() {
  local p="$1"
  git ls-files -z -- "$p" | xargs -0r git rm --cached -r --ignore-unmatch >/dev/null || true
}

# 디렉토리 단위
untrack_path "fieldmap"
untrack_path ".build"
untrack_path "bin"
untrack_path "**/__pycache__"
untrack_path "runmanager/module/__pycache__"

# 확장자 단위 (전체에서 골라내기)
git ls-files -z | \
  grep -zE '\.(o|obj|pcm|rootmap|so|dylib|a|d|pyc|log|tmp)$' | \
  xargs -0r git rm --cached --ignore-unmatch >/dev/null || true

# macOS 잡파일
git ls-files -z -- ".DS_Store" | xargs -0r git rm --cached --ignore-unmatch >/dev/null || true

# CMake 흔적(루트 빌드 잔재)
untrack_path "CMakeFiles"
git ls-files -z -- "CMakeCache.txt" "cmake_install.cmake" "install_manifest.txt" \
                   "module.modulemap" "Makefile" | \
  xargs -0r git rm --cached --ignore-unmatch >/dev/null || true

echo "✅ 인덱스 정리 완료 (로컬 파일은 그대로 남아있음)"
echo

# === 3) 요약 및 다음 단계 안내 ===
echo "📊 현재 상태 요약:"
git status --short

cat <<'MSG'

다음 명령으로 커밋/푸시를 마무리하세요:

  # 만약 'All conflicts fixed but you are still merging.' 메시지가 보이면:
  git commit              # 기본 merge 메시지로 머지 종료
  # 아니면 그냥 메시지 직접 주고 싶으면:
  # git commit -m "chore: ignore build & fieldmap; untrack generated files (keep local)"

  git push origin main

TIP:
- 이후부터는 위 경로들이 자동 무시됩니다.
- 로컬의 fieldmap/과 빌드 산출물은 보존됩니다(단지 Git에서만 제외).
MSG

