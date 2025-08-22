#!/usr/bin/env bash
set -euo pipefail

# 1) 원격 변화가 있는지(충돌 방지) 확인
git fetch origin

# 2) 작업 현황 보여주기
echo "== git status =="
git status

# 3) 커밋 메시지 입력(기본값: 'sync')
read -p "Commit message (default: sync): " msg
msg=${msg:-sync}

# 4) 스테이징(무시 목록은 .gitignore가 자동 제외)
git add -A

# 5) 커밋(변화 없으면 오류 대신 메시지)
if git commit -m "$msg"; then
  echo "Committed with message: $msg"
else
  echo "Nothing to commit."
fi

# 6) 푸시
git push origin main
