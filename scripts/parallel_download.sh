#!/bin/bash
# Segmented parallel download of the 1.15GB Liver Cell Atlas zip (server supports Accept-Ranges).
set -u
URL="https://www.livercellatlas.org/data_files/toDownload/rawData_human.zip"
TOTAL=1149697366
N=16
DIR="$(cd "$(dirname "$0")/.." && pwd)/data/GSE192742"
PARTS="$DIR/parts"
OUT="$DIR/rawData_human.zip"
mkdir -p "$PARTS"
SEG=$(( (TOTAL + N - 1) / N ))

fetch_part() {
  local i=$1
  local start=$(( i * SEG ))
  local end=$(( start + SEG - 1 ))
  if [ $end -gt $(( TOTAL - 1 )) ]; then end=$(( TOTAL - 1 )); fi
  local expected=$(( end - start + 1 ))
  local f="$PARTS/part_$(printf '%02d' $i)"
  local tries=0
  while true; do
    local have=0
    [ -f "$f" ] && have=$(stat -f%z "$f" 2>/dev/null || echo 0)
    if [ "$have" = "$expected" ]; then return 0; fi
    tries=$(( tries + 1 ))
    if [ $tries -gt 8 ]; then echo "PART $i FAILED after 8 tries (have=$have want=$expected)"; return 1; fi
    curl -sL --max-time 1800 -r ${start}-${end} "$URL" -o "$f"
  done
}

echo "Starting 16-segment parallel download (seg=$SEG bytes each)..."
for i in $(seq 0 $(( N - 1 ))); do fetch_part $i & done
wait

# verify all parts then concatenate in order
echo "Verifying parts..."
ok=1
for i in $(seq 0 $(( N - 1 ))); do
  start=$(( i * SEG )); end=$(( start + SEG - 1 ))
  if [ $end -gt $(( TOTAL - 1 )) ]; then end=$(( TOTAL - 1 )); fi
  expected=$(( end - start + 1 ))
  f="$PARTS/part_$(printf '%02d' $i)"
  have=$(stat -f%z "$f" 2>/dev/null || echo 0)
  if [ "$have" != "$expected" ]; then echo "  part $i incomplete ($have/$expected)"; ok=0; fi
done
if [ $ok -ne 1 ]; then echo "ABORT: parts incomplete"; exit 1; fi

echo "Concatenating -> $OUT"
cat "$PARTS"/part_* > "$OUT"
final=$(stat -f%z "$OUT")
echo "final size: $final (expected $TOTAL)"
if [ "$final" = "$TOTAL" ]; then
  echo "DOWNLOAD_OK"
  rm -rf "$PARTS"
  echo "=== zip listing (head) ==="
  unzip -l "$OUT" | head -30
else
  echo "DOWNLOAD_SIZE_MISMATCH"; exit 1
fi
