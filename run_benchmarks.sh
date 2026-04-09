#!/bin/bash
cd /c/Users/User/Projects/CppProjects/graduate-work

EXE="./build/graph_reliability.exe"
TIMEOUT=300
OUTFILE="benchmark_results.csv"

echo "graph,s,t,d,method,reliability,time_sec,recursions" > $OUTFILE

run_test() {
  local label=$1 file=$2 s=$3 t=$4 d=$5 m=$6
  printf "%-20s m=%s d=%-3s ... " "[$label]" "$m" "$d"
  result=$(timeout $TIMEOUT $EXE --run $file $s $t $d $m 2>&1)
  R=$(echo "$result" | grep "Reliability:" | awk '{print $2}')
  T=$(echo "$result" | grep "Avg Time:" | awk '{print $3}')
  C=$(echo "$result" | grep "Avg Recs:" | awk '{print $3}')
  if [[ -z "$R" ]]; then
    echo "TIMEOUT/ERROR"
    echo "$label,$s,$t,$d,$m,TIMEOUT,>300s,-" >> $OUTFILE
  else
    echo "R=$R  t=${T}s  recs=$C"
    echo "$label,$s,$t,$d,$m,$R,$T,$C" >> $OUTFILE
  fi
}

# Сосисочные графы (d = d_min)
for m in 3 4 5; do
  run_test "2block-3x3"  graphs_data/2_3x3_blocks_kao.txt         -1 -1  8  $m
  run_test "3block-3x3"  graphs_data/3_blocks_sausage_3x3_kao.txt  -1 -1 10  $m
  run_test "3block-4x4"  graphs_data/3_blocks_sausage_4x4_kao.txt  -1 -1 13  $m
  run_test "4block-3x3"  graphs_data/4_blocks_sausage_3x3_kao.txt  -1 -1 14  $m
  run_test "5block-3x3"  graphs_data/5_blocks_sausage_3x3_kao.txt  -1 -1 18  $m
  run_test "6block-3x3"  graphs_data/6_blocks_sausage_3x3_kao.txt  -1 -1 22  $m
done

# Реальные сети
for m in 3 4 5; do
  run_test "Geant2004"   graphs_data/Geant2004_kao.txt                                            -1 -1  6  $m
  run_test "Geant2009"   graphs_data/Geant2009_kao.txt                                            -1 -1 12  $m
  run_test "IEEE-118"    graphs_data/IEEE-118-node_kao.txt                                         -1 -1 10  $m
  run_test "UPS-Russia"  graphs_data/UPS_of_Russia_composed_with_colored_cut_vertices_kao.txt      -1 -1 18  $m
done

echo ""
echo "Done. Results in $OUTFILE"
