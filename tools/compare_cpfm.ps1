# Compare CPFM: graduate-work vs Nesterov harness
# Run from project root or build/
# Graphs: 2_3x3_blocks (s=0,t=16,D=10), 3_blocks_sausage_3x3 (s=0,t=28,D=12), 3_blocks_sausage_4x4 (s=0,t=50,D=18)

$build = if (Test-Path "build\graph_reliability.exe") { "build" } else { "." }
$graphs = if (Test-Path "graphs_data\2_3x3_blocks_kao.txt") { "graphs_data" } else { "..\graphs_data" }

Write-Host "=== graduate-work (Method 4) ===" -ForegroundColor Cyan
Write-Host "2_3x3_blocks (s=0, t=16, D=10):"
& "$build\graph_reliability.exe" --run "$graphs\2_3x3_blocks_kao.txt" -1 -1 10 4 1

Write-Host "`n3_blocks_sausage_3x3 (s=0, t=28, D=12):"
& "$build\graph_reliability.exe" --run "$graphs\3_blocks_sausage_3x3_kao.txt" -1 -1 12 4 1

Write-Host "`n3_blocks_sausage_4x4 (s=0, t=50, D=18) [may take 10+ min]:"
& "$build\graph_reliability.exe" --run "$graphs\3_blocks_sausage_4x4_kao.txt" -1 -1 18 4 1

Write-Host "`n=== Nesterov harness ===" -ForegroundColor Cyan
Push-Location $build
& .\nesterov_harness.exe
Pop-Location
