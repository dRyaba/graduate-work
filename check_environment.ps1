# Скрипт проверки окружения для сборки проекта
Write-Host "========================================" -ForegroundColor Cyan
Write-Host "Проверка окружения для сборки проекта" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan
Write-Host ""

$tools = @{
    "CMake" = @("cmake", "--version")
    "C++ Compiler (GCC)" = @("g++", "--version")
    "C++ Compiler (Clang)" = @("clang++", "--version")
    "C++ Compiler (MSVC)" = @("cl", "")
}

$found = @{}
$allOk = $true

foreach ($tool in $tools.GetEnumerator()) {
    $name = $tool.Key
    $cmdName = $tool.Value[0]
    $versionArg = $tool.Value[1]
    
    $cmd = Get-Command $cmdName -ErrorAction SilentlyContinue
    if ($cmd) {
        Write-Host "[✓] $name найден: " -NoNewline -ForegroundColor Green
        Write-Host $cmd.Source -ForegroundColor White
        
        if ($versionArg -ne "") {
            try {
                $version = & $cmdName $versionArg 2>&1 | Select-Object -First 1
                Write-Host "    Версия: $version" -ForegroundColor Gray
            } catch {
                # Ignore version errors
            }
        }
        $found[$name] = $true
    } else {
        Write-Host "[✗] $name не найден" -ForegroundColor Red
        $found[$name] = $false
    }
}

Write-Host ""
Write-Host "========================================" -ForegroundColor Cyan
Write-Host "Результаты проверки" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan

# Проверяем наличие хотя бы одного компилятора
$hasCompiler = $found["C++ Compiler (GCC)"] -or $found["C++ Compiler (Clang)"] -or $found["C++ Compiler (MSVC)"]

if ($found["CMake"] -and $hasCompiler) {
    Write-Host ""
    Write-Host "✓ Все необходимые инструменты установлены!" -ForegroundColor Green
    Write-Host ""
    
    if ($found["C++ Compiler (GCC)"]) {
        Write-Host "Рекомендуемая команда сборки:" -ForegroundColor Yellow
        Write-Host "  cmake .. -DCMAKE_BUILD_TYPE=Debug -G `"MinGW Makefiles`"" -ForegroundColor White
    } elseif ($found["C++ Compiler (Clang)"]) {
        Write-Host "Рекомендуемая команда сборки:" -ForegroundColor Yellow
        Write-Host "  cmake .. -DCMAKE_BUILD_TYPE=Debug -G `"MinGW Makefiles`" -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++" -ForegroundColor White
    } elseif ($found["C++ Compiler (MSVC)"]) {
        Write-Host "Рекомендуемая команда сборки:" -ForegroundColor Yellow
        Write-Host "  cmake .. -DCMAKE_BUILD_TYPE=Debug" -ForegroundColor White
    }
    
    Write-Host ""
    Write-Host "Затем:" -ForegroundColor Yellow
    Write-Host "  cmake --build . -j8" -ForegroundColor White
} else {
    Write-Host ""
    Write-Host "✗ Не все инструменты установлены!" -ForegroundColor Red
    Write-Host ""
    
    if (-not $found["CMake"]) {
        Write-Host "Необходимо установить:" -ForegroundColor Yellow
        Write-Host "  1. CMake: https://cmake.org/download/" -ForegroundColor White
    }
    
    if (-not $hasCompiler) {
        Write-Host "Необходимо установить компилятор C++:" -ForegroundColor Yellow
        Write-Host "  - MinGW-w64 (рекомендуется): https://www.msys2.org/" -ForegroundColor White
        Write-Host "  - Или LLVM/Clang: https://github.com/llvm/llvm-project/releases" -ForegroundColor White
        Write-Host ""
        Write-Host "Подробные инструкции см. в SETUP.md" -ForegroundColor Cyan
    }
}

Write-Host ""
