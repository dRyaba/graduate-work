# Руководство по установке и настройке

## Автономная система сборки (рекомендуется)

Этот проект настроен для работы с автономными компиляторами, которые не требуют Visual Studio.

## Вариант 1: MinGW-w64 (рекомендуется для Windows)

### Шаг 1: Установка MSYS2 (включает MinGW-w64)

1. Скачайте MSYS2: https://www.msys2.org/
2. Установите MSYS2 (по умолчанию в `C:\msys64`)
3. Откройте **MSYS2 MSYS** (не MinGW64!)

4. Обновите пакеты:
```bash
pacman -Syu
# Если попросит закрыть терминал - закройте и откройте снова
pacman -Su
```

5. Установите инструменты разработки:
```bash
pacman -S --needed base-devel mingw-w64-x86_64-toolchain mingw-w64-x86_64-cmake mingw-w64-x86_64-gdb
```

6. Добавьте MinGW в PATH:
   - Откройте **Системные переменные среды** (Win+R → `sysdm.cpl` → Дополнительно → Переменные среды)
   - В **Системные переменные** найдите `Path` → **Изменить**
   - Добавьте: `C:\msys64\mingw64\bin`
   - Добавьте: `C:\msys64\usr\bin`
   - **ОК** и перезапустите терминал

### Шаг 2: Проверка установки

Откройте **новый** PowerShell и проверьте:

```powershell
gcc --version
g++ --version
cmake --version
```

Должны быть версии компиляторов и CMake.

### Шаг 3: Сборка проекта

```powershell
cd C:\Users\User\source\repos\graduate-work
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug -G "MinGW Makefiles"
cmake --build . -j8
```

---

## Вариант 2: LLVM/Clang (альтернатива)

### Установка LLVM

1. Скачайте LLVM: https://github.com/llvm/llvm-project/releases
   - Выберите **LLVM-XX.X.X-win64.exe** (последняя версия)
2. При установке выберите **"Add LLVM to system PATH"**
3. Перезапустите терминал

### Проверка

```powershell
clang --version
clang++ --version
```

### Сборка

```powershell
cd C:\Users\User\source\repos\graduate-work
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug -G "MinGW Makefiles" -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++
cmake --build . -j8
```

---

## Вариант 3: Visual Studio Build Tools (если нужен MSVC)

Если вы предпочитаете MSVC, но не хотите устанавливать полную Visual Studio:

1. Скачайте **Build Tools for Visual Studio**: https://visualstudio.microsoft.com/downloads/
2. При установке выберите **"Desktop development with C++"**
3. CMake автоматически найдет компилятор

---

## Автоматическая проверка окружения

Создайте скрипт `check_environment.ps1`:

```powershell
Write-Host "Проверка окружения для сборки..." -ForegroundColor Cyan

$tools = @{
    "CMake" = "cmake"
    "C++ Compiler" = "g++"
}

$allOk = $true
foreach ($tool in $tools.GetEnumerator()) {
    $cmd = Get-Command $tool.Value -ErrorAction SilentlyContinue
    if ($cmd) {
        Write-Host "[OK] $($tool.Key): $($cmd.Source)" -ForegroundColor Green
        & $tool.Value --version | Select-Object -First 1
    } else {
        Write-Host "[FAIL] $($tool.Key) не найден!" -ForegroundColor Red
        $allOk = $false
    }
}

if ($allOk) {
    Write-Host "`nВсе инструменты установлены! Можно собирать проект." -ForegroundColor Green
} else {
    Write-Host "`nУстановите недостающие инструменты (см. SETUP.md)" -ForegroundColor Yellow
}
```

Запустите: `.\check_environment.ps1`

---

## Рекомендации

**Для Windows рекомендуется MinGW-w64 через MSYS2**, потому что:
- ✅ Полностью автономный (не требует Visual Studio)
- ✅ Кроссплатформенный (тот же код работает на Linux/Mac)
- ✅ Легко устанавливается и обновляется
- ✅ Поддерживает все современные стандарты C++
- ✅ Хорошая производительность

**Избегайте зависимости от Visual Studio**, если:
- Хотите кроссплатформенность
- Работаете в команде с разными ОС
- Нужна автономная система сборки

---

## Быстрый старт (MinGW-w64)

После установки MSYS2 и добавления в PATH:

```powershell
# 1. Очистить старую сборку (если была)
Remove-Item -Recurse -Force build -ErrorAction SilentlyContinue

# 2. Создать новую сборку
mkdir build
cd build

# 3. Конфигурировать с MinGW
cmake .. -DCMAKE_BUILD_TYPE=Debug -G "MinGW Makefiles"

# 4. Собрать
cmake --build . -j8

# 5. Запустить
.\graph_reliability.exe --help
```

---

## Устранение проблем

### "cmake: command not found"
- Убедитесь, что CMake добавлен в PATH
- Или используйте полный путь: `& "C:\Program Files\CMake\bin\cmake.exe"`

### "No CMAKE_CXX_COMPILER could be found"
- Установите MinGW-w64 через MSYS2 (см. выше)
- Или установите LLVM/Clang
- Добавьте компилятор в PATH

### "CMakeCache.txt was created in different directory"
- Удалите `build` директорию полностью
- Создайте новую и переконфигурируйте

### "make: command not found" (при использовании MinGW)
- Используйте генератор: `-G "MinGW Makefiles"`
- Или используйте: `cmake --build .` вместо `make`
