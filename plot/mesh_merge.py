
# できるだけ傘モジュールで呼ぶと、モジュール差異に強いです
import glob
import vtk

files = sorted(glob.glob("result_fluid_sups/online_10-4-4/rom_result_000001.vtk.*"))

# 1つ読んで型判定
probe_r = vtk.vtkDataSetReader()
probe_r.SetFileName(files[0])
probe_r.Update()
is_poly = probe_r.GetOutput().IsA("vtkPolyData")

# Appender を用意
appender = vtk.vtkAppendPolyData() if is_poly else vtk.vtkAppendFilter()

# まとめて投入
for fn in files:
    r = vtk.vtkDataSetReader()
    r.SetFileName(fn)
    r.Update()
    appender.AddInputData(r.GetOutput())
appender.Update()

# --- クリーニング（あるものを使う） ---
clean = None
if is_poly:
    clean = vtk.vtkCleanPolyData()
    clean.SetInputData(appender.GetOutput())
    clean.Update()
    output = clean.GetOutput()
else:
    # UnstructuredGrid 用：候補を順に試す
    CleanUG = None
    for cand in ("vtkCleanUnstructuredGrid", "vtkStaticCleanUnstructuredGrid", "vtkCleanGrid"):
        CleanUG = getattr(vtk, cand, None)
        if CleanUG is not None:
            break
    if CleanUG is not None:
        clean = CleanUG()
        clean.SetInputData(appender.GetOutput())
        # 近接許容誤差を調整したい場合は次を有効化（例：1e-12 など）
        # if hasattr(clean, "SetToleranceIsAbsolute"): clean.SetToleranceIsAbsolute(True)
        # if hasattr(clean, "SetAbsoluteTolerance"): clean.SetAbsoluteTolerance(0.0)
        clean.Update()
        output = clean.GetOutput()
    else:
        # クリーニングなしでそのまま出力（必要なら後で ParaView の "Clean to Grid" を使う）
        output = appender.GetOutput()

# --- 書き出し ---
if is_poly:
    w = vtk.vtkXMLPolyDataWriter()
    w.SetFileName("merged.vtp")
else:
    w = vtk.vtkXMLUnstructuredGridWriter()
    w.SetFileName("merged.vtu")

w.SetInputData(output)
w.Write()
print("Done.")
