import shutil as sh
import os

findir = "C:/Users/cuent/Dropbox (Empresa)/FPL Shared Folder/Daniel/Perple_x_noMn"

if os.path.exists(findir):
    sh.rmtree(findir)

sh.copytree("C:/Users/cuent/Desktop/Perple_x",findir)