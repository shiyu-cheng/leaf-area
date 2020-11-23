import os
from shutil import copyfile

src_dir = "F:\\PATH2-TLS\\LAI_repo\\VS\\x64\\Release\\"

configs = [
	{"dir_name": "E:\\Elli-6typical-20191125\\TreeIdeal-test\\1\\",
	 "cmd": "67 90 23 -img_debug -cache_pathlen -n_divisions 4 -use_gap_p -refine -ptx \".\\PointsTREEelliTest-1-cs1.ptx\" -envelope \".\\PointsTREEelliTest-1-cs1 - Cloud_CONCAVE_alphaShape_0.5.obj\" -otxt -resolution_degree 0.03"}
    #,{"dir_name": "E:\\Elli-6typical-20191125\\TreeIdeal-6000-sphericalEllipsoid-random\\1\\",
	 #"cmd": "-zenith_ranges 48 90 42 -img_debug -cache_pathlen -n_divisions 4 -use_gap_p -refine -ptx \".\PointsTREEelli6-1-cs1.ptx\" -envelope \".\PointsTREEelli6-1-cs1 - Cloud_CONCAVE_alphaShape_0.5.obj\" -otxt -resolution_degree 0.02"}
	]

for config in configs:
	dir_name = config["dir_name"]
	copyfile(src_dir + "LAIpathTLS.exe", dir_name + "LAIpathTLS.exe")
	copyfile(src_dir + "FreeImage.dll", dir_name + "FreeImage.dll")
	os.chdir(dir_name )
	os.system("LAIpathTLS.exe" + " " + config["cmd"])
