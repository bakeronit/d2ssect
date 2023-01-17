from setuptools import Extension,setup
import subprocess
import os, sys

def remove_suffix(input_string, suffix):
    if suffix and input_string.endswith(suffix):
        return input_string[:-len(suffix)]
    return input_string

# If specified these override any inferred values
JELLYFISH_LIBRARY_DIR = os.environ.get("JELLYFISH_LIBRARY_DIR", "/usr/local/lib")
JELLYFISH_INCLUDE_DIR = os.environ.get("JELLYFISH_INCLUDE_DIR", None)

jellyfish_base_path = None
jellyfish_major_version = None

if not JELLYFISH_INCLUDE_DIR:    
    # Attempt to infer Jellyfish library location
    #
    jellyfish_check_result = subprocess.run(["jellyfish","--version"],capture_output=True,text=True)
    if (jellyfish_check_result.returncode!=0):
        print("Could not determine jellyfish version. Use environment variables JELLYFISH_LIBRARY_DIR and JELLYFISH_INCLUDE_DIR to specify paths to Jellyfish 2 libraries on your system")
    else:
        jellyfish_version_string=jellyfish_check_result.stdout.split()[1]
        jellyfish_major_version = int(jellyfish_version_string.split(".")[0])
        if jellyfish_major_version<2:
            print("Detected jellyfish version 1 on path. Unable to infer location of jellyfish libraries. Use environment variables JELLYFISH_LIBRARY_DIR and JELLYFISH_INCLUDE_DIR to specify paths to Jellyfish 2 libraries on your system")
        else:
            which_jellyfish_result = subprocess.run(["which","jellyfish"],capture_output=True,text=True)
            if (which_jellyfish_result.returncode==0):
                jellyfish_path = which_jellyfish_result.stdout.strip()
                jellyfish_base_path = remove_suffix(jellyfish_path,"bin/jellyfish")

            JELLYFISH_INCLUDE_DIR = jellyfish_base_path+"include/jellyfish"+"-"+jellyfish_version_string
            JELLYFISH_LIBRARY_DIR = jellyfish_base_path+"lib"   


sys.stdout.write("INCLUDE:"+JELLYFISH_INCLUDE_DIR+"\n")
sys.stdout.write("LIB:"+JELLYFISH_LIBRARY_DIR+"\n")

setup(
    name='d2ssect',
    version='0.0.8',
	ext_modules=[
        Extension(
            name="d2ssect.jellyfish",
            sources=["d2ssect/jellyfish.cpp"], # all sources are compiled into a single binary file
            include_dirs  = [JELLYFISH_INCLUDE_DIR],
            library_dirs = [JELLYFISH_LIBRARY_DIR],
            libraries = ['jellyfish-2.0'],
            language='c++',
        ),
    ]
)