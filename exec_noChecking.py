import subprocess
from decimal import Decimal
import sys

#this scrip executes mc without checking statistics
def format_using_decimal(value):
    # Convert the float to a Decimal using string conversion to avoid precision issues
    decimal_value = Decimal(str(value))
    # Remove trailing zeros and ensure fixed-point notation
    formatted_value = decimal_value.quantize(Decimal(1)) if decimal_value == decimal_value.to_integral() else decimal_value.normalize()
    return str(formatted_value)


if (len(sys.argv)!=3):
    print("wrong number of arguments")
    exit()

T=float(sys.argv[1])

unitCellNum=int(sys.argv[2])
TStr=format_using_decimal(T)
#############################################
#launch mc, i.e., giving initial conditions

launchResult=subprocess.run(["python3", "launch_one_run.py", "./dataAllUnitCell"+str(unitCellNum)+"/row0/T"+TStr+"/run_T"+TStr+".mc.conf"])
print(launchResult.stdout)
if launchResult.returncode!=0:
    print("error in launch one run: "+str(launchResult.returncode))
#############################################


#############################################
#make mc
targetName="run_mc"
compileErrCode=10
cmake_process=subprocess.Popen(["cmake","."], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
while True:
    output = cmake_process.stdout.readline()
    if output == '' and cmake_process.poll() is not None:
        break
    if output:
        print(output.strip())
stdout, stderr = cmake_process.communicate()
if stdout:
    print(stdout.strip())
if stderr:
    print(stderr.strip())

make_process=subprocess.Popen(["make",targetName], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
while True:
    output = make_process.stdout.readline()
    if output == '' and make_process.poll() is not None:
        break
    if output:
        print(output.strip())
stdout, stderr = make_process.communicate()
if stdout:
    print(stdout.strip())
if stderr:
    print(stderr.strip())
#############################################

#############################################
#run executable
cppExecutable="./run_mc"
process = subprocess.Popen([cppExecutable, "./dataAllUnitCell"+str(unitCellNum)+"/row0/T"+TStr+"/cppIn.txt"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
while True:
    output = process.stdout.readline()
    if output == '' and process.poll() is not None:
        break
    if output:
        print(output.strip())
stdout, stderr = process.communicate()
if stdout:
    print(stdout.strip())
if stderr:
    print(stderr.strip())

#############################################
