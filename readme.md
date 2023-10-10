# Simplified Planar Inverse Kinematics Simulator 
This tool helps visualize and debug planar inverse kinematics for a simplified robot arm, and then provides a corrected function for a specified arm, based on the RR IK process. To run this code, do the following: 

1. Clone repository and setup python virtual environment 
```bash 
git clone git@github.com:njanne19/planar_arm_ik.git
cd planar_arm_ik
python -m venv env 
source ./env/bin/activate # On linux/macos 
./env/Scripts/activate # On Windows
```

2. Install required dependancies 
```bash 
pip install numpy matplotlib PyQt5
```

3. Run Simulator! 
```bash
python ik_interactive.py
```