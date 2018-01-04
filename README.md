# Localization of Kidnapped Vehicle
Project 3 of term 2 of Udacity self-driving car nanodegree

## Setup

This project involves the Term 2 Simulator which can be downloaded [here](https://github.com/udacity/self-driving-car-sim/releases).

Follow instructions [here](https://github.com/udacity/CarND-Kidnapped-Vehicle-Project) to set up Ubuntu and Mac environments. Also, if you are using
one of these environments then copy src/main\_mac\_linux.cpp to src/main.cpp . 

Instructions for 
setting up a native Windows environment and Visual Studio project:

* follow steps 1-4 of the instructions [here](https://github.com/fkeidel/CarND-Term2-ide-profile-VisualStudio/tree/master/VisualStudio)
* change install-windows.bat :
    - replace "EKF" with "Kidnapped-Vehicle"
    - delete row  "copy main.cpp ..\..\src" 
* run install-windows.bat 
* open solution in VisualStudio, remove from the solution all files from project 1, add all files from project 3
* run steps 6 and 7 from [here](https://github.com/fkeidel/CarND-Term2-ide-profile-VisualStudio/tree/master/VisualStudio)

## Files
* src - folder of source files.
* data - folder with the dataset of landmark coordinates.
* ide_profiles - folder with IDE profiles for Ecliple and XCode.

## Usage
1. Set up project environment for your operating system.
2. Compile the files to create the main executable of the project (the name of executable depends on the operating system).
3. Run Term 2 Simulator
4. In the simulator, press "Next" button twice to switch to Project 3, then press "Select" button.
5. Run the main executable of the project.
6. Press "Start" button in the simulator to start the localization process.

## Troubleshooting
Localization process succeds when simulator completes 2443 time steps in 100 seconds. In slow computers simulator will be able to complete less than 2443 steps and will give "out of time" error message. In this case please reduce the number of particles by changing line 36 in main.cpp. By default the code uses 100 particles.
