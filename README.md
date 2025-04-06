# GPS-Denied UAV Localization Using Feature Matching in Matlab
## Overview
This repository provides a **MATLAB-based implementation** of a GPS-denied UAV localization algorithm using **KAZE feature matching**, **image correction**, and **inertial sensor simulation**. The system is designed for **Software-In-the-Loop (SIL)** simulation and enables UAVs to estimate their position without GPS, using only **aerial imagery** and **IMU (INS) data**.

- **Visual feature extraction and matching** using KAZE
- **Image registration and correction** against a reference aerial map
- **Sensor simulation** for IMU feedback (accelerometer + gyroscope)
- **Evaluation loop** that simulates asynchronous data flow between UAV and estimator

This MATLAB version is designed for fast prototyping and visualization of the perception pipeline.

## Requirements
- MATLAB R2019a or newer

## Getting Started
### 1. Clone the Repository
```
bash
git clone https://github.com/Howard1181/GPS-Denied-UAV-Localization-Using-Feature-Matching_Matlab
cd GPS-Denied-UAV-Localization-Using-Feature-Matching_Matlab
```
### 2. Open Matlab
### 3. Run the Simulation
