# 🧠 Modelling and Analysis of Physiological Systems

This repository includes simulations, visualizations, and tools developed for studying various physiological systems — such as **cardiac**, **respiratory**, **neural**, **dendritic tree models**, and **compartmental systems** — using **MATLAB**, **CircAdapt**, and **Python GUI-based applications**.

---

## 🔍 Click to Explore:

- [🫀 Analysis of Cardiac Physiology](#-analysis-of-cardiac-physiology)
- [🌬️ Respiratory Mechanics Simulation](#-simulation-of-respiratory-mechanics)
- [🌳 Dendritic Tree Modeling](#-branched-cylinder-dendritic-tree-approximation)
- [🔁 Compartmental System Modeling](#-compartmental-system-modeling)
- [⚡ Properties of Hodgkin-Huxley Equations (MATLAB)](#-properties-of-the-hodgkin-huxley-equations---matlab)
- [🖥️ Hodgkin-Huxley Stimulator GUI (Python)](#-hodgkin-huxley-stimulator-gui-python)
- [📸 Project Snapshots](#-project-snapshots)

---

## 🫀 Analysis of Cardiac Physiology

Cardiac physiology is simulated and analyzed using **MATLAB** and **CircAdapt**. Key aspects include:

- Ventricular pressure-volume loop simulation
- Cardiac cycle visualization
- Arterial-ventricular coupling
- Simulation of disease/physiological states

> 🔧 Tools: `MATLAB`, `CircAdapt` (open-source cardiovascular model)

---

## 🌬️ Simulation of Respiratory Mechanics

Simulates the respiratory system using lumped-parameter models and includes:

- Lung compliance/resistance models
- Volume-pressure dynamics
- Breathing cycle simulation

> 🔧 Tools: `MATLAB`

---

## 🌳 Branched Cylinder Dendritic Tree Approximation

Dendritic structures are modeled using compartmental methods with branched cylinders to analyze:

- Passive membrane properties
- Signal attenuation along dendrites

> 🔧 Tools: `MATLAB`

---

## 🔁 Compartmental System Modeling

This section includes models for **physiological transport** and **exchange dynamics** using compartmental systems. These are useful for:

- Tracer kinetics
- Drug distribution
- Diffusion across membranes
- Multi-compartment exchange modeling

Each compartment represents a well-mixed volume with flow rates and exchange parameters defined by differential equations.
 
> 🔧 Tools: `MATLAB`

✅ Designed for education and experimentation with system-level physiology.

---

## ⚡ Properties of the Hodgkin-Huxley Equations – MATLAB

This section includes **MATLAB simulations** demonstrating the behavior of the Hodgkin-Huxley (HH) model:

- Activation/inactivation curves
- Gating variable dynamics (`m`, `h`, `n`)
- Conductance vs voltage graphs
- Action potential response


---

## 🖥️ Hodgkin-Huxley Stimulator GUI – Python

A **Python-based GUI application** (in development) to simulate and stimulate HH model dynamics with user-defined parameters.

  👩‍🎓 Built for students and researchers to study HH properties without writing code
  
  🧑‍🏫 User-friendly interface designed for ease of learning and analysis

### Features:
- GUI built with `PySimpleGUI` (custom installer provided)
- Real-time membrane potential plots
- Adjustable stimulus amplitude/duration
- Plot export and parameter control


> 🧪 Status: **Under construction**  
> 🧾 Executable available in [Releases](../../releases) — just run `.exe` to use.

---


> ![HH MATLAB](images/hh_matlab_example.png)

> **Python HH Stimulator App GUI**  
> ![HH GUI](Assets/hh_gui_example.png)

> **CircAdapt Cardiac Simulation**  
> ![Cardiac CircAdapt](images/circadapt_example.png)



---

## 🛠️ Tech Stack

- `Python` (NumPy, Matplotlib, PySimpleGUI)
- `MATLAB`
- `CircAdapt` model

---

## 🚧 To-Do

- [ ] Add calcium ion dynamics
- [ ] Improve GUI design
- [ ] Add signal export to CSV



