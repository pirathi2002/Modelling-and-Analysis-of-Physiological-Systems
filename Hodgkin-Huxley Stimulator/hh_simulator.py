import sys
from PyQt5.QtWidgets import (QApplication, QMainWindow, QVBoxLayout, QWidget, 
                            QPushButton, QComboBox, QSpinBox, QLabel, 
                            QHBoxLayout, QFileDialog, QMessageBox, QTextEdit,
                            QStatusBar)
from PyQt5.QtCore import QTimer, QThread, pyqtSignal, Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import re
from scipy.integrate import odeint

class HodgkinHuxley:
    def __init__(self):
        self.C_m = 1.0
        self.g_Na = 120.0
        self.g_K = 36.0
        self.g_L = 0.3
        self.E_Na = 55.0
        self.E_K = -72.0
        self.E_L = -49.387
        self.custom_ions = {"Ca": {"g": 0.0, "E": 120.0, "co": 2.4, "ci": 0.0001, "z": 1, "gate": 0.0}}
        self.temp = 6.3
        self.Q10 = 3.0
        self.stimuli = [{"start": 5.0, "duration": 1.0, "amplitude": 7.0}]
        self.t_max = 100.0
        self.dt = 0.01
        self.V = -65.0
        self.m = 0.05
        self.h = 0.6
        self.n = 0.32
        self.results = {}

    def add_custom_ion(self, name, g, E, co, ci, z=1, gate=0.0):
        self.custom_ions[name] = {"g": g, "E": E, "co": co, "ci": ci, "z": z, "gate": gate}

    def add_stimulus(self, start, duration, amplitude):
        self.stimuli.append({"start": start, "duration": duration, "amplitude": amplitude})

    def clear_stimuli(self):
        self.stimuli = []

    def alpha_m(self, V):
        return 0.1 * (25 - V) / (np.exp((25 - V) / 10) - 1)

    def beta_m(self, V):
        return 4.0 * np.exp(-V / 18)

    def alpha_h(self, V):
        return 0.07 * np.exp(-V / 20)

    def beta_h(self, V):
        return 1.0 / (np.exp((30 - V) / 10) + 1)

    def alpha_n(self, V):
        return 0.01 * (10 - V) / (np.exp((10 - V) / 10) - 1)

    def beta_n(self, V):
        return 0.125 * np.exp(-V / 80)

    def temperature_factor(self):
        return self.Q10**((self.temp - 6.3) / 10)

    def stimulus_current(self, t):
        current = 0.0
        for stim in self.stimuli:
            if stim['start'] <= t <= stim['start'] + stim['duration']:
                current += stim['amplitude']
        return current

    def calculate_nernst(self, co, ci, z, temp):
        R = 8.314
        F = 96485
        T = 273 + temp
        return (R * T) / (z * F) * np.log(co / ci) * 1000

    def dALLdt(self, X, t):
        V, m, h, n = X
        k = self.temperature_factor()
        I_Na = self.g_Na * m**3 * h * (V - self.E_Na)
        I_K = self.g_K * n**4 * (V - self.E_K)
        I_L = self.g_L * (V - self.E_L)
        I_custom = sum(params['g'] * params['gate'] * (V - params['E']) for params in self.custom_ions.values())
        I_stim = self.stimulus_current(t)
        dVdt = (I_stim - I_Na - I_K - I_L - I_custom) / self.C_m
        dmdt = k * (self.alpha_m(V) * (1 - m) - self.beta_m(V) * m)
        dhdt = k * (self.alpha_h(V) * (1 - h) - self.beta_h(V) * h)
        dndt = k * (self.alpha_n(V) * (1 - n) - self.beta_n(V) * n)
        return dVdt, dmdt, dhdt, dndt

    def run_simulation(self):
        t = np.arange(0, self.t_max, self.dt)
        X0 = (self.V, self.m, self.h, self.n)
        solution = odeint(self.dALLdt, X0, t)
        V, m, h, n = solution[:, 0], solution[:, 1], solution[:, 2], solution[:, 3]
        I_Na = self.g_Na * m**3 * h * (V - self.E_Na)
        I_K = self.g_K * n**4 * (V - self.E_K)
        I_L = self.g_L * (V - self.E_L)
        I_custom = {ion: np.full_like(t, params['g'] * params['gate'] * (V - params['E'])[0]) for ion, params in self.custom_ions.items()}
        stim = np.array([self.stimulus_current(time) for time in t])
        g_Na = self.g_Na * m**3 * h
        g_K = self.g_K * n**4
        g_L = np.full_like(t, self.g_L)
        self.results = {'time': t, 'V': V, 'm': m, 'h': h, 'n': n, 'I_Na': I_Na, 'I_K': I_K, 'I_L': I_L, 'I_custom': I_custom, 'stim': stim, 'g_Na': g_Na, 'g_K': g_K, 'g_L': g_L}
        return self.results

class HHThread(QThread):
    update_signal = pyqtSignal(dict)

    def __init__(self, hh_model):
        super().__init__()
        self.hh_model = hh_model

    def run(self):
        results = self.hh_model.run_simulation()
        self.update_signal.emit(results)

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.hh_model = HodgkinHuxley()
        self.init_ui()
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.update_plot)
        self.thread = None

    def init_ui(self):
        self.setWindowTitle('Hodgkin-Huxley Simulator')
        self.setGeometry(100, 100, 1000, 800)

        # Central widget and layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout(central_widget)

        # Control panel
        control_widget = QWidget()
        control_layout = QHBoxLayout(control_widget)

        # Simulation Parameters
        sim_frame = QWidget()
        sim_layout = QVBoxLayout(sim_frame)
        sim_layout.addWidget(QLabel('Simulation Parameters'))
        self.tmax_spin = QSpinBox()
        self.tmax_spin.setValue(100)
        self.dt_spin = QSpinBox()
        self.dt_spin.setValue(1)
        self.temp_spin = QSpinBox()
        self.temp_spin.setValue(63)
        sim_layout.addWidget(QLabel('Duration (ms):'))
        sim_layout.addWidget(self.tmax_spin)
        sim_layout.addWidget(QLabel('Time step (ms):'))
        sim_layout.addWidget(self.dt_spin)
        sim_layout.addWidget(QLabel('Temperature (x10 °C):'))
        sim_layout.addWidget(self.temp_spin)
        run_button = QPushButton('Run')
        run_button.clicked.connect(self.run_simulation)
        reset_button = QPushButton('Reset')
        reset_button.clicked.connect(self.reset_simulation)
        sim_layout.addWidget(run_button)
        sim_layout.addWidget(reset_button)
        control_layout.addWidget(sim_frame)

        # Stimulation
        stim_frame = QWidget()
        stim_layout = QVBoxLayout(stim_frame)
        stim_layout.addWidget(QLabel('Stimulation'))
        self.amp_spin = QSpinBox()
        self.amp_spin.setValue(7)
        self.start_spin = QSpinBox()
        self.start_spin.setValue(5)
        self.dur_spin = QSpinBox()
        self.dur_spin.setValue(1)
        stim_layout.addWidget(QLabel('Amplitude (μA/cm²):'))
        stim_layout.addWidget(self.amp_spin)
        stim_layout.addWidget(QLabel('Start (ms):'))
        stim_layout.addWidget(self.start_spin)
        stim_layout.addWidget(QLabel('Duration (ms):'))
        stim_layout.addWidget(self.dur_spin)
        add_stim_button = QPushButton('Add Stimulus')
        add_stim_button.clicked.connect(self.add_stimulus)
        clear_stim_button = QPushButton('Clear')
        clear_stim_button.clicked.connect(self.clear_stimuli)
        stim_layout.addWidget(add_stim_button)
        stim_layout.addWidget(clear_stim_button)
        control_layout.addWidget(stim_frame)

        # Ion Channels (simplified)
        ion_frame = QWidget()
        ion_layout = QVBoxLayout(ion_frame)
        ion_layout.addWidget(QLabel('Ion Channels'))
        self.ion_combo = QComboBox()
        self.ion_combo.addItems(self.hh_model.custom_ions.keys())
        ion_layout.addWidget(QLabel('Custom Ion:'))
        ion_layout.addWidget(self.ion_combo)
        self.ion_g_spin = QSpinBox()
        ion_layout.addWidget(QLabel('Conductance:'))
        ion_layout.addWidget(self.ion_g_spin)
        update_ion_button = QPushButton('Update')
        update_ion_button.clicked.connect(self.update_ion)
        ion_layout.addWidget(update_ion_button)
        control_layout.addWidget(ion_frame)

        layout.addWidget(control_widget)

        # Plot area
        self.fig = Figure(figsize=(8, 6), dpi=100)
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.canvas, self)
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)

        # Console
        self.console = QTextEdit()
        self.console.setReadOnly(True)
        layout.addWidget(self.console)

        # Status bar
        self.statusBar = QStatusBar()
        self.setStatusBar(self.statusBar)

        self.show()

    def log(self, message):
        timestamp = datetime.now().strftime("%H:%M:%S")
        self.console.append(f"[{timestamp}] {message}")
        self.statusBar.showMessage(message)

    def run_simulation(self):
        self.hh_model.t_max = self.tmax_spin.value()
        self.hh_model.dt = self.dt_spin.value() / 100.0  # Convert to ms
        self.hh_model.temp = self.temp_spin.value() / 10.0  # Convert from x10
        self.thread = HHThread(self.hh_model)
        self.thread.update_signal.connect(self.plot_results)
        self.thread.start()

    def plot_results(self, results):
        self.fig.clear()
        axs = [self.fig.add_subplot(411), self.fig.add_subplot(412), self.fig.add_subplot(413), self.fig.add_subplot(414)]
        axs[0].plot(results['time'], results['V'], 'b-')
        axs[0].set_title('Membrane Potential')
        axs[0].set_ylabel('V (mV)')
        axs[0].grid(True)
        axs[1].plot(results['time'], results['m'], 'r-', label='m')
        axs[1].plot(results['time'], results['h'], 'g-', label='h')
        axs[1].plot(results['time'], results['n'], 'b-', label='n')
        axs[1].set_title('Gating Variables')
        axs[1].set_ylabel('Value')
        axs[1].legend()
        axs[1].grid(True)
        axs[2].plot(results['time'], results['I_Na'], 'r-', label='I_Na')
        axs[2].plot(results['time'], results['I_K'], 'b-', label='I_K')
        axs[2].plot(results['time'], results['I_L'], 'g-', label='I_L')
        axs[2].plot(results['time'], results['stim'], 'k--', label='Stimulus')
        axs[2].set_title('Ionic Currents')
        axs[2].set_ylabel('Current (μA/cm²)')
        axs[2].legend()
        axs[2].grid(True)
        axs[3].plot(results['time'], results['g_Na'], 'r-', label='g_Na')
        axs[3].plot(results['time'], results['g_K'], 'b-', label='g_K')
        axs[3].plot(results['time'], results['g_L'], 'g-', label='g_L')
        axs[3].set_title('Ionic Conductances')
        axs[3].set_xlabel('Time (ms)')
        axs[3].set_ylabel('Conductance (mS/cm²)')
        axs[3].legend()
        axs[3].grid(True)
        self.fig.tight_layout()
        self.canvas.draw()
        self.log("Simulation completed.")

    def reset_simulation(self):
        self.hh_model = HodgkinHuxley()
        self.tmax_spin.setValue(100)
        self.dt_spin.setValue(1)
        self.temp_spin.setValue(63)
        self.amp_spin.setValue(7)
        self.start_spin.setValue(5)
        self.dur_spin.setValue(1)
        self.ion_g_spin.setValue(0)
        self.ion_combo.clear()
        self.ion_combo.addItems(self.hh_model.custom_ions.keys())
        self.console.clear()
        self.fig.clear()
        self.canvas.draw()
        self.log("Simulation reset.")

    def add_stimulus(self):
        amp = self.amp_spin.value()
        start = self.start_spin.value()
        dur = self.dur_spin.value()
        self.hh_model.add_stimulus(start, dur, amp)
        self.log(f"Added stimulus: {amp} μA/cm² at {start} ms for {dur} ms")

    def clear_stimuli(self):
        self.hh_model.clear_stimuli()
        self.log("All stimuli cleared")

    def update_ion(self):
        ion = self.ion_combo.currentText()
        if ion:
            self.hh_model.custom_ions[ion]['g'] = self.ion_g_spin.value()
            self.log(f"Updated ion {ion}: g={self.ion_g_spin.value()} mS/cm²")

    def update_plot(self):
        # Placeholder for real-time updates (e.g., serial data)
        pass

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())